use std::collections::HashMap;
use bktree::BkTree;
use crate::bus::{CellIterator, BusRecord};
use crate::utils::{int_to_seq, CbUmi, write_to_csv, my_hamming};
use counter::Counter;
use crate::cb_umi_errors::{find_shadows};
use polars::prelude::*;

// use indicatif::ProgressBar;

#[cfg(test)]
#[test]
fn main(){
    run(&"/home/michi/output.corrected.sort.bus".to_string(), &"/tmp/cb.csv".to_string(), 100_000)
}


pub fn run(busfile: &String, outfile: &String, nmax: usize){
    // nmax: maximum number of barcodes to consider, should be on the order of several millions
    let cb_iter = CellIterator::new(&busfile);

    let mut df = DataFrame::default();
    // let bar = ProgressBar::new();

    for (i, records) in cb_iter.map(|(_cb, rec)| rec).enumerate(){
        // println!("Doing cell with #{} records", records.len());
        let df_single_cell = do_single_cb(records);
        if df.is_empty(){
            df = df_single_cell;
        }
        else{
            df = df.vstack(&df_single_cell).unwrap();
        }
        if i % 100 == 0{
            println!("Iteration {}", i);
            println!("current height {}", df.height());

        }
        if df.height() > nmax{
            break
        }
    }
    println!("final height{}", df.height());

    write_to_csv(&mut df, outfile.to_string());

    // let fh = File::create("example.parquet").expect("could not create file");
    // ParquetWriter::new(fh).finish(&mut df);

    write_to_csv(&mut df.sum(), "/tmp/cb_sum.csv".to_string());


}

fn bus_record_to_cbumi(r: &BusRecord) -> CbUmi{
    CbUmi {cb: int_to_seq(r.CB, 16), umi: int_to_seq(r.UMI, 12)}
}

fn do_single_cb(bus_records: Vec<BusRecord>) -> DataFrame{

    // count UMI frequencies
    let mut freq_map: Counter<CbUmi, u32> = Counter::new();

    for cbumi in bus_records.iter().map(|r| bus_record_to_cbumi(r) ){
        let c = freq_map.entry(cbumi).or_insert(0);
        *c +=1;
    }

    let correct_umis = find_correct_umis(&freq_map); 

    // find the shadows 
    let mut polars_data: HashMap<String, Vec<u32>> = HashMap::new();
    let mut cellnames: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();

    for mc in correct_umis{
        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &freq_map);

        for (position, n_shadows) in nshadows_per_position.iter(){
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = freq_map.get(&mc2).unwrap();
        polars_data.entry("n_real".into()).or_insert(vec![]).push(*total_counts);

        cellnames.push(mc2.cb);
        umis.push(mc2.umi);
    }

    // to polars dataframe
    // problem, we need the colums sorted, for later stacking
    let mut vec_series = 
        polars_data.into_iter()
            .map(|(name, values)| Series::new(&format!("{name}"), values))
            .collect::<Vec<_>>();
    vec_series.sort_by(|a, b| a.name().partial_cmp(b.name()).unwrap());
    
    let df = DataFrame::new(vec_series).unwrap();

    let df_cb = Series::new("CB", cellnames);
    let df_umi = Series::new("UMI", umis);

    let df_final = df.hstack(&[df_cb, df_umi]).unwrap();    
    df_final

}




pub fn find_correct_umis(counter: &Counter<CbUmi, u32>) -> Vec<CbUmi>{

    // TODO we can probably speed this up. All the elements in Counter have the SAME CB!! not useful in the BKtree and just slows it down
    // starting with the most frequent UMIs, blacklist all UMIs 1BP away (by putting them into the BKtree)
    // keep itearting the UMIs util were done with all
    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    let mut correct_umis: Vec<CbUmi> = Vec::new();

    let mut c = 0;
    // let bar = ProgressBar::new(counter.len() as u64);

    for (cbumi, _freq) in counter.most_common(){
        // bar.inc(1);

        let umi = cbumi.umi.clone();
        let matches = bk.find(umi.clone(), 1);
        if matches.len() == 0{
            correct_umis.push(cbumi.clone());  // add it to the topN list
            bk.insert(umi);
            c += 1;
        }
        else{
            bk.insert(umi);
        }

        // if c >= n{
            // break
        // }
        // if c % 1_000 == 0{
        //     println!("BKTree Iteration {c}");
        // }          
    }
    // bar.finish();

    correct_umis
}