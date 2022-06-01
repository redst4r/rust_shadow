use std::collections::HashMap;
use bktree::BkTree;
use crate::bus::{CellIterator, BusRecord};
use crate::utils::{int_to_seq, CbUmi, write_to_csv, my_hamming, seq_to_int};
use counter::Counter;
use crate::cb_umi_errors::find_shadows;
use polars::prelude::*;

use indicatif::{ProgressBar, ProgressStyle};

#[cfg(test)]
#[test]
fn main(){
    run(&"/home/michi/output.corrected.sort.bus".to_string(), &"/tmp/cb.csv".to_string(), 1000)
}


pub fn run(busfile: &String, outfile: &String, nmax: usize, aggregate: bool){
    // nmax: maximum number of barcodes to consider, should be on the order of several millions
    let cb_iter = CellIterator::new(&busfile);

    let mut df = DataFrame::default();
    let bar = ProgressBar::new(nmax as u64);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise} ETA {eta}] {bar:40.cyan/blue} {pos}/{len} {per_sec}")
        .progress_chars("##-"));

    let mut number_of_umis_seen = 0;
    for records in cb_iter
        .map(|(_cb, rec)| rec)
        .filter(|rec| rec.len()>100){

        let df_single_cell = do_single_cb(records);

        bar.inc(df_single_cell.height() as u64);
        number_of_umis_seen += df_single_cell.height();

        if aggregate{
            // aggregate the counts of errors across UMIs
            let df_single_cell = df_single_cell.groupby(["CB"]).unwrap().sum().unwrap();
        }
        if df.is_empty(){
            df = df_single_cell;
        }
        else{
            df = df.vstack(&df_single_cell).unwrap();
        }

        if number_of_umis_seen > nmax{
            bar.finish();
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

    // sorting the records?! might make BKtrees faster
    // records.sort_by(|a, b| b.UMI.cmp(&a.UMI));

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

    // TODO pretty stupid to do it like this, i.e iterating over all potential shadows
    // we could just build a BKTree out of all records and easily query for related seqs
    // onnly problem: we'd have to figure out the position of the mutation

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

    for (cbumi, _freq) in counter.most_common(){

        let umi = cbumi.umi.clone();
        let matches = bk.find(umi.clone(), 1);
        if matches.len() == 0{
            correct_umis.push(cbumi.clone());  // add it to the topN list
            bk.insert(umi);
        }
        else{
            bk.insert(umi);
        }     
    }
    correct_umis
}



//---------------------------------------------------------------------------------
pub fn find_correct_umis_faster(counter: &Counter<CbUmi, u32>) -> (Vec<CbUmi>, BkTree<String>){

    // TODO we can probably speed this up. All the elements in Counter have the SAME CB!! not useful in the BKtree and just slows it down
    // starting with the most frequent UMIs, blacklist all UMIs 1BP away (by putting them into the BKtree)
    // keep itearting the UMIs util were done with all
    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    let mut correct_umis: Vec<CbUmi> = Vec::new();

    for (cbumi, _freq) in counter.most_common(){

        let umi = cbumi.umi.clone();
        let matches = bk.find(umi.clone(), 1);
        if matches.len() == 0{
            correct_umis.push(cbumi.clone());  // add it to the topN list
            bk.insert(umi);
        }
        else{
            bk.insert(umi);
        }     
    }
    (correct_umis, bk)
}

fn do_single_cb_bktree(bus_records: Vec<BusRecord>) -> DataFrame{

    // sorting the records?! might make BKtrees faster
    // records.sort_by(|a, b| b.UMI.cmp(&a.UMI));

    // count UMI frequencies

    let mut freq_map: Counter<CbUmi, u32> = Counter::new();
    for cbumi in bus_records.iter().map(|r| bus_record_to_cbumi(r) ){
        let c = freq_map.entry(cbumi).or_insert(0);
        *c +=1;
    }
    let (correct_umis, bk) = find_correct_umis_faster(&freq_map); 

    // find the shadows 
    let mut polars_data: HashMap<String, Vec<u32>> = HashMap::new();
    let mut cellnames: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();


    for mc in correct_umis{
        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows_bktree(mc, &bk, &freq_map);
        // let nshadows_per_position = find_shadows(mc, &freq_map);

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


fn find_first_mismatch(a: &String, b: &String) -> Option<usize>{
    let mut counter = 0;
    for (c1, c2) in  std::iter::zip((*a).chars(), (*b).chars()){
        if c1!=c2{
            return Some(counter);
        }
        else{
            counter+=1;
        }
    }
    return None
}

pub fn find_shadows_bktree(correct_cbumi: CbUmi, bk: &BkTree<String>, filtered_map: &Counter<CbUmi, u32>) -> HashMap<usize, u32>{
    // for a given "true" CB/UMI find the number of shadowed reads (i.e. reads with a single subsitution)
    // and count their number (position specific)
    // we get a dcitionary with position -> #shadow reads
    //
    
    let query_umi = correct_cbumi.umi.clone();
    let mut n_shadows_per_pos: HashMap<usize, u32> = HashMap::new();
    // lets just initialize with 0, to make sure all positions are in the table
    for i in 0..28{
        n_shadows_per_pos.insert(i, 0);
    }

    let matches = bk.find(query_umi, 1); // this should NOT include the query itself
    for (umi, _d) in matches{
        if _d == 0{
            continue;
        }
        let mutated_cbumi = CbUmi{cb: correct_cbumi.cb.clone(), umi: (*umi).clone()};
        let pos = find_first_mismatch(&correct_cbumi.umi, umi).unwrap();

        if let Some(f) = filtered_map.get(&mutated_cbumi){ // if we have seen this sequence. THIS MUST BE TRUE!!
            let current_counter = n_shadows_per_pos.entry(pos).or_insert(0);
            *current_counter += f;
        }
        else {
            panic!("Shouldnt happen. If the seq is in the BKtree, we must know its freq too")
        }
    };
    n_shadows_per_pos
}