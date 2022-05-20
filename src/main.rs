mod myfastq;
mod sqlite;
mod hset;
mod cb_umi_errors;
// use fastq::{parse_path, Record, Parser};
// use std::env::args;
use std::collections::{HashMap};
use counter::Counter;
use bktree::{BkTree, levenshtein_distance};
// use polars::frame::DataFrame;
// use polars::series::Series;
// use polars::prelude::NamedFrom;
use std::fs::File;
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};


fn top_n(counter: &Counter<(String, String), i32>, n: i32) -> Vec<(String, String)>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements
    let mut bk: BkTree<String> = BkTree::new(levenshtein_distance);

    let mut c = 0;
    for ((cb, umi), _freq) in counter.most_common(){

        let cb_umi = format!("{cb}_{umi}");
        if bk.find(cb_umi, 2).len() == 0{
            let cb_umi = format!("{cb}_{umi}");
            bk.insert(cb_umi);
            c += 1;
        }
        if c > n{
            break
        }
    }

    let mut top2: Vec<(String,String)> = Vec::new();

    for element in bk.into_iter(){
        let mut cb_umi = element.split('_');
        let cb = cb_umi.next().unwrap();
        let umi = cb_umi.next().unwrap();
        top2.push((String::from(cb), String::from(umi)));
    }
    top2
}


fn find_shadows(cb_umi: (String, String), filtered_map: &Counter<(String, String), i32>) -> HashMap<usize, i32>{

    let (cb, umi) = cb_umi;

    let mut all_muts: Vec<(String, String, usize)>= Vec::new();  // cb, umi, pos
    for pos in 0..umi.len() {
        for m in get_1bp_mutations(&umi, pos){
            let c2 = cb.clone();
            all_muts.push((c2, m, pos));
        }
    }
    // println!("{:?}", all_muts);

    let mut n_shadows_per_pos: HashMap<usize, i32> = HashMap::new();
    for (cb, umi, pos) in all_muts{
        if let Some(f) = filtered_map.get(&(cb, umi)){
            // for that shadow cb/umi, how often did we acutally see it
            n_shadows_per_pos.insert(pos, *f);
        }
        else{
            n_shadows_per_pos.insert(pos, 0);
        }
    }
    n_shadows_per_pos
}


fn get_1bp_mutations(seq: &String, pos: usize) -> Vec<String>{
    // return the three 1BP mutations of a sequence at the given position"
    let mut muts = Vec::with_capacity(3);

    for base in ['A', 'C', 'G', 'T']{
        let prefix = &seq[0..pos];
        let postfix = &seq[(pos+1)..];
        let m = format!("{prefix}{base}{postfix}");
        if m != *seq{
            muts.push(m);
        }
    }
    return muts
}

fn main() {
    // myfastq::run();
    // sqlite::run();
    // hset::run();

    // let fastq_file: String = "/home/michi/r1.fastq.gz".into();
    // let fastq_file: String = "/home/michi/Virus_barcode_CKDL220009934-1a_HN2MJDSX3_L3_1.fq.gz".into();
    let fastq_file: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();

    // let whitelist_file: String = "/home/michi/3M-february-2018.txt.gz".into();
    let whitelist_file: String = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz".into();
    let whitelist = cb_umi_errors::parse_whitelist_gz(whitelist_file);
    println!("Whitelist len {}", whitelist.len());

    
    let countmap = cb_umi_errors::count_cb_umi(fastq_file);
    // transform into shadow counter
    println!("len of counter {}", countmap.len());

    // filter for whitelist only entries
    // turn it into a counter
    println!("Filtering for whilelist");
    let mut filtered_map: Counter<(String, String), i32> = Counter::new();
    // let mut filtered_map: HashMap<(String, String), i32> = HashMap::new();

    for ((cb, umi), count) in countmap.iter(){
        if whitelist.contains(cb){
            filtered_map.insert((cb.clone(), umi.clone()), count.clone());
        }
    }
    println!("len of filtered counter {}", filtered_map.len());
    // println!("most common counter {:?}", filtered_map.most_common());

    // countmap.drain_filter(|k, v| {
    //     k.0 == k.1
    // });

    // now the hard part: group by CB, look at all UMIs therein
    println!("calculating most common");
    let most_common: Vec<(String,String)> = top_n(&filtered_map, 1000);
    println!("most common {:?}", most_common.len());

    // for each entry in most common, find potential shadows
    // add them to polars
    let mut polars_data: HashMap<String, Vec<i32>> = HashMap::new();

    for mc in most_common{
        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &filtered_map);
        for (position, n_shadows) in nshadows_per_position.iter(){
            if *n_shadows > 0 {
                println!("shadows {:?}", nshadows_per_position);
                // panic!("FOUND")

            }
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = filtered_map.get(&mc2).unwrap();
        polars_data.entry("total".into()).or_insert(vec![]).push(*total_counts)

    }
    
    // to polars dataframe
    let mut df = DataFrame::new(
        polars_data.into_iter()
            .map(|(name, values)| Series::new(&format!("{name}"), values))
            .collect::<Vec<_>>()).unwrap();
    
    println!("{:?}", df);

    // write to CSV
    let mut output_file: File = File::create("/tmp/out.csv").unwrap();
    CsvWriter::new(&mut output_file)
        .has_header(true)
        .finish(&mut df)
        .unwrap();    
}
