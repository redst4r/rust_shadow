mod myfastq;
mod sqlite;
mod hset;
mod cb_umi_errors;
// use fastq::{parse_path, Record, Parser};
// use std::env::args;
use std::collections::{HashMap};
use std::fs::File;
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};




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

    
    let mut countmap = cb_umi_errors::count_cb_umi(fastq_file);
    // transform into shadow counter
    println!("len of counter {}", countmap.len());

    // filter for whitelist only entries
    println!("Filtering for whilelist");

    let mut keys_to_remove: Vec<(String, String)> = Vec::new();
    for ((cb, umi), _count) in countmap.iter(){
        if !whitelist.contains(cb){
            keys_to_remove.push((cb.clone(), umi.clone()));
        }
    }

    for k in keys_to_remove{
        countmap.remove(&k);
    }

    println!("len of filtered counter {}", countmap.len());
    // unstable:
    // countmap.drain_filter(|k, v| {
    //     k.0 == k.1
    // });

    // now the hard part: group by CB, look at all UMIs therein
    println!("calculating most common");
    let most_common: Vec<(String,String)> = cb_umi_errors::top_n(&countmap, 1000);
    println!("most common {:?}", most_common.len());

    // for each entry in most common, find potential shadows
    // add them to polars
    let mut polars_data: HashMap<String, Vec<i32>> = HashMap::new();

    for mc in most_common{
        let mc2 = mc.clone();
        let nshadows_per_position = cb_umi_errors::find_shadows(mc, &countmap);
        for (position, n_shadows) in nshadows_per_position.iter(){
            if *n_shadows > 0 {
                println!("shadows {:?}", nshadows_per_position);
                // panic!("FOUND")

            }
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = countmap.get(&mc2).unwrap();
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
