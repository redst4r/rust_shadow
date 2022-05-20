// mod myfastq;
// mod sqlite;
// mod hset;
mod cb_umi_errors;
mod cb_errors;
// use fastq::{parse_path, Record, Parser};
// use std::env::args;


fn main() {
    // myfastq::run();
    // sqlite::run();
    // hset::run();

    // let fastq_file: String = "/home/michi/r1.fastq.gz".into();
    // let fastq_file: String = "/home/michi/Virus_barcode_CKDL220009934-1a_HN2MJDSX3_L3_1.fq.gz".into();
    let fastq_file: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();

    // let whitelist_file: String = "/home/michi/3M-february-2018.txt.gz".into();
    let whitelist_file: String = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz".into();
 
    cb_umi_errors::run(fastq_file, whitelist_file, "/tmp/out_cb_umi.csv".to_string());


    let fastq_file: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();
    let whitelist_file: String = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz".into();
    
    cb_errors::run(fastq_file, whitelist_file, "/tmp/out_cb.csv".to_string())


}
