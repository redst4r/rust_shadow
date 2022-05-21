use std::collections::HashSet;

// mod myfastq;
// mod sqlite;
// mod hset;
mod cb_umi_errors;
mod cb_umi_sketch;

mod cb_errors;
mod utils;
mod sketching;
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
 

    let whitelist_file: String = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz".into();
    let fastq_file1: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();
    let fastq_list = vec![fastq_file1];


    // let whitelist_file: String = "/home/mstrasse/TB4/resources/3M-february-2018.txt.gz".into();

    // let fastq_file1: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_001.fastq.gz".into();
    // let fastq_file2: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_002.fastq.gz".into();
    // let fastq_file3: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_003.fastq.gz".into();
    // let fastq_file4: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_004.fastq.gz".into();
    // let fastq_file5: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_005.fastq.gz".into();
    // let fastq_file6: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_006.fastq.gz".into();
    // let fastq_file7: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_007.fastq.gz".into();
    // let fastq_file8: String = "/home/mstrasse/TB4/tmp/fastq_tmp/L05Ai/bamtofastq_S1_L000_R1_008.fastq.gz".into();
    // let fastq_list = vec![fastq_file1, fastq_file2, fastq_file3, fastq_file4, fastq_file5, fastq_file6, fastq_file7, fastq_file8] ;

    // let fastq_file1: String = "/home/michi/01_Day2_GE_S1_L001_R1_001.fastq.gz".into();
    // let fastq_list = vec![fastq_file1];


    // sketching::run(fastq_list);
    // sketching::run_top(fastq_list);
    // sketching::run_gt1(fastq_list);
    // cb_umi_sketch::run_GB1(fastq_list, whitelist_file, "/tmp/out.csv".to_string());


    // use indicatif::ProgressBar;
    // let bar = ProgressBar::new_spinner();
    // for i in 0..100{
    //     bar.inc(1);
    // }
    // bar.finish();

    // cb_umi_sketch::run_topN(&fastq_list, whitelist_file, "/tmp/out.csv".to_string());
    cb_umi_errors::run(&fastq_list, whitelist_file, "/tmp/full_out.csv".to_string());


    if false{
    // cb_errors::run(fastq_list, whitelist_file, "/tmp/out_cb.csv".to_string())
    // cb_umi_errors::run(fastq_list, whitelist_file, "/tmp/out_cb_umi.csv".to_string());


    // cb_errors::run(fastq_file, whitelist_file, "/tmp/out_cb.csv".to_string())

    }
}



