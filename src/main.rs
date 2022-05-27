// mod myfastq;
// mod sqlite;
// mod hset;
mod cb_umi_errors;
mod cb_umi_sketch;

mod cb_errors;
mod utils;
mod sketching;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {

    // #[clap(parse(from_os_str), short ='o', long = "output")] 
    // output: std::path::PathBuf,    
    /// Path to output file
    #[clap(short ='o', long = "output")] 
    output: String,    

    // #[clap(parse(from_os_str), short = 'o', long= "whitelist")] 
    // whitelist: std::path::PathBuf,
    /// 10x CB whitelist
    #[clap(short = 'w', long= "whitelist")] 
    whitelist: String,

    #[clap(short = 'n', long= "ntop")] 
    topn: usize,

    /// List of fastq files
    #[clap()]
    fastq_list: Vec<String>,

    #[clap(short = 'c', long= "command")] 
    command: String,
    
}

fn main() {
    // myfastq::run();
    // sqlite::run();
    // hset::run();
    
    let args = Cli::parse();

    println!("Whitelist {:?}",args.whitelist);
    println!("Output {:?}",args.output);
    println!("Top N {:?}",args.topn);
    println!("FASTQ {:?}",args.fastq_list);

    if args.command == "cb"{
        println!("Doing CB only");
        cb_errors::run(&args.fastq_list, args.whitelist, args.output, args.topn)
    }
    else if args.command == "cb_umi_sketch" {
        println!("Doing CB_UMI sketch");
        cb_umi_sketch::run_topN(&args.fastq_list, args.whitelist, args.output, args.topn)
    }
    else if args.command == "cb_umi_exact" {
        println!("Doing CB_UMI sketch");
        println!("WARNING: MEMORY INTENSIVE!!");
        cb_umi_errors::run(&args.fastq_list, args.whitelist, args.output, args.topn)
    }    
    else{
        panic!("unknown command")
    }
    // cb_umi_sketch::run_topN(&args.fastq_list, args.whitelist, args.output, topn);
    // cb_umi_errors::run(&args.fastq_list, args.whitelist, "/tmp/full_out.csv".to_string(), topn);



    // let args: Vec<String> = env::args().collect();
    // let fastq_list = &args[1..];

    // let fastq_file: String = "/home/michi/r1.fastq.gz".into();
    // let fastq_file1: String = "/home/michi/01_Day2_GE_S1_L001_R1_001.fastq.gz".into();
    // let fastq_file: String = "/home/michi/Virus_barcode_CKDL220009934-1a_HN2MJDSX3_L3_1.fq.gz".into();
    // let whitelist_file: String = "/home/michi/3M-february-2018.txt.gz".into();
    // let fastq_list = vec![fastq_file1];

    // let whitelist_file: String = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz".into();
    // let fastq_file1: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();
    // let fastq_list = vec![fastq_file1];

    // /home/michi/.cargo/bin/cargo run --release -- -w /home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt.gz --ntop 10000 --output /tmp/cb_only.csv --command cb /home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz



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
    // cb_umi_errors::run(&fastq_list, whitelist_file, "/tmp/full_out.csv".to_string());


    if false{
    // cb_umi_errors::run(fastq_list, whitelist_file, "/tmp/out_cb_umi.csv".to_string());


    // cb_errors::run(fastq_file, whitelist_file, "/tmp/out_cb.csv".to_string())

    }
}



