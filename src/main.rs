// mod myfastq;
// mod sqlite;
// mod hset;
mod cb_umi_errors;
mod cb_umi_sketch;
mod cb_errors;
mod utils;
// mod sketching;

use clap::{self, Parser, Subcommand, Args};
mod bus;
mod bus_multi;
// mod merger;
mod busmerger;
mod tso_error;
mod cb_umi_per_cell;
mod cb_umi_per_cell_gene;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {  
    /// Path to output file
    #[clap(short ='o', long = "output")] 
    output: String,    

    #[clap(subcommand)]
    command: MyCommand
}

#[allow(non_camel_case_types)]
#[derive(Subcommand)]
enum MyCommand {
    cb(FastqArgs),
    cb_umi_exact(FastqArgs),
    cb_umi_sketch(FastqArgs),
    cb_umi_cell(BusArgs),
    cb_umi_cell_aggr(BusArgs),
    cb_umi_cell_gene(BusArgs2),
    busmerge(BusMergeArgs),
    tso_error(TSOArgs)
}

#[derive(Args)]
struct FastqArgs{
    /// 10x CB whitelist
    #[clap(short = 'w', long= "whitelist")] 
    whitelist: String,

    /// topN CM/UMI to consider
    #[clap(short = 'n', long= "ntop")] 
    topn: usize,

    /// List of fastq files
    #[clap()]
    fastq_list: Vec<String>,
}

#[derive(Args)]
struct BusArgs{
    /// Busfile to read the CB/UMIs from
    #[clap()]
    busfile: String,
    
    /// Max #entries to consider in the bus file
    #[clap(short = 'n', long= "nmax")] 
    nmax: usize
}


#[derive(Args)]
struct BusArgs2{
    /// Busfolder to read the CB/UMIs from
    /// contains the busfile and ec.matrix, transcripts
    #[clap()]
    busfolder: String,

    #[clap()]
    t2gfile: String,
    
    /// Max #entries to consider in the bus file
    #[clap(short = 'n', long= "nmax")] 
    nmax: usize
}

#[derive(Args)]
struct TSOArgs{
    /// List of fastq files
    #[clap()]
    fastq_list: Vec<String>,
}

#[derive(Args)]
struct BusMergeArgs{
    #[clap(long= "i1")] 
    inbus1: String,
    #[clap(long= "i2")] 
    inbus2: String,

    #[clap(long= "o1")] 
    outbus1: String,
    #[clap(long= "o2")] 
    outbus2: String,  
}

fn main() {
    // myfastq::run();
    // sqlite::run();
    // hset::run();
    
    let cli = Cli::parse();

    match cli.command{
        MyCommand::cb(args) => {
            println!("Doing CB only");
            cb_errors::run(&args.fastq_list, args.whitelist, cli.output, args.topn)           
        },
        MyCommand::cb_umi_exact(args) => {
            println!("Doing CB_UMI sketch");
            println!("WARNING: MEMORY INTENSIVE!!");
            cb_umi_errors::run(&args.fastq_list, args.whitelist, cli.output, args.topn)
        },
        MyCommand::cb_umi_sketch(args) => {
            println!("Doing CB_UMI sketch");
            cb_umi_sketch::run_top_n(&args.fastq_list, args.whitelist, cli.output, args.topn)
            
        },
        MyCommand::cb_umi_cell(args) => {
            println!("Doing CB_UMI via single cells");
            cb_umi_per_cell::run(&args.busfile, &cli.output, args.nmax, false)  
        }
        MyCommand::cb_umi_cell_aggr(args) => {
            println!("Doing AGGREGATE CB_UMI via single cells");
            cb_umi_per_cell::run(&args.busfile, &cli.output, args.nmax, true)
        }
        MyCommand::busmerge(args) => {
            println!("Doing bus merging");
            busmerger::merge_busfiles_on_overlap(args.inbus1, args.inbus2, args.outbus1, args.outbus2)      
        }
        MyCommand::tso_error(args) => {
            println!("Doing TSO error");
            tso_error::run(&args.fastq_list, cli.output)      
        }              
        MyCommand::cb_umi_cell_gene(args) => {
            println!("Doing CUG error");
            cb_umi_per_cell_gene::run(args.busfolder, &cli.output, args.nmax, false, args.t2gfile)      
        }        
    };


    // println!("Whitelist {:?}",cli.whitelist);
    // println!("Output {:?}",cli.output);
    // println!("Top N {:?}",cli.topn);
    // println!("FASTQ {:?}",cli.fastq_list);

    // cb_umi_sketch::run_topN(&cli.fastq_list, cli.whitelist, cli.output, topn);
    // cb_umi_errors::run(&cli.fastq_list, cli.whitelist, "/tmp/full_out.csv".to_string(), topn);



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
}



