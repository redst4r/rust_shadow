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
mod phred_counter;
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
    cb_umi_cell_gene(BusArgs2),
    busmerge(BusMergeArgs),
    tso_error(TSOArgs),
    phred(PhredArgs)
}

#[derive(Args)]
struct PhredArgs{
    /// List of fastq files
    #[clap()]
    fastq_list: Vec<String>,
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
    #[clap(long= "busfile")] 
    busfile: String,
    
    /// Max #entries to consider in the bus file
    #[clap(short = 'n', long= "nmax")] 
    nmax: usize,

    /// whether to aggregate the entries on a cell level (one row per cell), or ouput one row per UMI
    #[clap(long= "aggr")] 
    aggregate: bool,
}


#[derive(Args)]
struct BusArgs2{
    /// Busfolder to read the CB/UMIs from
    /// contains the busfile and ec.matrix, transcripts
    #[clap(long= "busfolder")] 
    busfolder: String,

    /// transcript to gene filename
    #[clap(long= "t2g")] 
    t2gfile: String,
    
    /// whether to aggregate the entries on a cell level (one row per cell), or ouput one row per UMI
    #[clap(long= "aggr")] 
    aggregate: bool,

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
            cb_umi_per_cell::run(&args.busfile, &cli.output, args.nmax, args.aggregate)  
        }
        MyCommand::busmerge(args) => {
            println!("Doing bus merging");
            busmerger::merge_busfiles_on_overlap(args.inbus1, args.inbus2, args.outbus1, args.outbus2)      
        }
        MyCommand::tso_error(args) => {
            println!("Doing TSO error");
            tso_error::run(&args.fastq_list, cli.output)     
        } 
        MyCommand::phred(args) => {
            println!("Doing CUG error");
            phred_counter::run(&args.fastq_list, cli.output)      
        }                   
        MyCommand::cb_umi_cell_gene(args) => {
            println!("Doing CUG error");
            cb_umi_per_cell_gene::run(args.busfolder, &cli.output, args.nmax, args.aggregate, args.t2gfile)      
        }        
    };
}



