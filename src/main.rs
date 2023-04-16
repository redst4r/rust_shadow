
use std::io::Write;
use std::fs::File;
use std::time::Instant;
use rustbustools::busmerger;
use clap::{self, Parser, Subcommand, Args};
use rustfastq::{cb_errors, cb_umi_errors, cb_umi_sketch, cb_umi_per_cell, phred_counter, tso_error, cb_umi_per_cell_gene, io::test_filter};

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
    phred(PhredArgs),
    count(CountArgs),
    filter,
}

#[derive(Args)]
struct CountArgs{
    /// List of fastq files
    #[clap()]
    fastq_list: Vec<String>,
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
            busmerger::merge_busfiles_on_overlap(&args.inbus1, &args.inbus2, &args.outbus1, &args.outbus2)      
        }
        MyCommand::tso_error(args) => {
            println!("Doing TSO error");
            tso_error::run(&args.fastq_list, cli.output)     
        } 
        MyCommand::phred(args) => {
            println!("Doing Phred Counter");
            phred_counter::run(&args.fastq_list, cli.output)      
        }                   
        MyCommand::cb_umi_cell_gene(args) => {
            println!("Doing CUG error");
            cb_umi_per_cell_gene::run(args.busfolder, &cli.output, args.nmax, args.aggregate, args.t2gfile)      
        }    
        MyCommand::count(args) => {
            println!("Doing counting");
            let mut file_handle = File::create(cli.output).unwrap();

            for filename in args.fastq_list{
                println!("Counting {}", filename.clone());

                let now = Instant::now();
                let c = count_fastq_reads(filename.clone());
                let elapsed_time = now.elapsed();
                println!("Counted {}, took {} minutes.", filename.clone(), elapsed_time.as_secs()/60);                

                // write result to filen
                file_handle.write_all(format!("{}\t{}\n", filename, c).as_bytes()).unwrap();
            }
        }
        MyCommand::filter => {
            test_filter();
        },
    };
}

pub fn count_fastq_reads(filename: String) -> usize{
    // count the nubmer of entries (not lines!) in the fastq
    let count = rustfastq::io::fastq_list_iter(&[filename], 1).count();
    count
}