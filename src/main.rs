
use std::io::Write;
use std::fs::File;
use std::time::Instant;
use clap::{self, Parser, Subcommand, Args};
use rustfastq::{phred_counter, io::quality_filter};

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
    phred(PhredArgs),
    count(CountArgs),
    qcfilter(QCFilterArgs),
}

#[derive(Args)]
struct QCFilterArgs{
    /// List of fastq files
    #[clap()]
    fastq_file: String,
    /// QC score threhold: any read with less will be dropped
    #[clap(short = 'q', long= "qcscore")] 
    qcscore: f32, 
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

fn main() {   
    let cli = Cli::parse();

    match cli.command{
        MyCommand::phred(args) => {
            println!("Doing Phred Counter");
            phred_counter::run(&args.fastq_list, cli.output)      
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
        MyCommand::qcfilter(args) => {
            quality_filter(&args.fastq_file, &cli.output, args.qcscore)
        },
    };
}

pub fn count_fastq_reads(filename: String) -> usize{
    // count the nubmer of entries (not lines!) in the fastq
    let count = rustfastq::io::fastq_list_iter(&[filename]).count();
    count
}

use std::collections::HashMap;
use itertools::izip;

pub fn paired_index_counter(i1file: String, i2file: String) {

    let i1_list = [i1file];
    let i2_list = [i2file];
    
    let i1 = rustfastq::io::fastq_list_iter(&i1_list);
    let i2 = rustfastq::io::fastq_list_iter(&i2_list);

    let mut counter: HashMap<(String, String), usize> = HashMap::new();

    for (f1, f2) in izip!(i1, i2) { //.take(1_000_000)
        let c = counter.entry((f1.seq, f2.seq)).or_insert(0);
        *c += 1;
    }

    // Get a sorted (by field 0 ("count") in reversed order) list of the
    // most frequently used indices:
    let mut count_vec: Vec<((String, String), usize)> = counter.into_iter().collect();

    count_vec.sort_by(|a,b| b.1.cmp(&a.1));

    for ((s1, s2), c) in  count_vec.iter().take(20) {
        // if c > 10000 {
            println!("{}_{}:{}", s1,s2,c)
        // }
    }
}

#[test]
fn test_paired() {
    paired_index_counter(
        "/home/michi/mounts/myDrive/230601_VH00715_118_AACVG5JM5_fastq/Undetermined_S0_L001_I1_001.fastq.gz".to_string(), 
        "/home/michi/mounts/myDrive/230601_VH00715_118_AACVG5JM5_fastq/Undetermined_S0_L001_I2_001.fastq.gz".to_string(), 
    );
}