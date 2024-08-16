
use std::io::BufWriter;
use std::io::Write;
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use clap::{self, Parser, Subcommand, Args};
use rustfastq::demultiplex;
use rustfastq::demultiplex::demux_dual_index;
use rustfastq::demultiplex::samplesheet_to_hashmap;
use rustfastq::demultiplex::samplesheet_to_hashmap_2;
use rustfastq::demultiplex::DualIndex;
use rustfastq::demultiplex::Samplename;
use rustfastq::demultiplex::Samplesheet;
use rustfastq::utils::get_spinner;
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
    count_sampleix(SampleIxArgs),
    demux_dual(DemuxDualArgs),
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
struct SampleIxArgs{
    /// List of fastq files
    #[clap(long= "i1")]
    i1_list: Vec<String>,
    #[clap(long= "i2")]
    i2_list: Vec<String>,
}

#[derive(Args)]
struct DemuxDualArgs{
    /// List of fastq files
    #[clap(long= "i1")]
    i1_list: Vec<String>,
    #[clap(long= "i2")]
    i2_list: Vec<String>,
    #[clap(long= "r1")]
    r1_list: Vec<String>,
    #[clap(long= "r2")]
    r2_list: Vec<String>,    

    #[clap(long= "samplesheet")]
    samplesheet: PathBuf,        
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

        /*
        cargo run --release -- -o /tmp/six.txt  count-sampleix  --i1 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_I1_001.fastq.gz  --i2 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_I2_001.fastq.gz 
         */
        MyCommand::count_sampleix(args) => {
            let count_map = paired_index_counter(args.i1_list, args.i2_list);
            let mut count_vec: Vec<((String, String), usize)> = count_map.into_iter().collect();
            count_vec.sort_by(|a,b| b.1.cmp(&a.1));

            // write it to the file
            let mut fh = BufWriter::new(File::create(cli.output).unwrap());

            for ((s1, s2), c) in  count_vec.iter() {
                    writeln!(fh, "{},{},{}", s1,s2,c).unwrap();
            }
        },
        /*
        cargo run --release -- -o /tmp/six.txt  demux-dual  
        --i1 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_I1_001.fastq.gz  
        --i2 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_I2_001.fastq.gz  
        --r1 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_R1_001.fastq.gz 
        --r2 /home/michi/tuba_mount/sequencing_data/IR-BL-001/fastq/01.RawData/Undetermined/Undetermined_Undetermined_22C7YHLT4_S0_L004_R2_001.fastq.gz 
        --samplesheet ttt
         */
        MyCommand::demux_dual(args) => {

            let outdir = Path::new(&cli.output);
            assert!(outdir.exists() && outdir.is_dir());

            // let sheet: HashMap<_,_> = vec![
            //     (DualIndex("GGGGGGGGGG".to_string(), "AGATCTCGGT".to_string()), Samplename("AGATCTCGGT".to_string())),
            //     (DualIndex("GGGGGGGGGG".to_string(), "ATAGATGCTC".to_string()), Samplename("ATAGATGCTC".to_string())),
            //     (DualIndex("GGGGGGGGGG".to_string(), "GTAACAGGAA".to_string()), Samplename("GTAACAGGAA".to_string())),
            // ].into_iter().collect();
            // let samplesheet = crate::demultiplex::Samplesheet::new(sheet);

            // let samplesheet = samplesheet_to_hashmap_2();
            let samplesheet = Samplesheet::from_csv(&args.samplesheet);

            let unassigned = "UA".to_string();
            demultiplex::demux_dual_index_2(samplesheet,unassigned, 
                args.i1_list, args.i2_list, args.r1_list, args.r2_list, outdir);


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

/// iterate through the index1/index2 reads and count the frequency of sample-barcode-pairs
pub fn paired_index_counter(i1_list: Vec<String>, i2_list: Vec<String>) -> HashMap<(String, String), usize> {
   
    let i1 = rustfastq::io::fastq_list_iter(&i1_list);
    let i2 = rustfastq::io::fastq_list_iter(&i2_list);

    let mut counter: HashMap<(String, String), usize> = HashMap::new();

    let bar = get_spinner();

    for (i,(f1, f2)) in izip!(i1, i2).enumerate() {
        let c = counter.entry((f1.seq, f2.seq)).or_insert(0);
        *c += 1;

        if i % 1_000_000 ==0 {
            bar.inc(1_000_000)
        }
    }
    counter
}

#[test]
fn test_paired() {
    let count_map = paired_index_counter(
        vec!["/home/michi/mounts/myDrive/230601_VH00715_118_AACVG5JM5_fastq/Undetermined_S0_L001_I1_001.fastq.gz".to_string()], 
        vec!["/home/michi/mounts/myDrive/230601_VH00715_118_AACVG5JM5_fastq/Undetermined_S0_L001_I2_001.fastq.gz".to_string()], 
    );

        // Get a sorted (by field 0 ("count") in reversed order) list of the
    // most frequently used indices:
    let mut count_vec: Vec<((String, String), usize)> = count_map.into_iter().collect();
    count_vec.sort_by(|a,b| b.1.cmp(&a.1));

    for ((s1, s2), c) in  count_vec.iter().take(20) {
        // if c > 10000 {
            println!("{}_{}:{}", s1,s2,c)
        // }
    }
}


/*
cargo run --release -- 

*/

#[test]
fn test_bamfile_RG() {
    use std::io::{BufRead, BufReader};

    // let my_headerline = "@RG     ID:E14C_2-747406_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0:0:1::        SM:E14C_2-747406_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0      LB:0.1  PU:E14C_2-747406_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0:0:1::        PL:ILLUMINA";
    let fh = BufReader::new(
        File::open("/tmp/header.sam").unwrap()
    );
    for line in fh.lines()
        .map(|r| r.unwrap())
        .filter(|l| l.starts_with("@RG")) {
        
        println!("Line: {}", line);

        let y = parse_rg_line(&line);
        println!("Parse: {:?}", y);

    }
    use regex::Regex;
    use std::str::FromStr;
    // from 10x
    fn parse_rg_line(line: &str) -> Option<(String, (String, u32))> {
        let mut entries = line.split('\t');
        entries.next()?; // consume @RG entry

        let mut tags = entries
            .map(|entry| entry.split_once(':').unwrap())
            .collect::<HashMap<_, _>>();

        println!("tags: {:?}", tags);


        let v = tags.remove("ID")?;
        let (rg, lane) = v.rsplit_once(':')?;
        println!("rg:{:?} lane:{:?}", rg, lane);
        // v.rsplitn(2, ':');
        println!("{:?}", u32::from_str(lane));

        let result = match u32::from_str(lane) {
            Ok(n) => Some((v.to_string(), (rg.to_string(), n))),
            Err(_) => {
                // Handle case in ALIGNER pipeline prior to 2.1.3 -- samtools merge would append a unique identifier to each RG ID tags
                // Detect this condition and remove from lane
                let re = Regex::new(r"^([0-9]+)-[0-9A-F]+$").unwrap();
                let cap = re.captures(lane)?;
                let lane_u32 = u32::from_str(cap.get(1).unwrap().as_str()).unwrap();
                Some((v.to_string(), (rg.to_string(), lane_u32)))
            }
        };
        println!("Final return: {:?}", result);

        result

    }

}
