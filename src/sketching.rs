use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bgzf;
use std::collections::{HashMap};
use counter::Counter;
use bktree::{BkTree, levenshtein_distance};
use crate::cb_umi_errors::{parse_r1, get_1bp_mutations, parse_whitelist_gz};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use std::fs::File;
use streaming_algorithms::{CountMinSketch};


//
//  n = 400M
//  c - c_true ~ 1   +> epsilon = 1/400M ~ 1e-8
//  delta = 1%
//
// It follows
//  m = 2.71/eps = 2.71/1e-8  = 10**9
//  k = ln(1/delta) = 100

pub fn count_cb_filelist(fname_list: Vec<String>) -> Counter<String, i32> {
    // coutns the CB/UMI pairs in the fastqs


    // reading the fastq.gz
    // we chain all those files together into a single iterator
    let file_iterators = fname_list.into_iter()
        .map(|fname|{
            let decoder = bgzf::Reader::from_path(fname).unwrap();
            let reader = BufReader::new(decoder);
            let my_iter = reader.lines()
                .enumerate().filter(|x| x.0 % 4 == 1)
                .map(|x| x.1);
            my_iter
        }
        );

    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);


    // use CountMinSketch for approximate freqs
    let mut ccc:CountMinSketch<String, u32> = CountMinSketch::new(0.01, 1e-4, {});

    // parsing the lines, counting
    let mut countermap: Counter<String, i32> = Counter::new();

    for (i, l) in my_iter.enumerate(){
        if let Ok(line) = l{
            if let Some((cb, _umi)) = parse_r1(line){
                let counter = countermap.entry(cb).or_insert(0);
                *counter += 1
            }
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}


pub fn run(fastq_list: Vec<String>){


    // let a = vec!["MS", "MS", "MS", "AK", "AK"];

    // let fastq_file1: String = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Fresh1/Fresh1_CKDL210025651-1a-SI_TT_C2_HVWMHDSX2_S2_L001_R1_001.fastq.gz".into();
    // let fastq_list = vec![fastq_file1] ;

    let file_iterators = fastq_list.into_iter()
        .map(|fname|{
            let decoder = bgzf::Reader::from_path(fname).unwrap();
            let reader = BufReader::new(decoder);
            let my_iter = reader.lines()
                .enumerate().filter(|x| x.0 % 4 == 1)
                .map(|x| x.1).take(10_000_000);
            my_iter
        }
        );

    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);



    let mut ccc:CountMinSketch<String, u32> = CountMinSketch::new(0.0001, 1e-9, {});
    let mut countermap: Counter<String, u32> = Counter::new();

    for (i, l) in my_iter.enumerate(){
        if let Ok(line) = l{
            if let Some((cb, _umi)) = parse_r1(line){
                let cb2 = cb.clone();
                let counter = countermap.entry(cb).or_insert(0);
                *counter += 1;

                ccc.push(&cb2, &1);
            }
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }        

    }

    let mut n_correct = 0;
    let mut n_incorrect = 0;

    let mut true_freqs: Vec<u32> = Vec::new();
    let mut approx_freqs: Vec<u32> = Vec::new();


    for (seq, freq) in countermap.iter(){
        let f_approx = ccc.get(seq);

        if f_approx != *freq{
            println!("Mismatch {seq} {f_approx}-{freq}");
            n_incorrect +=1;
        }
        else{
            n_correct+=1;
        }
        true_freqs.push(*freq);
        approx_freqs.push(f_approx)
    }


    println!("Correct {n_correct} Incorrect {n_incorrect}");

    
    let df_true = Series::new("true_frequency", true_freqs);
    let df_approx = Series::new("approx_frequency", approx_freqs);
    
    let mut df_final = DataFrame::new(vec![df_true, df_approx]).unwrap();
    println!("{:?}", df_final);

    // write to CSV
    let mut output_file: File = File::create("/tmp/CMS.csv".to_string()).unwrap();
    CsvWriter::new(&mut output_file)
        .has_header(true)
        .finish(&mut df_final)
        .unwrap();        
    
}