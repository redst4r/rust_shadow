use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bgzf;
use std::collections::{HashSet, HashMap};
use counter::Counter;
use crate::cb_umi_errors::{parse_r1, get_1bp_mutations, parse_whitelist_gz};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use std::fs::File;
use streaming_algorithms::{CountMinSketch, Top};


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

pub fn run_top(fastq_list: Vec<String>){

    // determinig the topN barcodes in the fastqs

    let file_iterators = fastq_list.into_iter()
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

    let TOPN = 10000; 
    let tolerance = 1e-8;
    let probability = 0.00001;


    let mut ccc:Top<String, u32> = Top::new(TOPN, probability, tolerance, {});
    let mut countermap: Counter<String, u32> = Counter::new();

    for (i, l) in my_iter.enumerate(){
        if let Ok(line) = l{
            if let Some((cb, _umi)) = parse_r1(line){
                let cb_umi = format!("{cb}_{_umi}");
                let cb_umi2 = format!("{cb}_{_umi}");

                // let cb2 = cb.clone();
                let counter = countermap.entry(cb_umi).or_insert(0);
                *counter += 1;

                ccc.push(cb_umi2, &1);
            }
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }        
    }

    let approx_top: HashSet<String> = ccc.iter()
        .map(|(seq, _freq)| (*seq).clone())
        .collect();

    // compare true to estimated topN
    let true_top: HashSet<String> = countermap.most_common().iter().take(TOPN).map(|(seq, _freq)| (*seq).clone()).collect();

    set_comparison(&true_top, &approx_top);

    // let aaa: Vec<(String, u32)> = countermap.most_common().into_iter().take(TOPN).collect();
    // println!("True top N {:?}", aaa);

    // println!("True top N {:?}", true_top);
    // println!("appr top N {:?}", approx_top);

    //write to file
    // first merge the cbs
    let mut true_top_dict: HashMap<String, u32> = HashMap::new();
    
    for (k,v) in ccc.iter(){
        true_top_dict.insert(k.clone() ,*v);
    }


    let mut true_freqs: Vec<u32> = Vec::new();
    let mut approx_freqs: Vec<u32> = Vec::new();

    for item in true_top.union(&approx_top).into_iter(){
        let f1 = countermap.get(item).cloned().unwrap_or(0);
        let f2 = true_top_dict.get(item).cloned().unwrap_or(0);
        true_freqs.push(f1);
        approx_freqs.push(f2)
    }

    let df_true = Series::new("true_frequency", true_freqs);
    let df_approx = Series::new("approx_frequency", approx_freqs);
    
    let mut df_final = DataFrame::new(vec![df_true, df_approx]).unwrap();
    println!("{:?}", df_final);

    // write to CSV
    let mut output_file: File = File::create("/tmp/topN.csv".to_string()).unwrap();
    CsvWriter::new(&mut output_file)
        .has_header(true)
        .finish(&mut df_final)
        .unwrap();  


}


struct GreaterThan1Bloom {
    minsketch : CountMinSketch<String, u32>,  // for storing items seen once
    greater_than_1_items : HashSet<String>
}
impl GreaterThan1Bloom {
    fn add_item(&mut self, item:&String){

        let current_freq = self.minsketch.get(item);
        // println!("Adding item {item}, current freq = {current_freq}");

        if  current_freq == 0{
            self.minsketch.push(item, &1);
        }
        else{
            self.greater_than_1_items.insert(item.clone());
        }

        // println!("Current #>1 {}", self.greater_than_1_items.len())
    }
}

pub fn run_gt1(fastq_list: Vec<String>){

    // get all CB_UMI that occur more than once

    let file_iterators = fastq_list.into_iter()
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

    let tolerance = 1e-8;
    let probability = 0.00001;


    let mut ccc:CountMinSketch<String, u32> = CountMinSketch::new(probability, tolerance, {});
    let mut v: HashSet<String> = HashSet::new();

    let mut GE = GreaterThan1Bloom {
        minsketch : ccc,
        greater_than_1_items : v
    };

    let mut countermap: Counter<String, u32> = Counter::new();

    for (i, l) in my_iter.enumerate(){
        if let Ok(line) = l{
            if let Some((cb, _umi)) = parse_r1(line){

                let cb_umi = format!("{cb}_{_umi}");
                let cb_umi2 = format!("{cb}_{_umi}");

                let counter = countermap.entry(cb_umi).or_insert(0);
                *counter += 1;

                GE.add_item(&cb_umi2);
            }
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }        
    }

    // true number of >1 elemets
    println!("Memory savings len(counter) {}", countermap.len());


    let f: Vec<String> = countermap.into_iter()
        .filter(|(k,v)| v>&1).map(|(k, v)| k)
        .collect();
    let fsize = f.len();

    println!("True #>1 {}  Estimated #>1 {}",fsize, GE.greater_than_1_items.len());


    let x: HashSet<String> = f.into_iter().collect();
    let y: HashSet<String> = GE.greater_than_1_items.into_iter().collect();
    set_comparison(&x, &y)

    // println!("True #>1 {:?} ", f);
    // println!("Est #>1 {:?} ", GE.greater_than_1_items);

}

use core::hash::Hash;
fn set_comparison<T>(set_a: &HashSet<T>,set_b: &HashSet<T>)
    where T: Eq + Hash
{
    // warning: this instantiates all the sets
    // instead, just iterate over them
    let intersect = set_a.intersection(set_b).collect::<Vec<&T>>();
    let AminusB = set_a.difference(set_b).collect::<Vec<&T>>();
    let BminusA = set_b.difference(set_a).collect::<Vec<&T>>();
    println!("|A&B| {} |A-B| {}  |B-A|{}", intersect.len(), AminusB.len(), BminusA.len());

    // let n_intersect = set_a.intersection(set_b).fold(0, |accum, _item| accum+1);
    // let n_AminusB = set_a.difference(set_b).fold(0, |accum, _item| accum+1);
    // let n_BminusA = set_b.difference(set_a).fold(0, |accum, _item| accum+1);

    let n_intersect = set_a.intersection(set_b).count();
    let n_AminusB = set_a.difference(set_b).count();
    let n_BminusA = set_b.difference(set_a).count();

    println!("|A&B| {} |A-B| {}  |B-A|{}", n_intersect, n_AminusB, n_BminusA);

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