// use std::io::{self, Read};
use std::io::BufReader;
use std::io::BufRead;
// use flate2::read::GzDecoder;
use rust_htslib::bgzf;
use std::collections::{HashMap, HashSet};
// use rusqlite::{Connection, Result};
// use rusqlite::NO_PARAMS;


pub fn parse_whitelist_gz(fname: String) -> HashSet<String>{
    // loading 10x CB whilelist from file
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines();//.take(10_000_000);
    let mut hset: HashSet<String> = HashSet::new();
    for l in my_iter{
        if let Ok(line) = l{
            hset.insert(line);
        }
    }
    hset
}

fn parse_r1(seq: String) -> Option<(String, String)>{
    // if seq.len() == 16 + 12{
        let cb = (&seq[0..16]).into();
        let umi = (&seq[16..28]).into();
        Some((cb, umi))
    // }
    // else
    // {
        // None
    // }
}

pub fn count_cb_umi(fname: String) -> HashMap<(String, String), i32> {
    // coutns the CB/UMI pairs in the fastq

    // reading the fastq.gz
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines()
        .enumerate().filter(|x| x.0 % 4 == 1)
        .map(|x| x.1)
        ;//.take(10_000_000);


    // parsing the lines, counting
    let mut countermap: HashMap<(String, String), i32> = HashMap::new();
    for (i, l) in my_iter.enumerate(){
        if let Ok(line) = l{
            if let Some((cb, umi)) = parse_r1(line){
                let counter = countermap.entry((cb, umi)).or_insert(0);
                *counter += 1
            }
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}

