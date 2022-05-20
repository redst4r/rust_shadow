// use std::io::{self, Read};
use std::io::BufReader;
use std::io::BufRead;
// use flate2::read::GzDecoder;
use rust_htslib::bgzf;
use std::collections::{HashMap, HashSet};
// use rusqlite::{Connection, Result};
// use rusqlite::NO_PARAMS;
use counter::Counter;
use bktree::{BkTree, levenshtein_distance};

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
    // TODO: size check
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

pub fn count_cb_umi(fname: String) -> Counter<(String, String), i32> {
    // coutns the CB/UMI pairs in the fastq

    // reading the fastq.gz
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines()
        .enumerate().filter(|x| x.0 % 4 == 1)
        .map(|x| x.1)
        ;//.take(10_000_000);


    // parsing the lines, counting
    // let mut countermap: HashMap<(String, String), i32> = HashMap::new();
    let mut countermap: Counter<(String, String), i32> = Counter::new();

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

pub fn top_n(counter: &Counter<(String, String), i32>, n: i32) -> Vec<(String, String)>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements
    let mut bk: BkTree<String> = BkTree::new(levenshtein_distance);

    let mut c = 0;
    for ((cb, umi), _freq) in counter.most_common(){

        let cb_umi = format!("{cb}_{umi}");
        if bk.find(cb_umi, 2).len() == 0{
            let cb_umi = format!("{cb}_{umi}");
            bk.insert(cb_umi);
            c += 1;
        }
        if c >= n{
            break
        }
    }

    // for now, all the frequent CB/UMIs are the elements of the BKTree
    // put them into a vector
    let mut top2: Vec<(String,String)> = Vec::new();

    for element in bk.into_iter(){
        let mut cb_umi = element.split('_');
        let cb = cb_umi.next().unwrap();
        let umi = cb_umi.next().unwrap();
        top2.push((String::from(cb), String::from(umi)));
    }
    top2
}


pub fn find_shadows(cb_umi: (String, String), filtered_map: &Counter<(String, String), i32>) -> HashMap<usize, i32>{
    // for a given "true" CB/UMI find the number of shadowed reads (i.e. reads with a single subsitution)
    // and count their number (position specific)
    // we get a dcitionary with position -> #shadow reads
    //
    let (cb, umi) = cb_umi;

    let mut all_muts: Vec<(String, String, usize)>= Vec::new();  // cb, umi, pos
    for pos in 0..umi.len() {
        for m in get_1bp_mutations(&umi, pos){
            let c2 = cb.clone();
            all_muts.push((c2, m, pos));
        }
    }
    // println!("{:?}", all_muts);

    let mut n_shadows_per_pos: HashMap<usize, i32> = HashMap::new();
    for (cb, umi, pos) in all_muts{
        if let Some(f) = filtered_map.get(&(cb, umi)){
            // for that shadow cb/umi, how often did we acutally see it
            n_shadows_per_pos.insert(pos, *f);
        }
        else{
            n_shadows_per_pos.insert(pos, 0);
        }
    }
    n_shadows_per_pos
}


fn get_1bp_mutations(seq: &String, pos: usize) -> Vec<String>{
    // return the three 1BP mutations of a sequence at the given position"
    let mut muts = Vec::with_capacity(3);

    for base in ['A', 'C', 'G', 'T']{
        let prefix = &seq[0..pos];
        let postfix = &seq[(pos+1)..];
        let m = format!("{prefix}{base}{postfix}");
        if m != *seq{
            muts.push(m);
        }
    }
    return muts
}