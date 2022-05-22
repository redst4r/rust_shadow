// use std::io::{self, Read};
use std::io::BufReader;
use std::io::BufRead;
// use flate2::read::GzDecoder;
use rust_htslib::bgzf;
use std::collections::{HashMap};
use counter::Counter;
use bktree::{BkTree, levenshtein_distance};
use polars::prelude::{DataFrame, NamedFrom, Series};
use crate::utils::{parse_whitelist_gz,parse_r1_struct, write_to_csv, all_mutations_for_cbumi, CbUmi};


pub fn count_cb_filelist(fname_list: &Vec<String>) -> Counter<CbUmi, u32> {
    // coutns the CB/UMI pairs in the fastq

    // reading the fastq.gz
    // we chain all those files together into a single iterator
    let file_iterators = fname_list.into_iter()
        .map(|fname|{
            let decoder = bgzf::Reader::from_path(fname).unwrap();
            let reader = BufReader::new(decoder);
            let my_iter = reader.lines()
                .enumerate().filter(|x| x.0 % 4 == 1)
                .map(|x| x.1)
                .filter_map(|line| line.ok()); //takes care of errors in file reading
            my_iter
        }
        );

    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);


    // parsing the lines, counting
    // let mut countermap: HashMap<CbUmi, i32> = HashMap::new();
    let mut countermap: Counter<CbUmi, u32> = Counter::new();

    for (i, line) in my_iter.enumerate(){
        if let Some(cbumi) = parse_r1_struct(line){
            let counter = countermap.entry(cbumi).or_insert(0);
            *counter += 1
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}

pub fn top_n(counter: &Counter<CbUmi, u32>, n: usize) -> Vec<CbUmi>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements
    let mut bk: BkTree<String> = BkTree::new(levenshtein_distance);

    let mut c = 0;
    for (cbumi, _freq) in counter.most_common(){

        let cb_umi_str = cbumi.to_string();
        let matches = bk.find(cb_umi_str, 1);
        if matches.len() == 0{
            let cb_umi_str =  cbumi.to_string();
            bk.insert(cb_umi_str);
            c += 1;
        }

        if c >= n{
            break
        }
        if c % 1_000 == 0{
            println!("Iteration {c} of {n}");
        }          
    }

    // for now, all the frequent CB/UMIs are the elements of the BKTree
    // put them into a vector
    let mut top2: Vec<CbUmi> = Vec::new();

    for element in bk.into_iter(){
        let cb_umi = CbUmi::from_string(&element);
        top2.push(cb_umi);
    }
    top2
}

pub fn find_shadows(cb_umi: CbUmi, filtered_map: &Counter<CbUmi, u32>) -> HashMap<usize, u32>{
    // for a given "true" CB/UMI find the number of shadowed reads (i.e. reads with a single subsitution)
    // and count their number (position specific)
    // we get a dcitionary with position -> #shadow reads
    //
    
    // let (cb_orig, umi_orig) = (cb_umi.cb, cb_umi.umi);//todo why is this not linted as unused

    let all_muts = all_mutations_for_cbumi(cb_umi);

    // let mut all_muts: Vec<(String, String, usize)>= Vec::new();  // cb, umi, pos
    // for pos in 0..umi_orig.len() {
    //     for m in get_1bp_mutations(&umi_orig, pos){
    //         let c2 = cb_orig.clone();
    //         all_muts.push((c2, m, pos));
    //     }
    // }


    let mut n_shadows_per_pos: HashMap<usize, u32> = HashMap::new();
    for (cbumi, mutated_pos) in all_muts{

        if let Some(f) = filtered_map.get(&cbumi){
            // for that shadow cb/umi, how often did we acutally see it
            // remember, theres 3 mutation at each position so we have to sum up
            let current_freq = n_shadows_per_pos.entry(mutated_pos).or_insert(0);
            *current_freq += f;
        }
        else{
            n_shadows_per_pos.entry(mutated_pos).or_insert(0);
        }
    }
    n_shadows_per_pos
}


pub fn run(fastq_list: &Vec<String>, whitelist_file: String, output_csv_file: String, topn:usize){


    // parse whitelist
    let whitelist = parse_whitelist_gz(whitelist_file);
    println!("Whitelist len {}", whitelist.len());
    
    let mut countmap = count_cb_filelist(fastq_list);
    // transform into shadow counter
    println!("len of counter {}", countmap.len());

    // filter for whitelist only entries
    println!("Filtering for whilelist");

    // let mut keys_to_remove: Vec<CbUmi> = Vec::new();
    // for (cbumi, _count) in countmap.iter(){
    //     if !whitelist.contains(&cbumi.cb){
    //         keys_to_remove.push(cbumi.clone());
    //     }
    // }

    // for k in keys_to_remove{
    //     countmap.remove(&k);
    // }

    println!("len of filtered counter {}", countmap.len());
    // unstable:
    // countmap.drain_filter(|k, v| {
    //     k.0 == k.1
    // });

    // now the hard part: group by CB, look at all UMIs therein
    println!("calculating most common");
    let most_common: Vec<CbUmi> = top_n(&countmap, topn);
    println!("most common {:?}", most_common.len());

    // for each entry in most common, find potential shadows
    // add them to polars
    let mut polars_data: HashMap<String, Vec<u32>> = HashMap::new();
    let mut cellnames: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();

    for mc in most_common{

        if !whitelist.contains(&mc.cb){
            continue
        }

        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &countmap);
        for (position, n_shadows) in nshadows_per_position.iter(){
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = countmap.get(&mc2).unwrap();
        polars_data.entry("total".into()).or_insert(vec![]).push(*total_counts);

        cellnames.push(mc2.cb);
        umis.push(mc2.umi);
    }
    
    // to polars dataframe
    let df = DataFrame::new(
        polars_data.into_iter()
            .map(|(name, values)| Series::new(&format!("{name}"), values))
            .collect::<Vec<_>>()).unwrap();

    let df_cb = Series::new("CB", cellnames);
    let df_umi = Series::new("UMI", umis);

    let mut df_final = df.hstack(&[df_cb, df_umi]).unwrap();    
    println!("{:?}", df_final);

    // write to CSV
    write_to_csv(&mut df_final, output_csv_file);    
}
