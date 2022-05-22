use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bgzf;
use std::collections::{HashMap};
use counter::Counter;
use bktree::{BkTree, levenshtein_distance};
use crate::utils::{parse_r1, get_1bp_mutations, parse_whitelist_gz, write_to_csv};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use std::fs::File;


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
                .map(|x| x.1)
                .filter_map(|line| line.ok()) //takes care of errors in file reading
                ;
            my_iter
        }
        );

    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);

    // parsing the lines, counting
    let mut countermap: Counter<String, i32> = Counter::new();

    for (i, line) in my_iter.enumerate(){
        if let Some((cb, _umi)) = parse_r1(line){
            let counter = countermap.entry(cb).or_insert(0);
            *counter += 1
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}


pub fn top_n(counter: &Counter<String, i32>, n: usize) -> Vec<String>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements
    let mut bk: BkTree<String> = BkTree::new(levenshtein_distance);

    let mut c = 0;

    for (cb, _freq) in counter.most_common(){
        let cb2 = cb.clone();  // TODO stupid cloing to be able to insert
        if bk.find(cb, 2).len() == 0{
            // let cb_umi = format!("{cb}_{umi}");
            bk.insert(cb2);   // TODO stupid cloing to be able to insert
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
    let mut top2: Vec<String> = Vec::new();

    for cb in bk.into_iter(){
        top2.push(cb);
    }
    top2
}

pub fn find_shadows(cb: String, filtered_map: &Counter<String, i32>) -> HashMap<usize, i32>{
    // for a given "true" CB find the number of shadowed reads (i.e. reads with a single subsitution)
    // and count their number (position specific)
    // we get a dcitionary with position -> #shadow reads
    //

    let mut all_muts: Vec<(String, usize)>= Vec::new();  // cb, umi, pos
    for pos in 0..cb.len() {
        for m in get_1bp_mutations(&cb, pos){
            all_muts.push((m, pos));
        }
    }

    let mut n_shadows_per_pos: HashMap<usize, i32> = HashMap::new();
    for (cb, pos) in all_muts{
        if let Some(f) = filtered_map.get(&cb){
            // update our count of how many shadows we've seen (remember, there's 3 mutants per site)
            let current_freq = n_shadows_per_pos.entry(pos).or_insert(0);
            *current_freq += f;
        }
        else{
            // just fill with 0 if not there,
            n_shadows_per_pos.entry(pos).or_insert(0);
        }
    }
    n_shadows_per_pos
}


pub fn run(fastq_list: Vec<String>, whitelist_file: String, output_csv_file: String, topn:usize){

    // parse whitelist
    let whitelist = parse_whitelist_gz(whitelist_file);
    println!("Whitelist len {}", whitelist.len());

    
    let countmap = count_cb_filelist(fastq_list);
    println!("len of counter {}", countmap.len());


    // now the hard part: group by CB, look at all UMIs therein
    println!("calculating most common");
    let most_common: Vec<String> = top_n(&countmap, topn);
    println!("most common {:?}", most_common.len());

    // filter the most common for whitelist
    let most_common_filtered = most_common.into_iter().filter(|x| whitelist.contains(x)).collect::<Vec<String>>();


    // for each entry in most common, find potential shadows
    // add them to polars
    let mut polars_data: HashMap<String, Vec<i32>> = HashMap::new(); // building up the columns of the DataFrame: name->values
    let mut cellnames: Vec<String> = Vec::new();

    for mc in most_common_filtered{
        let mc2 = mc.clone();

        // finding shadoes per position
        let nshadows_per_position = find_shadows(mc, &countmap);

        for (position, n_shadows) in nshadows_per_position.iter(){
            // if *n_shadows > 0 {
            //     println!("shadows {:?}", nshadows_per_position);
            // }
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = countmap.get(&mc2).unwrap();
        polars_data.entry("total".into()).or_insert(vec![]).push(*total_counts);
        cellnames.push(mc2)
    }
    
    // to polars dataframe
    let df = DataFrame::new(
        polars_data.into_iter()
            .map(|(name, values)| Series::new(&format!("{name}"), values))
            .collect::<Vec<_>>()).unwrap();
    
    let df_cb = Series::new("CB", cellnames);

    let mut df_final = df.hstack(&[df_cb]).unwrap();
    println!("{:?}", df_final);

    // write to CSV
    write_to_csv(&mut df_final, output_csv_file);    

}