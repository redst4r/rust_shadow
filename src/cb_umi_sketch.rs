// tracking the freqs of all CB-UMI is too much memory-wise
// lets use some probabilistic data structs
//
// We only store frequently seen CB/UMI (those are the candidates for "real" reads, not shadows)
// 1. either just looks for top10k CB/UMIs: problem those might not have an shadows
// 2. just keep every CB/UMI that we see more than once as a candidate for shadows

// For c in Approx-candidates
//
//      For a in possible_shadows(c)
//          count_freq(a)  // this includes iterating over the fastq again!
// 
//
use counter::Counter;
use std::collections::HashMap;
use crate::utils::{write_to_csv, parse_r1_struct, parse_whitelist_gz, all_mutations_for_cbumi, fastq_iter, CbUmi} ;
use crate::sketching::GreaterThan1Bloom;
use crate::cb_umi_errors::{top_n, find_shadows};
use polars::prelude::{DataFrame, NamedFrom, Series};
use streaming_algorithms::Top;
// use indicatif::ProgressIterator;
// use indicatif::ProgressBar;

const TOTAL_READS: usize = 50_000_000;


pub fn run_topN(fastq_list: &Vec<String>, whitelist_file: String, output_csv_file: String, topn:usize){
    let whitelist = parse_whitelist_gz(whitelist_file);
    println!("Whitelist len {}", whitelist.len());

    let my_iter = fastq_iter(fastq_list);
    //--------------------------
    // first pass over data
    // keepign track of the topN elements in  the approximate counter
    // however, only add XB/UMI that are whitelisted
    let tol = 1e-8;
    let prob = 0.000001;
    let mut ccc:Top<String, u32> = Top::new(topn, prob, tol, {});    

    // let bar = ProgressBar::new_spinner();
    for (i, cbumi) in my_iter.enumerate(){

        if !whitelist.contains(&cbumi.cb){
            continue;
        }

        let cm_umi_str = cbumi.to_string();  //convert to CB_UMI string for hashing
        ccc.push(cm_umi_str, &1);
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000);
        }
        // bar.inc(1)
    }
    // bar.finish();

    println!("Done with first pass");

    // now we have a candidate list, all are valid CB accodinf to the whitelist
    // get their ACTUAL frequences and the frequencies of their shadows
    // first, lets build a list of all candidates and their shadows
    let mut candidates_and_shadows: Counter<CbUmi, u32> = Counter::new(); //::with_capacity( TOPN * 28 * 4);  // around 1M entries for N=10k
    
    for (seq, _approx_freq) in ccc.iter(){
        let seq_cbumi = CbUmi::from_string(seq);
        
        // insert the sequence itself 
        // let qqq = seq.clone(); //qqq is of the form String: CB_UMI
        candidates_and_shadows.insert( seq_cbumi, 0);
        
        // insert all its potential shadows across positions
        let seq_cbumi = CbUmi::from_string(seq);
        let shadows = all_mutations_for_cbumi(seq_cbumi);
        for (s, _mutated_pos) in shadows{
            candidates_and_shadows.insert( s, 0);
        }
    }

    // now go thorugh the fastqs again, recoding the TRUE frequencies of those items
    let my_iter = fastq_iter(fastq_list);

    println!("Second pass");

    for (i, cbumi) in my_iter.enumerate(){
        if candidates_and_shadows.contains_key(&cbumi){
            let c = candidates_and_shadows.entry(cbumi).or_insert(0);  // the INSERT SHOULD NEVER HAPPEN
            // note that this will contain a 0 counter if we see the first read
            // this is not from the or_insert!! but from out init
            *c+=1;
        }
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000);
        }               
    }
    println!("{}", candidates_and_shadows.len());

    // now we have all the actual counts!
    // identify the ACTUAL REAL reads () removing possible frequent shadows
    println!("convertin ");

    // collect the topN CUs and their ACTUAL counts (ccc only stores the approximate counts)
    // this should be done by a toolz.valmap!!
    let mut countmap: Counter<CbUmi, u32> = Counter::new();
    
    let mut i = 0;
    // this is the actual 10k topN
    for (seq, _approx_freq) in ccc.iter(){
        let cbumi = CbUmi::from_string(seq);
        let real_freq = candidates_and_shadows.get(&cbumi).unwrap();

        // let cbumi = CbUmi::from_string(seq);
        // let cbumi = parse_r1_struct((*seq).clone()).unwrap();
        if whitelist.contains(&cbumi.cb){
            countmap.insert(cbumi, *real_freq);
        }
        else{
            panic!("this shouldnt happen, its already filtered")
        }
        i+=1;
    }
    assert_eq!(i, topn);


    println!("calculating most common; {}", countmap.len());
    let most_common: Vec<CbUmi> = top_n(&countmap, topn);
    println!("most common {:?}", most_common.len());


    // for each entry in most common, find potential shadows
    // add them to polars

    let mut polars_data: HashMap<String, Vec<u32>> = HashMap::new();
    let mut cellnames: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();

    for mc in most_common{
        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &candidates_and_shadows);

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



pub fn run_GB1(fastq_list: &Vec<String>, whitelist_file: String, output_csv_file: String){
    
    let whitelist = parse_whitelist_gz(whitelist_file);
    println!("Whitelist len {}", whitelist.len());

    let my_iter = fastq_iter(fastq_list);

    // first pass over data, storing all elements #>1
    let tol = 1e-9;
    let prob = 1e-4;
    let mut GB1 :GreaterThan1Bloom<> = GreaterThan1Bloom::new(prob, tol);

    for (i, cbumi) in my_iter.enumerate(){
        let cm_umi_str = cbumi.to_string();
        GB1.add_item(&cm_umi_str);
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000);
            GB1.status();
        }        
    }
}
