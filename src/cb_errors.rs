use std::collections::HashMap;
use counter::Counter;
use bktree::BkTree;
use crate::{utils::{get_1bp_mutations, write_to_csv, my_hamming}, io::fastq_seq_iter};
use polars::prelude::{DataFrame, NamedFrom, Series};
use crate::io::{parse_whitelist_gz};

pub fn count_cb_filelist(fname_list: &[String]) -> Counter<String, i32> {
    // coutns the CB/UMI pairs in the fastqs
    // reading the fastq.gz
    let my_iter = fastq_seq_iter(fname_list);

    // parsing the lines, counting
    let mut countermap: Counter<String, i32> = Counter::new();

    for (i, cbumi) in my_iter.enumerate(){
        let counter = countermap.entry(cbumi.cb).or_insert(0);
        *counter += 1;

        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}

pub fn top_n(counter: &Counter<String, i32>, n: usize) -> Vec<String>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements

    // ======================================================================
    // another problem: we only add to the BKTree if no related element is present
    // However, consider three items A:1000,B:100,C:10 in descending frequency and
    //  A-1->B-1->C  (A and B are one distance, B and C are one distance, A,C are two distance)
    // - A gets added
    // - we encounter B and skip it(since its a shadow of A, which is more frequent)
    // - we encounter C (which is a shadow of B). Since B was never added to the BKTree
    //   we would consider C a true molecule (even though a MORE FREQUENT seq within distance 1 exists!!)
    // 
    // Hence we really have to add the shadows to the BKTree itself.
    // 
    // This changes the criterion for the final list to:
    // - any item A in this list DOES NOT have a 1-distance neigbour B that is more frequent than A

    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    let mut top2: Vec<String> = Vec::new();

    let mut c = 0;
    for (cb, _freq) in counter.most_common(){
        let cb2 = cb.clone();  // TODO stupid cloing to be able to insert
        if bk.find(cb, 1).is_empty() {
            // let cb_umi = format!("{cb}_{umi}");
            top2.push(cb2.clone());  // add it to the topN list
            c += 1;
        }
        bk.insert(cb2); // add it to the BKtree anyways, see function header for explain

        if c >= n{
            break
        }
        if c % 1_000 == 0{
            println!("Iteration {c} of {n}");
        }          
    }
    assert_eq!(top2.len(), n);
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


pub fn run(fastq_list: &[String], whitelist_file: String, output_csv_file: String, topn:usize){

    // parse whitelist
    let whitelist = parse_whitelist_gz(&whitelist_file);
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
        polars_data.entry("n_real".into()).or_insert(vec![]).push(*total_counts);
        cellnames.push(mc2)
    }
    
    // to polars dataframe
    let df = DataFrame::new(
        polars_data.into_iter()
            .map(|(name, values)| Series::new(&name, values))
            .collect::<Vec<_>>()).unwrap();
    
    let df_cb = Series::new("CB", cellnames);

    let mut df_final = df.hstack(&[df_cb]).unwrap();
    println!("{:?}", df_final);

    // write to CSV
    write_to_csv(&mut df_final, output_csv_file);    

}