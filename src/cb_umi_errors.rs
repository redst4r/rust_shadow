use std::collections::HashMap;
use counter::Counter;
use bktree::BkTree;
use polars::prelude::{DataFrame, NamedFrom, Series};
use crate::utils::{parse_whitelist_gz, write_to_csv, all_mutations_for_cbumi, CbUmi, fastq_iter, my_hamming};

pub fn count_cb_filelist(fname_list: &Vec<String>) -> Counter<CbUmi, u32> {
    // coutns the CB/UMI pairs in the fastq

    // reading the fastq.gz
    let my_iter = fastq_iter(fname_list);

    // parsing the lines, counting
    let mut countermap: Counter<CbUmi, u32> = Counter::new();

    for (i, cbumi) in my_iter.enumerate(){
        let counter = countermap.entry(cbumi).or_insert(0);
        *counter += 1;
        
        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/1_000_000)
        }
    }
    countermap
}

pub fn top_n(counter: &Counter<CbUmi, u32>, n: usize) -> Vec<CbUmi>{

    // gets the Top_n elements from the counter, making sure that
    // they are not shadows of other frequent elements
    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    let mut top2: Vec<CbUmi> = Vec::new();

    let mut c = 0;
    for (cbumi, _freq) in counter.most_common(){

        let matches = bk.find(cbumi.to_string(), 1);
        if matches.len() == 0{
            top2.push(cbumi.clone());  // add it to the topN list
            c += 1;
        }
        bk.insert(cbumi.to_string());

        if c >= n{
            break
        }
        if c % 1_000 == 0{
            println!("Iteration {c} of {n}");
        }          
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
    
    let countmap = count_cb_filelist(fastq_list);
    // transform into shadow counter
    println!("len of counter {}", countmap.len());

    // filter for whitelist only entries
    // println!("Filtering for whilelist");

    // let mut keys_to_remove: Vec<CbUmi> = Vec::new();
    // for (cbumi, _count) in countmap.iter(){
    //     if !whitelist.contains(&cbumi.cb){
    //         keys_to_remove.push(cbumi.clone());
    //     }
    // }

    // for k in keys_to_remove{
    //     countmap.remove(&k);
    // }

    // println!("len of filtered counter {}", countmap.len());
    // unstable:
    // countmap.drain_filter(|k, v| {
    //     k.0 == k.1
    // });

    // debugging: save the counter
    println!("saving counter to file /tmp/test.csv");
    let mut cbs: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();
    let mut freqs: Vec<u32> = Vec::new();
    for (seq, freq) in countmap.iter(){
        cbs.push(seq.cb.clone());
        umis.push(seq.umi.clone());
        freqs.push(*freq);
    }
    let df_cb = Series::new("CB", cbs);
    let df_umi = Series::new("UMI", umis);
    let df_freq = Series::new("frequency", freqs);
    let mut df_ = DataFrame::new(vec![df_cb,df_umi, df_freq]).unwrap();
    write_to_csv(&mut df_, "/tmp/exact_counter_cbumi.csv".to_string());


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

        // if !whitelist.contains(&mc.cb){
        //     continue
        // }

        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &countmap);
        for (position, n_shadows) in nshadows_per_position.iter(){
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = countmap.get(&mc2).unwrap();
        polars_data.entry("n_real".into()).or_insert(vec![]).push(*total_counts);

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
