/// esimating the error fromth TSO
/// 
use std::collections::HashMap;
use counter::Counter;
use bktree::BkTree;
use crate::utils::{write_to_csv, my_hamming, CbUmi};
use crate::io::{fastq_iter, fastq_seq_iter};
use polars::prelude::{DataFrame, NamedFrom, Series};
use crate::cb_umi_errors::find_shadows;
use indicatif::{ProgressBar, ProgressStyle};


#[test]
fn testing(){
    use crate::test_files::TEST_FASTQ_R1;
    let s = vec![TEST_FASTQ_R1.to_string()];
    run(&s, "/tmp/tso.csv".to_string());
}

pub fn run(fast_files: &Vec<String>, output_csv_file:String){

    // check anything that is close to the TSO sequence, or shifted by a few base pairs
    let tso1 = "AAGCAGTGGTATCAAC_GCAGAGTACATG".to_string();  //note the _ sepaation of CB and UMI
    let tso2 = "AGCAGTGGTATCAACG_CAGAGTACATGG".to_string();
    let tso3 = "GCAGTGGTATCAACGC_AGAGTACATGGG".to_string();
    let tso4 = "CAGTGGTATCAACGCA_GAGTACATGGGG".to_string();
    let tso_list = vec![tso1, tso2, tso3, tso4];
    let my_iter = fastq_seq_iter(fast_files);
    let mut counter: Counter<CbUmi, u32> = Counter::new();

    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    for seq in tso_list.iter(){
        bk.insert(seq.clone())
    }

    // let top_x = 1_000_000;
    let bar = ProgressBar::new_spinner();
    // let bar = ProgressBar::new(top_x);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {spinner} {per_sec}")
        .progress_chars("##-"));


    // go over all records, see if anything is close to the TSO
    for seq in my_iter{  // this will do the same _ separation
    // for seq in my_iter.take(10000000 as usize){  // this will do the same _ separation
        bar.inc(1);

        // using a BKtree
        if true{
            let hits = bk.find(seq.to_string(), 1);
            if hits.len() > 0{
                // we found something thats close to the TSO
                let counter = counter.entry(seq).or_insert(0);
                *counter += 1;
            }
        }
        else{
            let max_distance = tso_list.iter().map(|tso| my_hamming(&tso, &seq.to_string())).min().unwrap();
            if max_distance <=1{

                let counter = counter.entry(seq).or_insert(0);
                *counter += 1;
            }
        }
    }
    bar.finish();
    // having found all TSO like sequences, lets detmine shadows
    
    // for each entry in TSO, find potential shadows
    // add them to polars
    let mut polars_data: HashMap<String, Vec<u32>> = HashMap::new();
    let mut cellnames: Vec<String> = Vec::new();
    let mut umis: Vec<String> = Vec::new();

    for mc in tso_list.into_iter().map(|c| CbUmi::from_string(&c)){

        if !counter.contains_key(&mc){
            continue
        }

        let mc2 = mc.clone();
        let nshadows_per_position = find_shadows(mc, &counter);
        for (position, n_shadows) in nshadows_per_position.iter(){
            polars_data.entry(format!("position_{position}")).or_insert(vec![]).push(*n_shadows)
        }
        let total_counts = counter.get(&mc2).unwrap();
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