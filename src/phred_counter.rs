use counter::Counter;
use crate::utils::{write_to_csv, };
use polars::prelude::*;
use crate::io::phred_iter;
use indicatif::{ProgressBar, ProgressStyle, };


#[cfg(test)]
#[test]
fn main(){
    use crate::test_files::TEST_FASTQ_R1;
    run(&vec![TEST_FASTQ_R1.to_string()],
    "/tmp/phred.csv".to_string(),
)
}


pub fn run(fastq_files: &Vec<String>, output_csv_file:String){

    let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts

    let my_iter = phred_iter(fastq_files);

    let bar = ProgressBar::new_spinner();
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {pos} {per_sec}")
        .progress_chars("##-"));

    for phred_string in my_iter{  //.take(1_000_0)
        for (position, phred_score) in phred_string.chars().enumerate(){
            let counter = phred_counter.entry((phred_score, position)).or_insert(0);
            *counter += 1;      
        }
        bar.inc(1);
    }
    bar.finish();

    // unwrap the whole thing int a dataframe with three cols: phred-char, position, freq
    let mut phred_scores: Vec<String> = Vec::new();
    let mut positions: Vec<u64> = Vec::new();
    let mut freqs: Vec<u64> = Vec::new();
    for ((phred_char, pos ), freq) in phred_counter.into_iter(){
        phred_scores.push(phred_char.to_string());
        positions.push(pos as u64);
        freqs.push(freq);
    }

    // let df_cb = Series::new("PHRED", phred_scores);
    // let series_phred = Series::new("phred", phred_scores);
    // let series_freq = Series::new("frequency", freqs);
    // let series_pos = Series::new("position", positions);


    let mut df = df!("phred" => phred_scores, "frequency" => freqs, "position" => positions).unwrap();
    
    write_to_csv(&mut df, output_csv_file);    

}
