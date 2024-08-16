use counter::Counter;
use itertools::izip;
use crate::io::fastq_phred_iter;
use indicatif::{ProgressBar, ProgressStyle, };


#[cfg(test)]
#[test]
fn main(){
    use crate::test_files::TEST_FASTQ_R1;
    run(&vec![TEST_FASTQ_R1.to_string()],"/tmp/phred.csv".to_string())
}


pub fn run(fastq_files: &[String], output_csv_file:String){

    let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts

    let my_iter = fastq_phred_iter(fastq_files);

    let bar = ProgressBar::new_spinner();
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {pos} {per_sec}").unwrap()
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


    // let mut df = df!("phred" => phred_scores, "frequency" => freqs, "position" => positions).unwrap();
    // write_to_csv_polars(&mut df, output_csv_file);    

    write_to_csv_simple(phred_scores,positions, freqs, output_csv_file).unwrap();
}

pub fn write_to_csv_simple(phred_scores: Vec<String>, positions: Vec<u64>, freqs: Vec<u64>, output_csv_file: String) -> Result<(), csv::Error>{
    let mut wtr = csv::Writer::from_path(output_csv_file)?;

    wtr.write_record(["phred","frequency","position"])?;
    for (ph, pos, freq) in izip!(phred_scores, positions, freqs) {
        wtr.write_record(&[ph, pos.to_string(), freq.to_string()])?;
    }
    wtr.flush()?;
    Ok(())

}
// pub fn write_to_csv_polars(df_final: &mut DataFrame, output_csv_file: String){
//     let mut output_file: File = File::create(output_csv_file).unwrap();
//     CsvWriter::new(&mut output_file)
//         .include_header(true)
//         .finish(df_final)
//         .unwrap();    
// }
