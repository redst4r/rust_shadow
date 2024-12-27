use core::str;
use std::io::BufReader;

use counter::Counter;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fastq::{parse_path, Record};
use rust_htslib::bgzf;
use rustfastq::{io::fastq_list_iter, test_files::TEST_FASTQ_R1};
// use noodles_fastq;
use noodles;

fn plain_and_compressed_iterator_speed(c: &mut Criterion){

    fn do_mine(fq: Vec<String>) -> usize{
        let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts

        for entry in fastq_list_iter(&fq).take(100_000) {
            
            let phred_string = entry.phred;
            for (position, phred_score) in phred_string.chars().enumerate(){
                let counter = phred_counter.entry((phred_score, position)).or_insert(0);
                *counter += 1;      
            }
        }
        phred_counter.len()
    }

    let fastq_list = vec![TEST_FASTQ_R1.to_string()];


    c.bench_function("my iterator",
     |b| b.iter(|| {
        do_mine(fastq_list.clone())
     }
    ));
    
    fn do_noodles(fq: Vec<String>) -> usize{
        let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts
        
        let decoder = bgzf::Reader::from_path(fq[0].clone()).unwrap();
        let inner = BufReader::with_capacity(800 * 1024, decoder);
        let mut reader = noodles::fastq::io::Reader::new(inner);
        for entry in reader.records().take(100_000) {
            
            let fq_entry = entry.unwrap();

            let phred_bytes = fq_entry.quality_scores();
            let phred_string = str::from_utf8(phred_bytes).unwrap();

            for (position, phred_score) in phred_string.chars().enumerate(){
                let counter = phred_counter.entry((phred_score, position)).or_insert(0);
                *counter += 1;      
            }
        }
        phred_counter.len()
    }

    c.bench_function("noodles iterator",
    |b| b.iter(|| {
        do_noodles(fastq_list.clone())
    }
   ));

    fn do_fastq(fname: &str) -> usize {
        let h = parse_path(Some(fname), |mut parser| {
            let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts
            let mut counter = 0;
            parser.each(|record| {
                let phred_bytes = record.qual();
                let phred_string = str::from_utf8(phred_bytes).unwrap();
                for (position, phred_score) in phred_string.chars().enumerate(){
                    let counter = phred_counter.entry((phred_score, position)).or_insert(0);
                    *counter += 1;      
                }
                counter += 1;
                if counter < 100000 {
                    true
                } else {
                    false
                }
            }).expect("Invalid fastq file");
            phred_counter
        }).expect("some error");
        h.len()
    }
    c.bench_function("fastq iterator",
    |b| b.iter(|| {
        do_fastq(&fastq_list[0].clone())
    }
   ));


   fn do_fastq_multicore(fname: &str) -> usize {
    let h = parse_path(Some(fname), |mut parser| {
        let mut counter = 0;
        let nthreads = 4;
        let results: Vec<_> = parser.parallel_each(nthreads, |record_sets| {
            let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts
            let mut thread_total = 0;
            for record_set in record_sets {
                for record in record_set.iter() {
                    let phred_bytes = record.qual();
                    let phred_string = str::from_utf8(phred_bytes).unwrap();
                    for (position, phred_score) in phred_string.chars().enumerate(){
                        let counter = phred_counter.entry((phred_score, position)).or_insert(0);
                        *counter += 1;      
                    }                    
                    thread_total += 1;
                }
                if thread_total > 100000 {
                    break
                }
            }
            phred_counter
        }).expect("Invalid fastq file");

        let mut phred_counter: Counter<(char, usize), u64> = Counter::new();  // phred, position -> #counts
        for h in results {
            phred_counter.extend(h);
        }
        phred_counter 
    }).expect("some error");
    h.len()
    }

    c.bench_function("fastq paralell iterator",
    |b| b.iter(|| {
        do_fastq_multicore(&fastq_list[0].clone())
    }
   ));

}

criterion_group!(benches, plain_and_compressed_iterator_speed);

criterion_main!(benches);