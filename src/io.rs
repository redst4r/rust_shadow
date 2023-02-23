use std::io::BufReader;
use std::io::BufRead;

use rust_htslib::bgzf;
use std::collections::HashSet;
use crate::utils::CbUmi;
use crate::utils::parse_r1_struct;

pub fn fastq_iter(fastqname: &str, line: usize) -> impl Iterator<Item=String> + '_  {
    // iterates over the fastq file
    // each fastq entry is four lines, we only pick out a specific one using the
    // line arg:  
    // line==0 -> header
    // line==1 -> fastq seq
    // line==2 -> sep
    // line==3 -> phred    
    assert!(line < 4);
    // assert!(line >= 0);

    let decoder = bgzf::Reader::from_path(fastqname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines()
        .skip(line)
        .step_by(4)
        // .enumerate().filter( move |(line_id, _)| line_id % 4 == line_tmp)
        .filter_map(|line| line.ok()); //takes care of errors in file reading
    my_iter
}

pub fn fastq_list_iter(fastq_list: &[String], line: usize) -> impl Iterator<Item=String> + '_ {
    // iterates over the concatenation of fastq files specified
    // each fastq entry is four lines, we only pick out a specific one using the
    // line arg:  
    // line==0 -> header
    // line==1 -> fastq seq
    // line==2 -> sep
    // line==3 -> phred
    assert!(line < 4);
    // assert!(line >= 0);
    let my_iter = fastq_list.iter()
        .flat_map(move |fname| fastq_iter(fname, line));

        // .map(move |fname| fastq_iter(fname, line))
        // chaining, flatmapping all the iterators into a single one
        // .flatten();     
    my_iter     
}

#[test]
fn test_fastq(){
    let fastq_entry1 = "@some read id
AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
";
    let fastq_entry2 = "@another read id
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
";
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;

    let fastqname = "/tmp/foo.fastq";
    let f = File::create(fastqname).expect("Unable to create file");
    let mut f = BufWriter::new(f);

    f.write_all(fastq_entry1.as_bytes()).expect("Unable to write data");
    f.write_all(fastq_entry2.as_bytes()).expect("Unable to write data");
    f.flush().unwrap();

    let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()],0).collect();
    assert_eq!(lines, vec!["@some read id", "@another read id"]);

    let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()],1).collect();
    assert_eq!(lines, vec!["AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"]);

    let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()],3).collect();
    assert_eq!(lines, vec!["FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"])    ;
}


pub fn fastq_seq_iter(fastq_list: &[String]) -> impl Iterator<Item=CbUmi> + '_ {
    // iterates the sequences of the the fast files
    fastq_list_iter(fastq_list, 1).filter_map(parse_r1_struct)
}

pub fn fastq_phred_iter (fastq_list: &[String]) -> impl Iterator<Item=String> + '_ {
    // instead if yielding the sequence, this one yields the PHRED ASCII scores of the reads
    fastq_list_iter(fastq_list, 3)
}

pub fn parse_whitelist_gz(fname: &String) -> HashSet<String>{
    // loading 10x CB whilelist from file
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines();//.take(10_000_000);
    let mut hset: HashSet<String> = HashSet::new();
    for line in my_iter.flatten(){  // flatten filters out the Error elements
        hset.insert(line);
    }
    hset
}