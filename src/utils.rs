use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bgzf;
use std::collections::HashSet;
use polars::prelude::{CsvWriter, DataFrame, SerWriter};
use std::fs::File;
use core::hash::Hash;


pub fn fastq_iter_bare(fastq_list: &Vec<String>, line: usize) -> impl Iterator<Item=String> + '_ {
    // iterates over the concatenation of fastq files specified
    // each fastq entry is four lines, we only pick out a specific one using the
    // line arg:  
    // line==0 -> header
    // line==1 -> fastq seq
    // line==2 -> sep
    // line==3 -> phred

    let file_iterators = fastq_list.into_iter()
        .map(move |fname|{
            let decoder = bgzf::Reader::from_path(fname).unwrap();
            let reader = BufReader::new(decoder);
            let line_tmp = line.clone();
            let my_iter = reader.lines()
                .enumerate().filter( move |(line_id, _)| line_id % 4 == line_tmp)
                .map(|(_, res)| res)
                .filter_map(|line| line.ok()); //takes care of errors in file reading
            my_iter
        }
        );
    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);
    my_iter
}


pub fn fastq_iter(fastq_list: &Vec<String>) -> impl Iterator<Item=CbUmi> + '_ {
    // iterates the sequences of the the fast files
    let my_iter = fastq_iter_bare(fastq_list, 1).filter_map(|line| parse_r1_struct(line));
    my_iter
}

pub fn phred_iter(fastq_list: &Vec<String>) -> impl Iterator<Item=String> + '_ {
    // instead if yielding the sequence, this one yields the PHRED ASCII scores of the reads
    let my_iter = fastq_iter_bare(fastq_list, 3);
    my_iter
}


pub fn parse_whitelist_gz(fname: &String) -> HashSet<String>{
    // loading 10x CB whilelist from file
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines();//.take(10_000_000);
    let mut hset: HashSet<String> = HashSet::new();
    for l in my_iter{
        if let Ok(line) = l{
            hset.insert(line);
        }
    }
    hset
}


pub fn parse_r1_struct(seq: String) -> Option<CbUmi>{
    // TODO: size check
    // if seq.len() == 16 + 12{
        let cb = (&seq[0..16]).into();
        let umi = (&seq[16..28]).into();
        Some(CbUmi{cb, umi})
    // }
    // else
    // {
        // None
    // }
}

#[derive(Hash, PartialEq, Eq, Clone)]
pub struct CbUmi{
    pub cb: String,
    pub umi: String,
}

impl CbUmi {
    pub fn to_string(&self) -> String{
        format!("{}_{}", self.cb, self.umi)
    }
    pub fn from_string(s:&String) -> CbUmi{
        let mut split = s.split('_');
        let cb = split.next().expect(s);
        let umi = split.next().expect(s);
        CbUmi {
            cb: cb.to_string(), 
            umi: umi.to_string()
        }
    } 
}


#[derive(Hash, PartialEq, Eq, Clone)]
pub struct CbUmiGene{
    pub cb: String,
    pub umi: String,
    pub gene: String
}
impl CbUmiGene {
    pub fn to_string(&self) -> String{
        format!("{}_{}_{}", self.cb, self.umi, self.gene)
    }
    pub fn from_string(s:&String) -> CbUmiGene{
        let mut split = s.split('_');
        let cb = split.next().expect(s);
        let umi = split.next().expect(s);
        let gene = split.next().expect(s);
        CbUmiGene {
            cb: cb.to_string(), 
            umi: umi.to_string(),
            gene: gene.to_string()
        }
    } 
}


pub fn get_1bp_mutations(seq: &String, pos: usize) -> Vec<String>{
    // return the three 1BP mutations of a sequence at the given position"
    let mut mutation = Vec::with_capacity(3);
    for base in ["A", "C", "G", "T"]{
        let mut m = seq.clone();
        m.replace_range(pos..(pos+1), base);
        if m != *seq{
            mutation.push(m);
        }
    }
    // for base in ['A', 'C', 'G', 'T']{
        // let prefix = &seq[0..pos];
        // let postfix = &seq[(pos+1)..];
        // let m = format!("{prefix}{base}{postfix}");
        // if m != *seq{
        //     mutation.push(m);
        // }
    // }
    return mutation
}

pub fn all_mutations_for_cbumi(cb_umi: CbUmi) -> Vec<(CbUmi, usize)>{
    // we have to mutate both CB und UMI
    // lets turn it into a plain CBUMI string, mutate and turn back to CB_UMI
    let total_len = cb_umi.cb.len() + cb_umi.umi.len();

    let seq_plain = format!("{}{}", cb_umi.cb, cb_umi.umi);  // turn into single string

    let mut shadows_plain: Vec<(String, usize)> = Vec::with_capacity(total_len*3);
    for pos in 0..total_len{
        let mutations_at_pos = get_1bp_mutations(&seq_plain, pos).into_iter().map(|x| (x,pos));
        for matp in mutations_at_pos{
            shadows_plain.push(matp);
        }
    }
    assert_eq!(shadows_plain.len(), total_len*3);

    // this looks like [("AAA", 0), ("BAA", 0), ... ("AAB", 2)]
    
    // convert back to Cb Umi
    let shadows = shadows_plain.iter() // convert back to CbUmi
        .map(|(cbumi_str, pos)| {
            let cbumi = CbUmi{ 
                cb: (&cbumi_str[0..16]).to_string(), 
                umi: (&cbumi_str[16..28]).to_string()
            };
            (cbumi, *pos)
        }
    ).collect();
    shadows
}

pub fn set_comparison<T>(set_a: &HashSet<T>,set_b: &HashSet<T>)
    where T: Eq + Hash
{
    // warning: this instantiates all the sets
    // instead, just iterate over them
    // let intersect = set_a.intersection(set_b).collect::<Vec<&T>>();
    // let AminusB = set_a.difference(set_b).collect::<Vec<&T>>();
    // let BminusA = set_b.difference(set_a).collect::<Vec<&T>>();
    // println!("|A&B| {} |A-B| {}  |B-A|{}", intersect.len(), AminusB.len(), BminusA.len());

    // let n_intersect = set_a.intersection(set_b).fold(0, |accum, _item| accum+1);
    // let n_A_minus_B = set_a.difference(set_b).fold(0, |accum, _item| accum+1);
    // let n_BminusA = set_b.difference(set_a).fold(0, |accum, _item| accum+1);

    let n_intersect = set_a.intersection(set_b).count();
    let n_a_minus_b = set_a.difference(set_b).count();
    let n_bminus_a = set_b.difference(set_a).count();

    println!("|A&B| {} |A-B| {}  |B-A|{}", n_intersect, n_a_minus_b, n_bminus_a);

}


pub fn write_to_csv(df_final: &mut DataFrame, output_csv_file: String){
    let mut output_file: File = File::create(output_csv_file).unwrap();
    CsvWriter::new(&mut output_file)
        .has_header(true)
        .finish(df_final)
        .unwrap();    
}


pub fn my_hamming(a: &String, b: &String) -> isize {
    // hamming distance for two strings of the same size
    assert_eq!(a.len(), b.len());
    let mut counter: isize = 0;
    // for (c1, c2) in  std::iter::zip((*a).chars(), (*b).chars()){  // todo: change to bytes, might be faster
    for (c1, c2) in  std::iter::zip((*a).bytes(), (*b).bytes()){  // todo: change to bytes, might be faster
        if c1 != c2{
            counter +=1 ;
        }
    };
    counter
}


pub fn sequence_composition(s: &String) -> (u32, u32,u32,u32){
    let mut counter_a = 0;
    let mut counter_c = 0;
    let mut counter_g = 0;
    let mut counter_t = 0;

    for c in s.chars(){
        match c{
            'A' => counter_a+=1,
            'C' => counter_c+=1,
            'G' => counter_g+=1,
            'T' => counter_t+=1,
            _ => panic!("Unknown char {}", c)
        }
    }

    return (counter_a, counter_c, counter_g, counter_t)

}

#[cfg(test)]
mod tests {
    // use std::io::W/rite;
    // use crate::bus::{BusRecord, BusHeader, CellIterator, BusIteratorBuffered, BusWriter};
    use crate::utils::*;

    #[test]
    fn test_my_hamming(){
        assert_eq!(my_hamming(&("AAAA".to_string()), &("AAAT".to_string())), 1);
        assert_eq!(my_hamming(&("AAAA".to_string()), &("TTTT".to_string())), 4);

    }
    #[test]
    fn test_get_1bp_mutations(){
        assert_eq!(get_1bp_mutations(&"AAAA".to_string(), 2), vec!["AACA".to_string(), "AAGA".to_string(), "AATA".to_string()]);
        assert_eq!(get_1bp_mutations(&"TTTT".to_string(), 1), vec!["TATT".to_string(), "TCTT".to_string(), "TGTT".to_string()]);
    }
}