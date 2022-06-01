use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bgzf;
use std::collections::{HashSet};
use polars::prelude::{CsvWriter, DataFrame, SerWriter};
use std::fs::File;
use core::hash::Hash;

pub fn fastq_iter(fastq_list: &Vec<String>) -> impl Iterator<Item=CbUmi> + '_ {
    let file_iterators = fastq_list.into_iter()
        .map(|fname|{
            let decoder = bgzf::Reader::from_path(fname).unwrap();
            let reader = BufReader::new(decoder);
            let my_iter = reader.lines()
                .enumerate().filter(|x| x.0 % 4 == 1)
                .map(|x| x.1)
                .filter_map(|line| line.ok()) //takes care of errors in file reading
                .filter_map(|line| parse_r1_struct(line))
                // .take(TOTAL_READS)
                ;
            my_iter
        }
        );
    // chaining, flatmapping all the iterators into a single one
    let my_iter = file_iterators.flat_map(|x| x);
    my_iter
}

pub fn parse_whitelist_gz(fname: String) -> HashSet<String>{
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

// pub fn parse_r1(seq: String) -> Option<(String, String)>{
//     // TODO: size check
//     // if seq.len() == 16 + 12{
//         let cb = (&seq[0..16]).into();
//         let umi = (&seq[16..28]).into();
//         Some((cb, umi))
//     // }
//     // else
//     // {
//         // None
//     // }
// }


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

pub fn seq_to_int(seq: String) -> u64{

    assert!(seq.len() <= 32); // cant handle longer sequences in a single 64bit integer!
    let s: String = seq
    .chars()
    .map(|x| match x {
        'A' => '0',
        'C' => '1',
        'G' => '2',
        'T' => '3',
        _ => panic!("unkown seq character"),
    })
    .collect();
    let int_seq = u64::from_str_radix(&s, 4).unwrap();
    int_seq
}


pub fn int_to_seq(i: u64, seq_len:u64) -> String{

    let mut q = i;
    let mut result: Vec<u64> = Vec::new();
    while q>= 4{
        let quotient = q / 4;
        let remainder = q % 4;
        result.push(remainder);
        q = quotient;
    }
    result.push(q);

    while result.len() < seq_len as usize{
        result.push(0);

    }

    result.reverse();

    let s: String = result.iter()
    .map(|x| match x {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => panic!("unkown seq character"),
    })
    .collect::<String>();

    // println!("{:?}", s);
    s
}   


pub fn my_hamming(a: &String, b: &String) -> isize {
    // hamming distance for two strings of the same size
    assert_eq!(a.len(), b.len());
    let mut counter: isize = 0;
    for (c1, c2) in  std::iter::zip((*a).chars(), (*b).chars()){
        if c1 != c2{
            counter +=1 ;
        }
    };
    counter
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
    #[test]
    fn encode_seq(){
        assert_eq!(seq_to_int("A".to_string()), 0);
        assert_eq!(seq_to_int("C".to_string()), 1);
        assert_eq!(seq_to_int("G".to_string()), 2);
        assert_eq!(seq_to_int("T".to_string()), 3);

        assert_eq!(seq_to_int("GCCA".to_string()), 148);

    }

    #[test]
    fn decode_seq(){
        let x = vec![0,1,2,3,4,5,6,7,8];
        // for i in x{
        //     println!("DECDODED {}:{:?}", i, int_to_seq(i, 5));
        // }


        //  base order
        assert_eq!(int_to_seq(0, 1), "A");
        assert_eq!(int_to_seq(1, 1), "C");
        assert_eq!(int_to_seq(2, 1), "G");
        assert_eq!(int_to_seq(3, 1), "T");
    // 
        // # padding leading A's
        assert_eq!(int_to_seq(0, 3), "AAA");
        assert_eq!(int_to_seq(1, 3), "AAC");
        assert_eq!(int_to_seq(2, 3), "AAG");
        assert_eq!(int_to_seq(3, 3), "AAT");
    // 
        assert_eq!(int_to_seq(4, 2), "CA");
        assert_eq!(int_to_seq(5, 2), "CC");
        assert_eq!(int_to_seq(6, 2), "CG");
        assert_eq!(int_to_seq(7, 2), "CT");
    // 
        assert_eq!(int_to_seq(148, 4), "GCCA");
    // 
        // # make sure to raise an error when the decoded string is actually longer
        // # then requested (since its probably a bug in the code calling _decode_int_to_ACGT)
        // with pytest.raises(AssertionError):
            // busio._decode_int_to_ACGT(148, seq_len=2)
        // with pytest.raises(AssertionError):
            // busio._decode_int_to_ACGT(-1, seq_len=1)
    }

}