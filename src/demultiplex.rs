//! Code to demultiplex (based on I1, I2)
//! 

use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter, Write}, path::Path};

use itertools::{izip, Itertools};
use rust_htslib::bgzf::{self, CompressionLevel};

use crate::utils::get_spinner;

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub  struct DualIndex(pub String, pub String);

// to denote a single sample by some prefix
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub  struct Samplename(pub String);

pub struct Samplesheet {
    sheet: HashMap<DualIndex, Samplename>, 
    empty_sample: Samplename, // a special samplename indicating anything not matchgin the indices
}

impl Samplesheet {

    pub fn new(sheet: HashMap<DualIndex, Samplename>) -> Self {
        Self { 
            sheet , 
            empty_sample: Samplename("Undetermined".to_owned())
        }
    }

    /// Constructs the Samplesheet from a csv table
    /// with three columns i7 index, i5 index, samplename
    pub fn from_csv(file: &Path) -> Self {
        let sheet = samplesheet_to_hashmap_2(file);
        Self::new(sheet)
    }

    /// creates the FastQ writers for the sample sheet,
    /// i.e. each sample has a writer for R1 and R2
    fn create_writers(&self, outdir: &Path, undetermined_prefix: &str) -> HashMap<Samplename, (Box<dyn Write>, Box<dyn Write>)> {
        let prefixes = self.sheet.values().unique().collect_vec();
        let mut writers: HashMap<Samplename, (Box<dyn Write>, Box<dyn Write>)> = HashMap::new();
    
        fn get_encoder(fname: &str) -> BufWriter<rust_htslib::bgzf::Writer> {
            let encoder = BufWriter::new(
                bgzf::Writer::from_path_with_level(fname, CompressionLevel::Fastest).unwrap()
            );
            encoder
        }

        // create a writer for each sample (which can have multiple sample-indices)
        for sname in prefixes {
            let fname_r1 = format!("{}/{}.R1.fq.gz", outdir.to_str().unwrap(), sname.0);
            let fname_r2 = format!("{}/{}.R2.fq.gz", outdir.to_str().unwrap(), sname.0);
            let encoder_r1 = get_encoder(&fname_r1);
            let encoder_r2 = get_encoder(&fname_r2);
            writers.insert(sname.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));       
        }

        // add the writer for unassigned
        let fname_r1 = format!("{}/{}.R1.fq.gz", outdir.to_str().unwrap(), undetermined_prefix);
        let fname_r2 = format!("{}/{}.R2.fq.gz", outdir.to_str().unwrap(), undetermined_prefix);
        let encoder_r1 = get_encoder(&fname_r1);
        let encoder_r2 = get_encoder(&fname_r2);
        let empty_sample = Samplename("Undetermined".to_owned());
        writers.insert(empty_sample.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));        
        writers
    }


    pub fn get_samplename_from_index(&self, dual_ix: DualIndex) -> &Samplename {
        let samplename = self.sheet.get(&dual_ix).unwrap_or(&self.empty_sample,);
        samplename
    }

}


pub  fn demux_dual_index_2(samplesheet: Samplesheet, undetermined_prefix: String, i1_list: Vec<String>, i2_list: Vec<String>, r1_list: Vec<String>, r2_list: Vec<String>, outfolder: &Path) {

    // let prefixes = sample_indices_fnames.values().unique().collect_vec();
    // let mut writers: HashMap<Samplename, (Box<dyn Write>, Box<dyn Write>)> = HashMap::new();

    // // create a writer for each sample (which can have multiple sample-indices)
    // for sname in prefixes {
    //     let fname_r1 = format!("{}/{}.R1.fq.gz", outfolder, sname.0);
    //     let fname_r2 = format!("{}/{}.R2.fq.gz", outfolder, sname.0);

    //     let encoder_r1 = BufWriter::new(
    //         bgzf::Writer::from_path_with_level(fname_r1, CompressionLevel::Fastest).unwrap()
    //     );
    //     let encoder_r2 = BufWriter::new(
    //         bgzf::Writer::from_path_with_level(fname_r2, CompressionLevel::Fastest).unwrap()
    //     );
    //     writers.insert(sname.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));       
    // }
    // // add the writer for unassigned
    // let fname_r1 = format!("{}/{}.R1.fq.gz", outfolder, undetermined_prefix);
    // let fname_r2 = format!("{}/{}.R2.fq.gz", outfolder, undetermined_prefix);
    // let encoder_r1 = BufWriter::new(
    //     bgzf::Writer::from_path_with_level(fname_r1, CompressionLevel::Fastest).unwrap()
    // );
    // let encoder_r2 = BufWriter::new(
    //     bgzf::Writer::from_path_with_level(fname_r2, CompressionLevel::Fastest).unwrap()
    // );    
    // let empty_sample = Samplename("Undetermined".to_owned());

    // writers.insert(empty_sample.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));

    
    let mut writers = samplesheet.create_writers(outfolder, &undetermined_prefix);

    let i1_iter = crate::io::fastq_list_iter(&i1_list);
    let i2_iter = crate::io::fastq_list_iter(&i2_list);
    let r1_iter = crate::io::fastq_list_iter(&r1_list);
    let r2_iter = crate::io::fastq_list_iter(&r2_list);

    let pbar = get_spinner();
    for (counter, (i1, i2, r1, r2)) in izip!(i1_iter, i2_iter, r1_iter, r2_iter).enumerate() {
        let key: DualIndex = DualIndex(i1.seq, i2.seq);

        // let samplename = sample_indices_fnames.get(&key).unwrap_or(&empty_sample,);
        let samplename = samplesheet.get_samplename_from_index(key);
        let (writer_r1, writer_r2) = writers.get_mut(samplename).unwrap();

        write!(writer_r1, "{}", r1.to_string(), ).unwrap();
        write!(writer_r2, "{}", r2.to_string(), ).unwrap();

        if counter % 1_000_000 == 0{
            pbar.inc(1_000_000);
        }
    }
}

pub  fn demux_dual_index(sample_indices_fnames: HashMap<DualIndex, (String, String)>, undetermined_fname: (String,String), i1_list: Vec<String>, i2_list: Vec<String>, r1_list: Vec<String>, r2_list: Vec<String>) {
    let i1_iter = crate::io::fastq_list_iter(&i1_list);
    let i2_iter = crate::io::fastq_list_iter(&i2_list);
    let r1_iter = crate::io::fastq_list_iter(&r1_list);
    let r2_iter = crate::io::fastq_list_iter(&r2_list);

    let empty_index = DualIndex("".to_string(), "".to_string());

    let mut writers: HashMap<DualIndex, (Box<dyn Write>, Box<dyn Write>)> = HashMap::new();
    for (ix,(fname_r1, fname_r2)) in sample_indices_fnames.iter() {
        let encoder_r1 = BufWriter::new(
            bgzf::Writer::from_path_with_level(fname_r1, CompressionLevel::Fastest).unwrap()
        );
        let encoder_r2 = BufWriter::new(
            bgzf::Writer::from_path_with_level(fname_r2, CompressionLevel::Fastest).unwrap()
        );        
        writers.insert(ix.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));
    }
    // add the writer for unassigned
    let encoder_r1 = BufWriter::new(
        bgzf::Writer::from_path_with_level(undetermined_fname.0, CompressionLevel::Fastest).unwrap()
    );
    let encoder_r2 = BufWriter::new(
        bgzf::Writer::from_path_with_level(undetermined_fname.1, CompressionLevel::Fastest).unwrap()
    );    
    writers.insert(empty_index.clone(), (Box::new(encoder_r1), Box::new(encoder_r2)));


    let pbar = get_spinner();
    for (counter, (i1, i2, r1, r2)) in izip!(i1_iter, i2_iter, r1_iter, r2_iter).enumerate() {
        let key = DualIndex(i1.seq, i2.seq);

        let (writer_r1, writer_r2) = match writers.get_mut(&key){ // if Index not present, default to the undetermined
            Some(writer) => writer,
            None => writers.get_mut(&empty_index).unwrap()
        };
        
        write!(writer_r1, "{}", r1.to_string(), ).unwrap();
        write!(writer_r2, "{}", r2.to_string(), ).unwrap();

        if counter % 1_000_000 == 0{
            pbar.inc(1_000_000);
        }
    }
}

/// From a csv with i5,i7, R1filename, R2filename
/// 
pub fn samplesheet_to_hashmap(fname: &str) -> HashMap<DualIndex, (String, String)> {
    
    let mut sheet = HashMap::new();
    let mut rdr = csv::Reader::from_reader(BufReader::new(File::open(fname).unwrap()));
    for result in rdr.records() {
        let record = result.expect("a CSV record");

        let i1 = record.get(0).unwrap();
        let i2 = record.get(1).unwrap();
        let f1 = record.get(2).unwrap();
        let f2 = record.get(3).unwrap();

        sheet.insert(
            DualIndex(i1.to_string(), i2.to_string()), 
            (f1.to_string(), f2.to_string())
        );
    };
    sheet
}

/// From a csv with i5,i7, Prefix
pub fn samplesheet_to_hashmap_2(fname: &Path) -> HashMap<DualIndex, Samplename> {
    
    let mut sheet = HashMap::new();
    let mut rdr = csv::Reader::from_reader(BufReader::new(File::open(fname).unwrap()));
    for result in rdr.records() {
        let record = result.expect("a CSV record");

        let i1 = record.get(0).unwrap();
        let i2 = record.get(1).unwrap();
        let prefix = record.get(2).unwrap();

        sheet.insert(
            DualIndex(i1.to_string(), i2.to_string()), 
            Samplename(prefix.to_string())
        );
    };
    sheet
}