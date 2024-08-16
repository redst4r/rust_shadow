use core::panic;
use rust_htslib::bgzf;
use rust_htslib::bgzf::CompressionLevel;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;

fn switch_base(base: char) -> char{
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        _ => panic!("unknown base")
    }
}

pub fn reverse_complement(seq: &str) -> String {
    let mut rc = String::with_capacity(seq.len());

    for c in seq.chars().rev() {
        rc.push(switch_base(c))
    }
    rc
}

/// A single FastQ entry, with header, sequence and quality scores
pub struct FastqEntry {
    pub header: String,
    pub seq: String,
    pub phred: String,
}

impl FastqEntry {
    /// Turns the FastQ entry intro a String representation that can directly be written
    /// to a fastq file 
    pub fn to_string(&self) -> String {
        // format is much slower!!
        // format!("{}\n{}\n+\n{}\n", self.header, self.seq, self.phred)
        let mut s = String::with_capacity(self.header.len() + self.seq.len() * 2 + 5); //4newlines and a +
        s.push_str(&self.header);
        s.push('\n');
        s.push_str(&self.seq);
        s.push_str("\n+\n");
        s.push_str(&self.phred);
        s.push('\n');
        s
    }
}
/// Iterator over a fastq.gz file, yielding [`FastEntry`]
pub struct FastIterator {
    reader: BufReader<bgzf::Reader>,
}

impl FastIterator {
    pub fn new(fastqname: &str) -> Self {
        let decoder = bgzf::Reader::from_path(fastqname).unwrap();
        let reader = BufReader::with_capacity(800 * 1024, decoder);
        FastIterator { reader }
    }
}

impl Iterator for FastIterator {
    type Item = FastqEntry;

    fn next(&mut self) -> Option<Self::Item> {
        let mut header = String::new();
        // try to read a header
        match self.reader.read_line(&mut header) {
            Ok(0) => None,
            Ok(_n) => {
                let mut seq = String::new();
                self.reader.read_line(&mut seq).unwrap();

                let mut dummy = String::new();
                self.reader.read_line(&mut dummy).unwrap();

                let mut phred = String::new();
                self.reader.read_line(&mut phred).unwrap();

                let fq = FastqEntry {
                    header: header.trim().to_string(),
                    seq: seq.trim().to_string(),
                    phred: phred.trim().to_string(),
                };
                Some(fq)
            }
            Err(e) => panic!("{}", e),
        }
    }
}


use once_cell::sync::Lazy;
pub static PHRED_LOOKUP: Lazy<PhredCache> = Lazy::new(|| {
    let lookup = PhredCache::new();
    lookup
});



/// Caches the Phred symbol to probability translation table
/// TODO: could be done using lazy_static
pub struct PhredCache {
    cache: Vec<f32>,
}
impl PhredCache {
    pub fn new() -> Self {
        let mut cache: Vec<f32> = Vec::new();
        for i in 33..76 {
            let c: char = i.into();
            let p = phred_symbol_to_prob(c);
            cache.push(p);
        }
        PhredCache { cache }
    }
    pub fn get_prob(&self, c: char) -> f32 {
        let i = (c as u32) - 33;
        self.cache[i as usize]
    }
}

// use cached::proc_macro::cached;
// #[cached]
fn phred_symbol_to_prob(phred: char) -> f32 {
    let q = (phred as u32) - 33;
    10_f32.powf(-(q as f32) / 10_f32)
}

// fn avg_phred(phred: &str) -> f32{
//     let n_chars = phred.len() as f32;
//     let summed_probs: f32 = phred.chars().map(|c| phred_symbol_to_prob(c)).sum();
//     summed_probs / n_chars
//     // 0.00000000001
// }

// def Phred2symbol(phred:str):
//     "phred score to ascii"
//     return str(chr(phred+33))

/// Filters a fastq-file for all reads having an aggregated PhredScore of > `threshold_qc`
/// # Parameters:
/// * fastqname: File to be filtered
/// * outname: File where to write the filtered records
/// * threshold_qc: minimum  (aggreated) Phred Score a read needs to pass to get written
pub fn quality_filter(fastqname: &str, outname: &str, threshold_qc: f32) {
    // reading the fastq
    let fastq_iter = FastIterator::new(fastqname);

    let cache = PhredCache::new();

    let encoder = bgzf::Writer::from_path_with_level(outname, CompressionLevel::Fastest).unwrap();
    let mut writer = BufWriter::new(encoder);

    let mut total_reads = 0;
    let mut passing_reads = 0;

    for fq in fastq_iter {
        total_reads += 1;

        let probs: f32 = fq.phred.chars().map(|c| cache.get_prob(c)).sum();
        let avg_qual = probs / (fq.phred.len() as f32);

        if avg_qual < threshold_qc {
            write!(writer, "{}", fq.to_string()).unwrap();
            passing_reads += 1;
        }
    }
    println!(
        "{}/{}({}) reads passed QC",
        passing_reads,
        total_reads,
        (passing_reads as f32) / (total_reads as f32)
    )
}

// zcat kraken_out.filtered.gz | awk '{ print $2}' | less
pub fn read_filter_whitelist(fastqname: &str, outname: &str, whitelist: &str) {
    let whitelist_reader = BufReader::new(File::open(whitelist).unwrap());
    let whitelist_header: HashSet<String> = whitelist_reader.lines().map(|f| f.unwrap()).collect();

    // reading the fastq
    let fastq_iter = FastIterator::new(fastqname);

    //writing the filtered
    let encoder = bgzf::Writer::from_path_with_level(outname, CompressionLevel::Fastest).unwrap();
    let mut writer = BufWriter::new(encoder);

    let mut total_reads = 0;
    let mut passing_reads = 0;

    for fq in fastq_iter {
        total_reads += 1;

        if whitelist_header.contains(&fq.header) {
            write!(writer, "{}", fq.to_string()).unwrap();
            passing_reads += 1;
        }
    }
    println!(
        "{}/{}({}) reads were whitelisted",
        passing_reads,
        total_reads,
        (passing_reads as f32) / (total_reads as f32)
    )
}

/// Chaining many fastq files into a single iterator
pub fn fastq_list_iter(fastq_list: &[String]) -> impl Iterator<Item = FastqEntry> + '_ {
    let my_iter = fastq_list
        .iter()
        .flat_map(move |fname| FastIterator::new(fname));
    my_iter
}

/// Chaining many fastq files into a single iterator
pub fn fastq_phred_iter(fastq_list: &[String]) -> impl Iterator<Item = String> + '_ {
    // instead if yielding the sequence, this one yields the PHRED ASCII scores of the reads
    fastq_list_iter(fastq_list).map(|fq| fq.phred)
}

/// Loading 10x CB whilelist from file
/// Returns a HashSet of CBs
pub fn parse_whitelist_gz(fname: &String) -> HashSet<String> {
    let decoder = bgzf::Reader::from_path(fname).unwrap();
    let reader = BufReader::new(decoder);
    let my_iter = reader.lines(); //.take(10_000_000);
    let mut hset: HashSet<String> = HashSet::new();
    for line in my_iter.flatten() {
        // flatten filters out the Error elements
        hset.insert(line);
    }
    hset
}

#[cfg(test)]
mod testing {
    use crate::io::reverse_complement;

    // #[test]
    use super::{fastq_list_iter, quality_filter, FastqEntry, PhredCache};
    use rust_htslib::bgzf;
    use rust_htslib::bgzf::CompressionLevel;
    use std::io::BufWriter;
    use std::io::Write;

    fn test_make_fastq() {
        let n = 1_000_000_usize;
        let out = "/tmp/test.fastq.gz";
        let encoder = bgzf::Writer::from_path_with_level(out, CompressionLevel::Fastest).unwrap();
        let mut writer = BufWriter::new(encoder);

        // let seq_len = 150;
        let dummyseq = "A".repeat(150);
        let dummphred = "F".repeat(150);

        for i in 0..n {
            let fq = FastqEntry {
                header: format!("@Read{i}"),
                seq: dummyseq.clone(),
                phred: dummphred.clone(),
            };
            write!(writer, "{}", fq.to_string()).unwrap();
        }
    }
    #[test]
    fn test_PhredCache() {
        let cache = PhredCache::new();
        assert_eq!(0.0001, cache.get_prob('I')); //Q40
        assert_eq!(0.001, cache.get_prob('?')); //Q30
        assert_eq!(0.01, cache.get_prob('5')); //Q20
        assert_eq!(0.1, cache.get_prob('+')); //Q10
        assert_eq!(1_f32, cache.get_prob('!'));
    }
    // #[test]
    pub fn test_filter() {
        // let file = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/DSP1/DSP1_CKDL210025651-1a-SI_TT_A2_HVWMHDSX2_S4_L001_R1_001.fastq.gz";
        let file = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/raw_data/Ice1/Ice1_CKDL210025651-1a-SI_TT_D2_HVWMHDSX2_S8_L001_R2_001.fastq.gz";
        // let file = "/tmp/test.fastq.gz";
        let out = "/tmp/filtered.fastq.gz";

        println!("Filtering!");
        use std::time::Instant;
        let now = Instant::now();
        quality_filter(file, out, 0.01);
        let elapsed_time = now.elapsed();
        println!("Running took {} sec.", elapsed_time.as_secs());
    }
    #[test]
    fn test_fastq() {
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
        // TODO there's an issue with trailing lines!!

        use std::fs::File;
        use std::io::BufWriter;
        use std::io::Write;

        let fastqname = "/tmp/foo.fastq";
        let f = File::create(fastqname).expect("Unable to create file");
        let mut f = BufWriter::new(f);

        f.write_all(fastq_entry1.as_bytes())
            .expect("Unable to write data");
        f.write_all(fastq_entry2.as_bytes())
            .expect("Unable to write data");
        f.flush().unwrap();

        let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()])
            .map(|fq| fq.header)
            .collect();
        assert_eq!(lines, vec!["@some read id", "@another read id"]);

        let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()])
            .map(|fq| fq.seq)
            .collect();
        assert_eq!(
            lines,
            vec![
                "AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC",
                "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"
            ]
        );

        let lines: Vec<_> = fastq_list_iter(&vec![fastqname.to_string()])
            .map(|fq| fq.phred)
            .collect();
        assert_eq!(
            lines,
            vec![
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
            ]
        );
    }
    #[test]
    fn test_rc(){
        assert_eq!(
            reverse_complement("AAGG"), "CCTT"
        );
        assert_eq!(
            reverse_complement("AAAA"), "TTTT"
        );
        assert_eq!(
            reverse_complement("ATGC"), "GCAT"
        );
    }
    
}
