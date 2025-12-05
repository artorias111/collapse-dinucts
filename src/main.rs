use std::fs::File;
use std::path::Path;
use std::io::{self, prelude::*, BufReader};
use flate2::read::MultiGzDecoder;
use clap::Parser;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)]
    reads: String,

    #[arg(short, long, default_value_t = 0.0)]
    threshold: f32,
}

fn main() {
    let args = Args::parse();
    let path = &args.reads;
    let threshold = &args.threshold;
    let _ = stream_fastq(path, *threshold);
}


fn stream_fastq(filepath: &str, threshold: f32) -> io::Result<()> {
    let is_gzipped:i32;
    let filepath = Path::new(filepath);

    if filepath.extension().unwrap() == "gz" {
        is_gzipped = 1;
    } else {
        is_gzipped = 0;
    }

    let _display = filepath.display();

    let file = match File::open(&filepath) {
        Err(why) => panic!("couldn't open the reads file! {why}"),
        Ok(file) => file,
    };

    let reader: Box<dyn BufRead> = if is_gzipped == 1 {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };


    let mut is_header_assigned = 0;
    let mut header = String::new();

    for line in reader.lines() {
        let fastq_line = match line {
            Err(_) => panic!("Empty fastq file, please try again"),
            Ok(file) => file,
        };
        if fastq_line.starts_with('@') {
            header = fastq_line;
            is_header_assigned = 1;
            continue;
        }

        if is_header_assigned == 1 && threshold == 0.0 {
            let collapsed_line = collapse_dinuct(&fastq_line);
            println!(">{header}_dc"); // dinucleotide collapsed
            println!("{collapsed_line}");
            is_header_assigned = 0;
        } else if is_header_assigned == 1 && threshold > 0.0 { 
            let retained_string = match collapse_dinuct_threshold(&fastq_line, threshold) {
                Some(line) => line,
                None => String::from("not retained")
            };

            if retained_string != "not retained" {
                println!(">{header}"); 
                println!("{retained_string}");
            }
            is_header_assigned = 0;
        }

    }

    Ok(())
}

fn collapse_dinuct(fastq_entry: &str) -> String {

    let mut chars_iter = fastq_entry.chars();
    let mut should_collapse = 0;
    let mut new_str = String::new();
    let mut prev_c1 = 'A';
    let mut prev_c2 = 'A';

    while let (Some(c1), Some(c2)) = (chars_iter.next(), chars_iter.next()) {
        if c1 == prev_c1 {
            if c2 == prev_c2 {
                if c1 != c2 {
                    should_collapse = 1;
                }
            }
        }

        prev_c1 = c1;
        prev_c2 = c2;

        if should_collapse == 0 {
            new_str.push(c1);
            new_str.push(c2);
        }
        should_collapse = 0;
    }

    new_str
}


fn collapse_specific_dincuelotide(query: &str, fastq_entry:&str) -> String {
    let rc_query = revcomp(query);
    let flipped_query = flip_dinuc(query);
    let rc_flipped_query = revcomp(&flipped_query);
    todo!("implement collapse specific dinuc");
}



fn collapse_dinuct_threshold(fastq_entry: &str, threshold: f32) -> Option<String> { // if the dinucleotide
                                                                                // percentage is above
                                                                                // a threshold, discard
                                                                                // the whole read
    let mut cur_count: f32 = 0.0; // current count of dincleotides
    let mut chars_iter = fastq_entry.chars();
    let mut prev_c1 = 'A';
    let mut prev_c2 = 'A';
    let mut len: f32 = 0.0;

    while let (Some(c1), Some(c2)) = (chars_iter.next(), chars_iter.next()) {
        len += 2.0;
        if c1 == prev_c1 {
            if c2 == prev_c2 {
                if c1 != c2 {
                    cur_count += 2.0;
                    
                }
            }
        }
        prev_c1 = c1;
        prev_c2 = c2
    }

    let percentage_count = cur_count*100.0/len;
    if percentage_count > threshold {
        return None
    }
    // println!("{cur_count}, {len}, {percentage_count}"); // debug print
    Some(String::from(fastq_entry))
}


fn revcomp(nucleotide: &str) -> String {
    let mut chars_iter = nucleotide.chars().rev();
    let mut revcomp_str = String::new();
    while let Some(c) = chars_iter.next() {
        let x = match c {
            'A' => 'T',
            'G' => 'C',
            'T' => 'A',
            'C' => 'G',
            'N' => 'N',
            _ => 'N'
        };
        revcomp_str.push(x);
    }
    revcomp_str
}

fn flip_dinuc(dinuc: &str) -> String { // 'AG' -> 'GA' and so on
    dinuc.chars().rev().collect()
}


// tests

#[cfg(test)]
mod tests {
    use crate::collapse_dinuct;

    #[test]
    fn test1_collapse_dinuct_dev() {
        let sequence = "ATATAGGGCGAGACCCCCGAGAGA";
        collapse_dinuct(sequence);
    }

    #[test]
    fn test2_collapse_dinuct_dev() {
        let sequence = "GCGCGCGC";
        collapse_dinuct(sequence);
    }

    #[test]
    fn test3_collapse_dinuct_dev() {
        let sequence = "TTAGGCTTTGCGCAGTAGCGCGCGCGCGCGCGAATATATTATATATATATATATATATATATATATATATATATATATATATATATGGGGGGGGGGGCGCGCGCGCGCGCGATATATATATATATAAGAGAGAGAGAGAGAGTCTCTCTCTCTCTCTCTCTC";
        collapse_dinuct(sequence);
    }

    // revcomp tests
    use crate::revcomp;
    
    #[test]
    fn revcomp_test1() {
        assert_eq!(revcomp("ATATAGATA"), "TATCTATAT");
    }

    // fliptests

    use crate::flip_dinuc;

    #[test]
    fn flip_test1() {
        assert_eq!("GA", flip_dinuc("AG"));
    }

    // threshold tests

    use crate::collapse_dinuct_threshold;

    #[test]
    fn threshold_test1() {
        let sequence = "ATATAGGGCGAGACCCCCGAGAGA";
        let x = collapse_dinuct_threshold(sequence, 22.0);
        assert_eq!(x, None);
    }

}

