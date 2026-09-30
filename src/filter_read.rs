use crate::cmd_build_db::base_is_nonambig;
use anyhow::Error;
use rust_htslib::bam::record::{Cigar, Record};

pub struct ReadCover {
    bases_aligned: usize,
    bases_matched: usize,
    len: usize,
}

impl ReadCover {
    pub fn frac_aligned(&self) -> f32 {
        (self.bases_aligned / self.len) as f32
    }

    pub fn frac_matched(&self) -> f32 {
        (self.bases_matched / self.len) as f32
    }
}

#[inline(always)]
pub fn filter_read(rec: &Record) -> Result<ReadCover, Error> {
    let mut bases_aligned: usize = 0;
    let mut bases_matched: usize = 0;

    let len = rec.seq_len();
    let seq = rec.seq();
    let mut seqi: usize = 0;

    for op in &rec.cigar().0 {
        match op {
            Cigar::Match(_) => {
                anyhow::bail!("Wrong SAM format! Need X/= instead of M. Use SAM format 1.4+")
            }

            Cigar::Equal(len) => {
                for i in seqi..seqi + (*len as usize) {
                    if base_is_nonambig(seq[i]) {
                        bases_matched += 1;
                    }

                    bases_aligned += 1;
                }
                seqi += *len as usize;
            }

            Cigar::Diff(len) => {
                bases_aligned += *len as usize;
                seqi += *len as usize;
            }

            Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::Pad(len) => {
                seqi += *len as usize;
            }

            _ => (),
        }
    }
    #[cfg(test)]
    eprintln!("bases matched: {} bases aligned: {}, total: {}", bases_matched, bases_aligned, len);

    Ok(ReadCover { bases_aligned, bases_matched, len })
}

#[cfg(test)]
mod test {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    #[test]
    fn test_filter_read1() {
        let mut rec = Record::new();
        let seq =
            "GGTCACTGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACGGAGTCTCACTCTGTCGCCCAGGCTGGAGTGCA";
        let cigar = "1X3=1I5=52X2I2=1X33=";
        let qual = vec![40; seq.len()];
        rec.set(b"test1", Some(&CigarString::try_from(cigar.as_bytes()).unwrap()), seq.as_bytes(), qual.as_slice());

        let cov = filter_read(&rec).unwrap();
        assert!(cov.frac_aligned() < 0.7 || cov.frac_matched() < 0.7);
    }
}
