
use crate::cmd_build_db::base_is_nonambig;
use rust_htslib::bam::record::{Record, Cigar, CigarString};
use anyhow::Error;

#[inline(always)]
pub fn filter_read(
    rec: &Record,
    min_frac_bases_aligned: f32,
    min_frac_bases_matched: f32,
) -> Result<bool, Error> {
    let mut bases_aligned: usize = 0;
    let mut bases_matched: usize = 0;

    let len = rec.seq_len();
    let mina = (min_frac_bases_aligned * len as f32) as usize;
    let minm = (min_frac_bases_matched * len as f32) as usize;
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

    Ok(bases_aligned >= mina && bases_matched >= minm)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_filter_read1() {
        let mut rec = Record::new();
        let seq = "GGTCACTGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACGGAGTCTCACTCTGTCGCCCAGGCTGGAGTGCA";
        let cigar = "1X3=1I5=52X2I2=1X33=";
        let qual = vec![40; seq.len()];
        rec.set(b"test1", Some(&CigarString::try_from(cigar.as_bytes()).unwrap()), seq.as_bytes(), qual.as_slice());

        assert_eq!(filter_read(&rec, 0.7, 0.7).unwrap(), false);
    }
}
