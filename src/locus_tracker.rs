use crate::cmd_cluster::Range;
use anyhow::Error;
use rust_htslib::bam::HeaderView;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::io::Write;

#[derive(Clone)]
pub struct Alignment {
    pub tid: i32,
    pub depth: usize,
    pub span: Range,
}

pub struct LocusTracker {
    pub map: HashMap<String, Vec<Alignment>>,
}

pub struct AlignmentReport {
    pub alns: HashMap<String, (Vec<Alignment>, HashSet<i32>)>,
}

#[derive(Serialize)]
struct AlignmentReportRow<'a> {
    source: &'a str,
    target: &'a str,
    depth: usize,
    loci: usize,
    references: String,
}

impl AlignmentReport {
    pub fn serialize<W: Write>(
        &self,
        writer: &mut csv::Writer<W>,
        header: &HeaderView,
        source: &str,
    ) -> Result<(), Error> {
        let mut targets: Vec<_> = self.alns.iter().collect();
        targets.sort_by_key(|(target, _)| *target);

        for (target, (alignments, tids)) in targets {
            let mut references: Vec<_> = tids
                .iter()
                .map(|tid| String::from_utf8_lossy(header.tid2name(*tid as u32)).into_owned())
                .collect();
            references.sort();

            writer.serialize(AlignmentReportRow {
                source,
                target,
                depth: alignments.iter().map(|alignment| alignment.depth).sum(),
                loci: alignments.len(),
                references: references.join(";"),
            })?;
        }

        writer.flush()?;
        Ok(())
    }
}

impl LocusTracker {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }

    pub fn add(&mut self, target: impl ToString, tid: i32, range: Range) {
        self.map
            .entry(target.to_string())
            .or_default()
            .push(Alignment {
                tid,
                span: range,
                depth: 1,
            });
    }

    pub fn resolve(&mut self, min_calls_per_loci: usize, read_len: usize) -> AlignmentReport {
        assert!(read_len > 0);

        let mut report = AlignmentReport {
            alns: HashMap::with_capacity(self.map.len()),
        };

        for (k, v) in self.map.iter_mut() {
            v.sort_by(|a, b| {
                a.span
                    .start
                    .cmp(&b.span.start)
                    .then(a.span.end.cmp(&b.span.end))
            });

            let mut dest: Vec<Alignment> = Vec::with_capacity(v.len());
            let mut prev = v[0].clone();
            let mut set = HashSet::from([prev.tid]);

            for cur in v[1..].iter() {
                if cur.span.start <= prev.span.end {
                    prev.span.end = std::cmp::max(cur.span.end, prev.span.end);
                    prev.depth += cur.depth
                } else {
                    dest.push(prev);
                    prev = cur.clone();
                }

                set.insert(cur.tid);
            }

            dest.push(prev);

            let total: usize = dest.iter().map(|d| d.span.len() / read_len).sum();

            if total >= min_calls_per_loci {
                report.alns.insert(k.clone(), (dest, set));
            }

            // *v = dest;
        }

        report
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rust_htslib::bam::{Header, header::HeaderRecord};

    fn assert_spans(report: &AlignmentReport, target: &str, expected: &[Range]) {
        let actual: Vec<_> = report.alns[target]
            .0
            .iter()
            .map(|alignment| alignment.span)
            .collect();

        assert_eq!(actual, expected);
    }

    #[test]
    pub fn test_resolve_1() {
        let mut lt = LocusTracker::new();

        let r1 = Range { start: 0, end: 350 };
        let r2 = Range { start: 5, end: 120 };
        let r3 = Range {
            start: 470,
            end: 530,
        };
        let r4 = Range { start: 0, end: 1 };

        lt.add(1, 0, r3);
        lt.add(1, 0, r2);
        lt.add(1, 0, r1);
        lt.add(2, 1, r4);

        let report = lt.resolve(0, 1);
        assert_spans(&report, "1", &[Range { start: 0, end: 350 }, r3]);
        assert_spans(&report, "2", &[Range { start: 0, end: 1 }]);
    }

    #[test]
    fn test_resolve_sorts_disjoint_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, 0, Range { start: 20, end: 30 });
        lt.add(1, 0, Range { start: 0, end: 10 });
        lt.add(1, 0, Range { start: 40, end: 50 });

        let report = lt.resolve(0, 1);

        assert_spans(
            &report,
            "1",
            &[
                Range { start: 0, end: 10 },
                Range { start: 20, end: 30 },
                Range { start: 40, end: 50 },
            ],
        );
    }

    #[test]
    fn test_resolve_merges_touching_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, 0, Range { start: 10, end: 20 });
        lt.add(1, 0, Range { start: 0, end: 10 });

        let report = lt.resolve(0, 1);

        assert_spans(&report, "1", &[Range { start: 0, end: 20 }]);
    }

    #[test]
    fn test_resolve_merges_nested_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, 0, Range { start: 5, end: 10 });
        lt.add(1, 0, Range { start: 0, end: 20 });
        lt.add(1, 0, Range { start: 5, end: 15 });

        let report = lt.resolve(0, 1);

        assert_spans(&report, "1", &[Range { start: 0, end: 20 }]);
    }

    #[test]
    fn test_resolve_merges_overlap_chain() {
        let mut lt = LocusTracker::new();
        lt.add(1, 0, Range { start: 8, end: 15 });
        lt.add(1, 0, Range { start: 0, end: 5 });
        lt.add(1, 0, Range { start: 4, end: 10 });

        let report = lt.resolve(0, 1);

        assert_spans(&report, "1", &[Range { start: 0, end: 15 }]);
    }

    #[test]
    fn test_resolve_handles_targets_independently() {
        let mut lt = LocusTracker::new();
        lt.add(1, 0, Range { start: 0, end: 10 });
        lt.add(
            2,
            1,
            Range {
                start: 100,
                end: 120,
            },
        );
        lt.add(1, 0, Range { start: 5, end: 15 });
        lt.add(
            2,
            1,
            Range {
                start: 110,
                end: 130,
            },
        );

        let report = lt.resolve(0, 1);

        assert_spans(&report, "1", &[Range { start: 0, end: 15 }]);
        assert_spans(
            &report,
            "2",
            &[Range {
                start: 100,
                end: 130,
            }],
        );
    }

    #[test]
    fn serializes_report_as_csv() {
        let mut lt = LocusTracker::new();
        lt.add("species", 0, Range { start: 0, end: 10 });
        lt.add("species", 1, Range { start: 20, end: 30 });
        let report = lt.resolve(0, 1);
        let header = test_header();
        let mut writer = csv::Writer::from_writer(Vec::new());

        report
            .serialize(&mut writer, &header, "sample.bam")
            .unwrap();
        report
            .serialize(&mut writer, &header, "sample-2.bam")
            .unwrap();

        let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output,
            "source,target,depth,loci,references\nsample.bam,species,2,2,ref-a;ref-b\nsample-2.bam,species,2,2,ref-a;ref-b\n"
        );
    }

    #[test]
    fn serializes_report_as_tsv() {
        let mut lt = LocusTracker::new();
        lt.add("species", 0, Range { start: 0, end: 10 });
        let report = lt.resolve(0, 1);
        let header = test_header();
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());

        report
            .serialize(&mut writer, &header, "sample.bam")
            .unwrap();

        let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output,
            "source\ttarget\tdepth\tloci\treferences\nsample.bam\tspecies\t1\t1\tref-a\n"
        );
    }

    fn test_header() -> HeaderView {
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "ref-a")
                .push_tag(b"LN", 100),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "ref-b")
                .push_tag(b"LN", 100),
        );
        HeaderView::from_header(&header)
    }
}
