use crate::Range;
use std::collections::HashMap;

#[derive(Clone)]
pub struct Alignment {
    pub depth: usize,
    pub span: Range,
}

pub struct LocusTracker {
    pub map: HashMap<String, Vec<Alignment>>,
}

impl LocusTracker {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }

    pub fn add(&mut self, target: impl ToString, range: Range) {
        self.map
            .entry(target.to_string())
            .or_default()
            .push(Alignment {
                span: range,
                depth: 1,
            });
    }

    pub fn resolve(&mut self) {
        for (_, v) in self.map.iter_mut() {
            v.sort_by(|a, b| {
                a.span
                    .start
                    .cmp(&b.span.start)
                    .then(a.span.end.cmp(&b.span.end))
            });

            let mut dest: Vec<Alignment> = Vec::with_capacity(v.len());
            let mut prev = v[0].clone();

            for cur in v[1..].iter() {
                if cur.span.start <= prev.span.end {
                    prev.span.end = std::cmp::max(cur.span.end, prev.span.end);
                    prev.depth += cur.depth
                } else {
                    dest.push(prev);
                    prev = cur.clone();
                }
            }

            dest.push(prev);
            *v = dest;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn assert_alignments(lt: &LocusTracker, target: &str, expected: &[(Range, usize)]) {
        let actual: Vec<_> = lt.map[target]
            .iter()
            .map(|alignment| (alignment.span, alignment.depth))
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

        lt.add(1, r3);
        lt.add(1, r2);
        lt.add(1, r1);
        lt.add(2, r4);

        assert_alignments(&lt, "1", &[(r3, 1), (r2, 1), (r1, 1)]);
        assert_alignments(&lt, "2", &[(r4, 1)]);

        lt.resolve();
        assert_alignments(
            &lt,
            "1",
            &[(Range { start: 0, end: 350 }, 2), (r3, 1)],
        );
        assert_alignments(&lt, "2", &[(Range { start: 0, end: 1 }, 1)]);
    }

    #[test]
    fn test_resolve_sorts_disjoint_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 20, end: 30 });
        lt.add(1, Range { start: 0, end: 10 });
        lt.add(1, Range { start: 40, end: 50 });

        lt.resolve();

        assert_alignments(
            &lt,
            "1",
            &[
                (Range { start: 0, end: 10 }, 1),
                (Range { start: 20, end: 30 }, 1),
                (Range { start: 40, end: 50 }, 1),
            ],
        );
    }

    #[test]
    fn test_resolve_merges_touching_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 10, end: 20 });
        lt.add(1, Range { start: 0, end: 10 });

        lt.resolve();

        assert_alignments(&lt, "1", &[(Range { start: 0, end: 20 }, 2)]);
    }

    #[test]
    fn test_resolve_merges_nested_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 5, end: 10 });
        lt.add(1, Range { start: 0, end: 20 });
        lt.add(1, Range { start: 5, end: 15 });

        lt.resolve();

        assert_alignments(&lt, "1", &[(Range { start: 0, end: 20 }, 3)]);
    }

    #[test]
    fn test_resolve_merges_overlap_chain() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 8, end: 15 });
        lt.add(1, Range { start: 0, end: 5 });
        lt.add(1, Range { start: 4, end: 10 });

        lt.resolve();

        assert_alignments(&lt, "1", &[(Range { start: 0, end: 15 }, 3)]);
    }

    #[test]
    fn test_resolve_handles_targets_independently() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 0, end: 10 });
        lt.add(
            2,
            Range {
                start: 100,
                end: 120,
            },
        );
        lt.add(1, Range { start: 5, end: 15 });
        lt.add(
            2,
            Range {
                start: 110,
                end: 130,
            },
        );

        lt.resolve();

        assert_alignments(&lt, "1", &[(Range { start: 0, end: 15 }, 2)]);
        assert_alignments(&lt, "2", &[(Range { start: 100, end: 130 }, 2)]);
    }
}
