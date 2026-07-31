use crate::Range;
use std::collections::HashMap;

pub struct LocusTracker {
    pub map: HashMap<String, Vec<Range>>,
}

impl LocusTracker {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }

    pub fn add(&mut self, target: impl ToString, range: Range) {
        self.map.entry(target.to_string()).or_default().push(range);
    }

    pub fn resolve(&mut self) {
        for (_, v) in self.map.iter_mut() {
            v.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));

            let mut dest: Vec<Range> = Vec::with_capacity(v.len());
            let mut prev = v[0];

            for cur in v[1..].iter() {
                if cur.start <= prev.end {
                    prev.end = std::cmp::max(cur.end, prev.end);
                } else {
                    dest.push(prev);
                    prev = *cur;
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

        assert_eq!(lt.map.get("1"), Some(&vec![r3, r2, r1]));
        assert_eq!(lt.map.get("2"), Some(&vec![r4]));

        lt.resolve();
        assert_eq!(
            lt.map.get("1"),
            Some(&vec![Range { start: 0, end: 350 }, r3])
        );

        assert_eq!(lt.map.get("2"), Some(&vec![Range { start: 0, end: 1 }]));
    }

    #[test]
    fn test_resolve_sorts_disjoint_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 20, end: 30 });
        lt.add(1, Range { start: 0, end: 10 });
        lt.add(1, Range { start: 40, end: 50 });

        lt.resolve();

        assert_eq!(
            lt.map.get("1"),
            Some(&vec![
                Range { start: 0, end: 10 },
                Range { start: 20, end: 30 },
                Range { start: 40, end: 50 },
            ])
        );
    }

    #[test]
    fn test_resolve_merges_touching_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 10, end: 20 });
        lt.add(1, Range { start: 0, end: 10 });

        lt.resolve();

        assert_eq!(lt.map.get("1"), Some(&vec![Range { start: 0, end: 20 }]));
    }

    #[test]
    fn test_resolve_merges_nested_ranges() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 5, end: 10 });
        lt.add(1, Range { start: 0, end: 20 });
        lt.add(1, Range { start: 5, end: 15 });

        lt.resolve();

        assert_eq!(lt.map.get("1"), Some(&vec![Range { start: 0, end: 20 }]));
    }

    #[test]
    fn test_resolve_merges_overlap_chain() {
        let mut lt = LocusTracker::new();
        lt.add(1, Range { start: 8, end: 15 });
        lt.add(1, Range { start: 0, end: 5 });
        lt.add(1, Range { start: 4, end: 10 });

        lt.resolve();

        assert_eq!(lt.map.get("1"), Some(&vec![Range { start: 0, end: 15 }]));
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

        assert_eq!(lt.map.get("1"), Some(&vec![Range { start: 0, end: 15 }]));
        assert_eq!(
            lt.map.get("2"),
            Some(&vec![Range {
                start: 100,
                end: 130,
            }])
        );
    }
}
