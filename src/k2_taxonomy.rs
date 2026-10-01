use crate::k2_utils::{K2ClassificationRow, K2ClassificationStatus, K2ReportRow};
use anyhow::{Error, bail};
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::BufRead;

#[derive(Debug, Clone)]
pub struct K2Taxon {
    pub parent: u32,
    pub row: K2ReportRow,
    pub children: Vec<u32>,
    pub index: usize,
}

impl std::fmt::Display for K2Taxon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let spacer = std::iter::repeat_n("-", self.row.depth * 2).collect::<String>();
        writeln!(f, "{spacer} taxon {}, {}", self.row.taxid, self.row.name)?;
        writeln!(f, "{spacer} children {:?}", self.children)?;
        writeln!(f, "{spacer} parent {}", self.parent)
    }
}

#[derive(Default, Debug)]
// NOTE: read names being stored twice, once in K2ClassificationRow and again in
// read_name_to_read_id
pub struct K2Taxonomy {
    pub tree: HashMap<u32, K2Taxon>,
    pub read_id_to_read: HashMap<u32, K2ClassificationRow>,
    pub taxid_to_direct_reads: HashMap<u32, Vec<u32>>,
    pub read_name_to_read_id: HashMap<Vec<u8>, u32>,
}

impl<'a> K2Taxonomy {
    pub fn ingest_report<C: BufRead>(&mut self, report: C) -> Result<(), Error> {
        let mut line_iter = report.lines();
        let mut queue: VecDeque<K2Taxon> = VecDeque::new();

        let mut parent = K2ReportRow::try_from(line_iter.next().unwrap().expect("read first line of report").as_str())?;

        // skip unclassified row.
        if parent.taxid == 0 {
            if let Some(line) = line_iter.next() {
                parent = K2ReportRow::try_from(line.expect("Reading second report line").as_str())?;
            } else {
                bail!("report empty");
            }
        }

        assert!(parent.taxid == 1);

        queue.push_back(K2Taxon { parent: 1, row: parent, children: vec![], index: 0 });

        for (i, l) in line_iter.enumerate() {
            let line = l.expect("read line from report file");
            let row = K2ReportRow::try_from(line.trim())?;

            // if we have lower depth than previous, we have read all children of a given level.
            while row.depth < queue.front().unwrap().row.depth {
                // we unwrap here because we expect the root (taxid = 1) to remain always
                let mut cur = queue.pop_front().unwrap();
                cur.parent = queue.front().unwrap().row.taxid;
                self.tree.insert(cur.row.taxid, cur);
            }

            let front = queue.front().unwrap();

            queue.push_front(K2Taxon { parent: front.row.taxid, row, children: vec![], index: i + 1 })
        }

        while let Some(mut cur) = queue.pop_front() {
            cur.parent = if let Some(parent) = queue.front() { parent.row.taxid } else { cur.row.taxid };
            self.tree.insert(cur.row.taxid, cur);
        }

        // second pass: create child arrays
        for (taxid, parent) in self.tree.values().map(|f| (f.row.taxid, f.parent)).collect::<Vec<(u32, u32)>>() {
            if taxid == 1 {
                continue;
            };
            let parent = self.tree.get_mut(&parent).expect("child to parent taxid lookup");
            parent.children.push(taxid);
        }

        Ok(())
    }

    pub fn rank_nodes(&'a self, rank: &str) -> impl Iterator<Item = &'a K2Taxon> {
        let mut nodes = self.tree.values();
        std::iter::from_fn(move || {
            while let Some(node) = nodes.next() {
                if node.row.rank == rank {
                    return Some(node);
                }
            }

            None
        })
    }

    pub fn get_kmer_calls(&self, read_name: &[u8]) -> Option<&K2ClassificationRow> {
        self.read_name_to_read_id.get(read_name).and_then(|id| self.read_id_to_read.get(id))
    }

    pub fn get_all_reads_under_node(&'a self, taxid: u32) -> impl Iterator<Item = &'a K2ClassificationRow> {
        let mut ids: &[u32] = &[];
        let mut id_index = 0;

        let mut traverse = if let Some(node) = self.tree.get(&taxid) {
            std::iter::chain(
                std::iter::once(node),
                self.descendants(taxid).flatten().map(|f| self.tree.get(f).expect("get child")),
            )
            .collect::<Vec<&K2Taxon>>()
        } else {
            vec![]
        }
        .into_iter();

        std::iter::from_fn(move || {
            loop {
                if id_index < ids.len() {
                    let ret = self.get_read(ids[id_index]).expect("getting read by ID");
                    id_index += 1;
                    return Some(ret);
                } else if let Some(node) = traverse.next() {
                    id_index = 0;
                    ids = self.taxid_to_direct_reads.get(&node.row.taxid).map(|v| v.as_slice()).unwrap_or(&[]);
                } else {
                    return None;
                }
            }
        })
    }

    #[allow(dead_code)]
    pub fn display_tree<W: std::io::Write>(&self, mut writer: W) -> Result<(), Error> {
        let mut nodes = self.tree.values().cloned().collect::<Vec<K2Taxon>>();
        nodes.sort_by(|a, b| a.index.cmp(&b.index));
        for node in nodes {
            write!(writer, "{}", node)?;
        }

        Ok(())
    }

    pub fn ingest_read_list<C: BufRead>(&mut self, reads: C) -> Result<(), Error> {
        let mut i = 0;
        for l in reads.lines() {
            let line = l.expect("reading k2 classification txt");
            let row = K2ClassificationRow::try_from(line.trim())?;
            if row.status == K2ClassificationStatus::Unclassified {
                continue;
            };

            let idx = u32::try_from(i)?;
            i += 1;

            assert!(self.read_name_to_read_id.insert(row.read_name.as_bytes().to_vec(), idx).is_none());
            self.taxid_to_direct_reads.entry(row.assigned_taxid).or_default().push(idx);

            self.read_id_to_read.insert(idx, row);
        }
        Ok(())
    }

    pub fn sanity_check_counts(&self) -> Result<(), Error> {
        for v in self.tree.values() {
            let direct_reads =
                if let Some(taxon) = self.taxid_to_direct_reads.get(&v.row.taxid) { taxon.len() as u64 } else { 0 };

            if direct_reads != v.row.exact_reads {
                bail!(
                    "Validation error: Discordant direct read count between classification and report ({} vs {}): Taxon ID: {}, Taxon name: {}",
                    direct_reads,
                    v.row.exact_reads,
                    v.row.taxid,
                    v.row.name
                );
            }
        }

        for k in self.taxid_to_direct_reads.keys() {
            if !self.tree.contains_key(k) {
                // if self.tree.get(k).is_none() {
                bail!("Anomalous data: read classification has taxon ID not found in report: {k}")
            }
        }

        Ok(())
    }

    pub fn build<C: BufRead, R: BufRead>(reads: C, report: R) -> Result<Self, Error> {
        let mut ret = Self::default();
        ret.ingest_report(report)?;
        ret.ingest_read_list(reads)?;
        ret.sanity_check_counts()?;
        Ok(ret)
    }

    pub fn get_read(&self, id: u32) -> Option<&K2ClassificationRow> {
        self.read_id_to_read.get(&id)
    }

    pub fn get(&self, taxid: u32) -> Option<&K2Taxon> {
        self.tree.get(&taxid)
    }

    pub fn descendants(&self, taxid: u32) -> impl Iterator<Item = &[u32]> {
        let mut deque: VecDeque<u32> = VecDeque::from([taxid]);
        let mut seen: HashSet<u32> = HashSet::new();

        std::iter::from_fn(move || {
            while let Some(cur) = deque.pop_front() {
                if seen.insert(cur) {
                    let children = self.tree.get(&cur).map(|v| v.children.as_slice()).unwrap_or(&[]);

                    deque.extend(children);
                    return Some(children);
                }
            }

            None
        })
    }

    pub fn descendants_taxa(&self, taxid: u32) -> impl Iterator<Item = &K2Taxon> {
        self.descendants(taxid).flatten().map(|t| self.tree.get(t).expect("getting descendant taxon"))
    }

    pub fn lineage(&self, tax_id: u32) -> Option<Vec<&K2Taxon>> {
        let mut lineage = Vec::new();
        let mut seen = HashSet::new();
        let mut current = tax_id;

        loop {
            if !seen.insert(current) {
                return None;
            }

            let taxon = self.tree.get(&current).unwrap();
            lineage.push(taxon);

            if taxon.parent == current {
                break;
            }

            current = taxon.parent;
        }

        lineage.reverse();
        Some(lineage)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const REPORT: &str = concat!(
        "33.33\t1\t1\tU\t0\tunclassified\n",
        "66.67\t2\t0\tR\t1\troot\n",
        "66.67\t2\t1\tG\t10\t  Genus\n",
        "33.33\t1\t1\tS\t11\t    Species\n",
    );
    const CLASSIFICATIONS: &str = concat!(
        "U\tunclassified-read\t0\t100\t0:10\n",
        "C\tgenus-read\t10\t100|100\t10:2 11:1 |:| A:1\n",
        "C\tspecies-read\t11\t100\t11:3 11:2 999:1\n",
    );

    #[test]
    fn builds_tree_with_unified_read_indices() {
        let taxonomy =
            K2Taxonomy::build(Cursor::new(CLASSIFICATIONS.as_bytes()), Cursor::new(REPORT.as_bytes())).unwrap();

        assert_eq!(taxonomy.read_id_to_read.len(), 2);
        let read0 = taxonomy.get_read(0).unwrap();
        let tax10 = taxonomy.get(10).unwrap();
        let tax11 = taxonomy.get(11).unwrap();

        assert_eq!(read0.read_name, "genus-read");

        assert_eq!(taxonomy.taxid_to_direct_reads.get(&tax10.row.taxid), Some(vec![0u32].as_ref()));

        assert_eq!(taxonomy.taxid_to_direct_reads.get(&tax11.row.taxid), Some(vec![1u32].as_ref()));

        assert_eq!(taxonomy.tree.get(&tax11.row.taxid).map(|node| node.children.as_slice()), Some(vec![].as_slice()));

        assert_eq!(
            taxonomy.tree.get(&tax10.row.taxid).map(|node| node.children.as_slice()),
            Some(vec![11_u32].as_slice())
        );

        eprintln!("{:?}", taxonomy);
        let lineage = taxonomy.lineage(11).unwrap();
        assert_eq!(lineage.iter().map(|node| node.row.taxid).collect::<Vec<_>>(), vec![1, 10, 11]);

        assert_eq!(taxonomy.descendants(1).flatten().copied().collect::<Vec<_>>(), vec![10, 11].as_slice());
    }

    #[test]
    fn rejects_report_count_mismatches() {
        let classifications =
            CLASSIFICATIONS.to_string() + "C\tgenus-read2\t10\t100\t10:2\nC\textra-read\t10\t100\t10:2\n";

        let error = K2Taxonomy::build(Cursor::new(classifications), Cursor::new(REPORT.as_bytes())).unwrap_err();

        eprintln!("{error}");
        assert!(error.to_string().contains("Genus"))
    }

    #[test]
    fn rejects_assignments_missing_from_report() {
        let classifications = CLASSIFICATIONS.to_string() + "C\tread\t99\t100\t99:1\n";
        let error = K2Taxonomy::build(Cursor::new(classifications), Cursor::new(REPORT.as_bytes())).unwrap_err();

        eprintln!("{error}");
        assert!(error.to_string().contains("report: 99"));
    }
}
