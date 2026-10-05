use rust_htslib::bam::record::Record;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ReadHash {
    name_hash: u64,
    source_hash: u64,
}

impl ReadHash {
    pub fn from_read_and_source(rec: &Record, source: &str) -> Self {
        let mut hash = std::hash::DefaultHasher::new();
        rec.qname().hash(&mut hash);
        let name_hash = hash.finish();

        let mut hash = std::hash::DefaultHasher::new();
        source.hash(&mut hash);
        let source_hash = hash.finish();

        Self { name_hash, source_hash }
    }
}

pub struct ReadBucket {
    pub hash: ReadHash,
    pub rec: Record,
    pub mate: Option<Record>,
}

#[derive(Default)]
pub struct SamMergeBuffer {
    q: Vec<ReadBucket>,
    read_ids: HashSet<u64>,
}

impl SamMergeBuffer {
    pub fn intake(&mut self, record: &Record, src: &str) -> Option<Vec<ReadBucket>> {
        let hash = ReadHash::from_read_and_source(record, src);

        if self.read_ids.contains(&hash.name_hash) {
            // already encountered this mate. So add and keep going.
            self.update_mate(record, hash);
            None
        } else {
            // we encountered a read QNAME we haven't seen before. Assuming a name-sorted BAM file,
            // this means we can emit the reads we've been holding in the queue.
            let ret = std::mem::take(&mut self.q);
            self.read_ids.clear();
            self.read_ids.insert(hash.name_hash);
            self.q.push(ReadBucket { hash, rec: record.clone(), mate: None });
            Some(ret)
        }
    }

    pub fn update_mate(&mut self, rec: &Record, hash: ReadHash) {
        for node in self.q.iter_mut() {
            if node.hash == hash {
                assert!(node.mate.is_none());
                node.mate = Some(rec.clone());
                return;
            }
        }

        // no mate found with same name and source, so update queue.
        self.q.push(ReadBucket { hash, rec: rec.clone(), mate: None });
    }

    pub fn flush(&mut self) -> Vec<ReadBucket> {
        std::mem::take(&mut self.q)
    }
}
