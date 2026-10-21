use crate::{
    cmd_cluster::Range,
    filter_read::{ReadCover, filter_read},
};
use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, record::Record};
use std::hash::{Hash, Hasher};

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
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

/// min amount of info to keep in merge buffer. Don't want to clone records to store info I don't
/// need.
#[derive(Debug)]
pub struct ClusterAlignment {
    hash: ReadHash,
    pub seq_len: usize,
    pub name: Vec<u8>,
    pub forward: bool,
    pub tid: i32,
    pub taxid: u32,
    pub map_interval: Range,
    pub cover: ReadCover,
}

impl ClusterAlignment {
    pub fn from_read_and_source(read: &Record, taxid: u32, source: &str) -> Result<Self, Error> {
        let hash = ReadHash::from_read_and_source(read, source);
        let map_interval = Range { start: read.pos(), end: read.reference_end() - 1 };
        let cover = filter_read(read)?;
        let seq_len = read.seq_len();
        Ok(Self {
            hash,
            name: read.qname().to_vec(),
            forward: !read.is_reverse(),
            tid: read.tid(),
            seq_len,
            taxid,
            cover,
            map_interval,
        })
    }
}

pub struct ReadBucket {
    pub rec: ClusterAlignment,
    pub mate: Option<ClusterAlignment>,
}

#[derive(Default)]
pub struct SamMergeBuffer {
    q: Vec<ReadBucket>,
    read_id: u64,
}

impl SamMergeBuffer {
    pub fn intake(&mut self, aln: ClusterAlignment) -> Option<Vec<ReadBucket>> {
        if self.read_id == aln.hash.name_hash {
            // already encountered this mate. So add and keep going.
            self.update_mate(aln);
            None
        } else {
            // we encountered a read QNAME we haven't seen before. Assuming a name-sorted BAM file,
            // this means we can emit the reads we've been holding in the queue.
            let ret = std::mem::take(&mut self.q);
            self.read_id = aln.hash.name_hash;
            self.q.push(ReadBucket { rec: aln, mate: None });
            Some(ret)
        }
    }

    pub fn update_mate(&mut self, aln: ClusterAlignment) {
        for node in self.q.iter_mut() {
            if node.rec.hash == aln.hash {
                assert!(node.mate.is_none());
                node.mate = Some(aln);
                return;
            }
        }

        // no mate found with same name and source, so update queue.
        self.q.push(ReadBucket { rec: aln, mate: None });
    }

    pub fn flush(&mut self) -> Vec<ReadBucket> {
        std::mem::take(&mut self.q)
    }
}
