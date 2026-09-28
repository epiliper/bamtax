use anyhow::{Context, Error, bail};

const TID_INDIC: &str = "TID_";

#[repr(usize)]
pub enum K2TxtCol {
    IsClassified = 0,
    ReadHeader = 1,
    TaxonId = 2,
    SequenceLength = 3,
    KmerIndices = 4,
}

impl TryFrom<usize> for K2TxtCol {
    type Error = Error;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::IsClassified),
            1 => Ok(Self::ReadHeader),
            2 => Ok(Self::TaxonId),
            3 => Ok(Self::SequenceLength),
            4 => Ok(Self::KmerIndices),
            other => bail!("Anomalous number of report columns {other}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct K2ReportRow {
    pub fraction: f32,
    pub total_reads: u64,
    pub exact_reads: u64,
    pub rank: String,
    pub taxid: u32,
    pub name: String,
    pub depth: usize,
}

impl TryFrom<&str> for K2ReportRow {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        const NCOLS: usize = 6;

        let fields = value.split('\t').collect::<Vec<_>>();
        if fields.len() != NCOLS {
            bail!(
                "Invalid K2 report; anomalous number of columns (tab-delimited): {}, expecting {NCOLS}",
                fields.len()
            );
        }

        let indented_name = fields[5].trim_end();
        let indentation = indented_name.len() - indented_name.trim_start_matches(' ').len();
        if indentation % 2 != 0 {
            bail!("Invalid K2 report indentation: expected two spaces per level");
        }

        Ok(Self {
            fraction: fields[0].trim().parse().context("invalid report fraction")?,
            total_reads: fields[1].trim().parse().context("invalid report clade read count")?,
            exact_reads: fields[2].trim().parse().context("invalid report exact read count")?,
            rank: fields[3].trim().to_owned(),
            taxid: fields[4].trim().parse().context("invalid report taxid")?,
            name: indented_name[indentation..].to_owned(),
            depth: indentation / 2,
        })
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum K2ClassificationStatus {
    Unclassified,
    Classified,
}

impl TryFrom<char> for K2ClassificationStatus {
    type Error = Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'U' => Ok(K2ClassificationStatus::Unclassified),
            'C' => Ok(K2ClassificationStatus::Classified),
            _ => bail!(
                "Invalid character string for classification status (column {}): {}",
                K2TxtCol::IsClassified as usize,
                value
            ),
        }
    }
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone)]
pub enum K2HitType {
    TaxonId(u32),
    NotInDatabase,
    AmbiguousBase,
}

impl TryFrom<&str> for K2HitType {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value == "A" {
            return Ok(Self::AmbiguousBase);
        }

        let intvalue = value.parse::<u32>()?;
        if intvalue == 0 {
            return Ok(Self::NotInDatabase);
        }

        Ok(Self::TaxonId(intvalue))
    }
}

impl std::fmt::Display for K2HitType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TaxonId(tid) => write!(f, "TaxonID {tid}"),
            Self::NotInDatabase => write!(f, "Not in DB"),
            Self::AmbiguousBase => write!(f, "Ambig base"),
        }
    }
}

// Describes one taxon/not-classified/ambiguous run in Kraken's hit list.
#[derive(Debug, Eq, PartialEq, PartialOrd, Ord, Clone)]
pub struct K2HitDesc {
    pub hit: K2HitType,
    pub n_kmers: usize,
}

impl TryFrom<&str> for K2HitDesc {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if let Some((tid, n_kmers)) = value.trim().split_once(':')
            && let (Ok(hit), Ok(n_kmers)) = (K2HitType::try_from(tid), n_kmers.parse::<usize>())
        {
            return Ok(K2HitDesc { hit, n_kmers });
        }

        bail!("Invalid tid/kmer string: {value}");
    }
}

impl std::fmt::Display for K2HitDesc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.hit {
            K2HitType::TaxonId(id) => write!(f, "{TID_INDIC}{id}"),
            K2HitType::NotInDatabase => write!(f, "not_in_database"),
            K2HitType::AmbiguousBase => write!(f, "ambiguous_base"),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct K2ClassificationRow {
    pub status: K2ClassificationStatus,
    pub read_name: String,
    pub assigned_taxid: u32,
    pub hits: Vec<K2HitDesc>,
}

impl TryFrom<&str> for K2ClassificationRow {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        const NCOLS: usize = 5;

        let fields = value.split('\t').collect::<Vec<_>>();
        if fields.len() != NCOLS {
            bail!(
                "Invalid K2 classification row; anomalous number of columns (tab-delimited): {}, expecting {NCOLS}",
                fields.len()
            );
        }

        let mut status_chars = fields[0].chars();
        let status = status_chars.next().context("missing classification status")?;
        if status_chars.next().is_some() {
            bail!("invalid classification status: {}", fields[0]);
        }

        let hits = fields[4]
            .split_whitespace()
            .filter(|token| *token != "|:|")
            .map(K2HitDesc::try_from)
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            status: K2ClassificationStatus::try_from(status)?,
            read_name: fields[1].to_owned(),
            assigned_taxid: fields[2].parse().context("invalid classification taxid")?,
            hits,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_report_metadata_and_indentation() {
        let row = K2ReportRow::try_from(" 50.00\t10\t4\tS\t123\t    Species name  ").unwrap();

        assert_eq!(row.fraction, 50.0);
        assert_eq!(row.total_reads, 10);
        assert_eq!(row.exact_reads, 4);
        assert_eq!(row.rank, "S");
        assert_eq!(row.taxid, 123);
        assert_eq!(row.name, "Species name");
        assert_eq!(row.depth, 2);
    }

    #[test]
    fn rejects_malformed_report_rows() {
        assert!(K2ReportRow::try_from("1\t2\t3").is_err());
        assert!(K2ReportRow::try_from("1\t2\t3\tS\t4\t odd").is_err());
    }

    #[test]
    fn parses_classification_and_paired_separator() {
        let row = K2ClassificationRow::try_from("C\tread-1\t123\t100|100\t123:4 0:2 |:| A:1 456:3").unwrap();

        assert_eq!(row.status, K2ClassificationStatus::Classified);
        assert_eq!(row.read_name, "read-1");
        assert_eq!(row.assigned_taxid, 123);
        assert_eq!(row.hits.len(), 4);
        assert_eq!(row.hits[0].hit, K2HitType::TaxonId(123));
        assert_eq!(row.hits[3].hit, K2HitType::TaxonId(456));
    }
}
