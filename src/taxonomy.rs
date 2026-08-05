use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

#[derive(Debug, PartialEq, Eq)]
pub struct Taxon {
    pub tax_id: u32,
    pub parent_tax_id: u32,
    pub rank: String,
    pub name: String,
}

#[derive(Debug)]
pub struct Taxonomy {
    nodes: HashMap<u32, Taxon>,
    memo: HashMap<u32, Option<u32>>,
}

impl Taxonomy {
    pub fn from_dir(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let nodes_path = path.join("nodes.dmp");
        let names_path = path.join("names.dmp");
        let nodes =
            File::open(&nodes_path).with_context(|| format!("open {}", nodes_path.display()))?;
        let names =
            File::open(&names_path).with_context(|| format!("open {}", names_path.display()))?;

        Self::from_readers(BufReader::new(nodes), BufReader::new(names))
    }

    pub fn get(&self, tax_id: u32) -> Option<&Taxon> {
        self.nodes.get(&tax_id)
    }

    /// Returns the lineage from the root through the requested taxon.
    pub fn lineage(&self, tax_id: u32) -> Option<Vec<&Taxon>> {
        let mut lineage = Vec::new();
        let mut seen = HashSet::new();
        let mut current = tax_id;

        loop {
            if !seen.insert(current) {
                return None;
            }

            let taxon = self.nodes.get(&current)?;
            lineage.push(taxon);
            if taxon.parent_tax_id == current {
                break;
            }
            current = taxon.parent_tax_id;
        }

        lineage.reverse();
        Some(lineage)
    }

    /// Finds the requested taxon's closest ancestor with rank `species`.
    pub fn species(&mut self, tax_id: u32) -> Option<&Taxon> {
        if let Some(memo) = self.memo.get(&tax_id) {
            return memo.and_then(|t| self.get(t));
        }

        let mut seen = HashSet::new();
        let mut current = tax_id;

        loop {
            if !seen.insert(current) {
                self.memo.insert(tax_id, None);
                return None;
            }

            let taxon = self.nodes.get(&current)?;

            if taxon.rank == "species" {
                self.memo.insert(tax_id, Some(taxon.tax_id));
                return Some(taxon);
            }
            if taxon.parent_tax_id == current {
                self.memo.insert(tax_id, None);
                return None;
            }
            current = taxon.parent_tax_id;
        }
    }

    fn from_readers(nodes: impl BufRead, names: impl BufRead) -> Result<Self> {
        let mut taxonomy = Self {
            nodes: HashMap::new(),
            memo: HashMap::new(),
        };

        for (line_number, line) in nodes.lines().enumerate() {
            let line = line.with_context(|| format!("read nodes.dmp line {}", line_number + 1))?;
            let mut fields = line.split("\t|\t");
            let tax_id = parse_tax_id(fields.next(), "tax_id", "nodes.dmp", line_number)?;
            let parent_tax_id =
                parse_tax_id(fields.next(), "parent tax_id", "nodes.dmp", line_number)?;
            let rank = fields
                .next()
                .with_context(|| format!("nodes.dmp line {} has no rank", line_number + 1))?;

            taxonomy.nodes.insert(
                tax_id,
                Taxon {
                    tax_id,
                    parent_tax_id,
                    rank: rank.to_owned(),
                    name: String::new(),
                },
            );
        }

        for (line_number, line) in names.lines().enumerate() {
            let line = line.with_context(|| format!("read names.dmp line {}", line_number + 1))?;
            let mut fields = line.split("\t|\t");
            let tax_id = parse_tax_id(fields.next(), "tax_id", "names.dmp", line_number)?;
            let name = fields
                .next()
                .with_context(|| format!("names.dmp line {} has no name", line_number + 1))?;
            let _unique_name = fields.next();
            let name_class = fields
                .next()
                .with_context(|| format!("names.dmp line {} has no name class", line_number + 1))?
                .trim_end_matches("\t|");

            if name_class == "scientific name"
                && let Some(taxon) = taxonomy.nodes.get_mut(&tax_id)
            {
                taxon.name = name.to_owned();
            }
        }

        Ok(taxonomy)
    }
}

fn parse_tax_id(
    field: Option<&str>,
    field_name: &str,
    file_name: &str,
    zero_based_line_number: usize,
) -> Result<u32> {
    let line_number = zero_based_line_number + 1;
    field
        .with_context(|| format!("{file_name} line {line_number} has no {field_name}"))?
        .parse()
        .with_context(|| format!("invalid {field_name} in {file_name} line {line_number}"))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn taxonomy() -> Taxonomy {
        let nodes = concat!(
            "1\t|\t1\t|\tno rank\t|\tignored\t|\n",
            "2\t|\t1\t|\tsuperkingdom\t|\tignored\t|\n",
            "10\t|\t2\t|\tspecies\t|\tignored\t|\n",
            "11\t|\t10\t|\tstrain\t|\tignored\t|\n",
        );
        let names = concat!(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n",
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n",
            "10\t|\tExample common name\t|\t\t|\tcommon name\t|\n",
            "10\t|\tExample species\t|\t\t|\tscientific name\t|\n",
            "11\t|\tExample strain\t|\t\t|\tscientific name\t|\n",
        );

        Taxonomy::from_readers(Cursor::new(nodes), Cursor::new(names)).unwrap()
    }

    #[test]
    fn gets_taxon_by_tax_id() {
        let taxonomy = taxonomy();

        assert_eq!(taxonomy.get(10).unwrap().name, "Example species");
        assert!(taxonomy.get(999).is_none());
    }

    #[test]
    fn returns_root_to_taxon_lineage() {
        let taxonomy = taxonomy();
        let lineage = taxonomy.lineage(11).unwrap();
        let tax_ids: Vec<_> = lineage.iter().map(|taxon| taxon.tax_id).collect();

        assert_eq!(tax_ids, vec![1, 2, 10, 11]);
        assert!(taxonomy.lineage(999).is_none());
    }

    #[test]
    fn resolves_species_for_descendant_tax_id() {
        let mut taxonomy = taxonomy();

        assert_eq!(taxonomy.species(11).unwrap().tax_id, 10);
        assert_eq!(taxonomy.species(10).unwrap().name, "Example species");
        assert!(taxonomy.species(2).is_none());
    }
}
