use clap::ValueEnum;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use std::cell::RefCell;

#[derive(Debug, Eq)]
pub struct Taxon {
    pub tax_id: u32,
    pub parent_tax_id: u32,
    pub rank: String,
    pub name: String,
}

impl PartialEq for Taxon {
    fn eq(&self, other: &Self) -> bool {
        self.tax_id == other.tax_id
    }
}

impl std::hash::Hash for Taxon {
    fn hash<H>(&self, state: &mut H)
    where
        H: std::hash::Hasher,
    {
        self.tax_id.hash(state)
    }
}

#[derive(Debug, Copy, Clone, ValueEnum)]
pub enum Rank {
    Species,
    Genus,
}

#[derive(Debug)]
// Do not use with multiple threads.
pub struct Taxonomy {
    nodes: HashMap<u32, Taxon>,
    species_memo: RefCell<HashMap<u32, Option<u32>>>,
    genus_memo: RefCell<HashMap<u32, Option<u32>>>,
    children: HashMap<u32, Vec<u32>>,
}

impl<'a> Taxonomy {
    pub fn from_dir(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let nodes_path = path.join("nodes.dmp");
        let names_path = path.join("names.dmp");
        let nodes = File::open(&nodes_path).with_context(|| format!("open {}", nodes_path.display()))?;
        let names = File::open(&names_path).with_context(|| format!("open {}", names_path.display()))?;

        Self::from_readers(BufReader::new(nodes), BufReader::new(names))
    }

    pub fn get(&self, tax_id: u32) -> Option<&Taxon> {
        self.nodes.get(&tax_id)
    }

    /// Returns the lineage from the root through the requested taxon.
    pub fn lineage(&'a self, tax_id: u32) -> Option<Vec<&'a Taxon>> {
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

    pub fn species(&self, tax_id: u32) -> Option<&Taxon> {
        self.lookup(tax_id, Rank::Species)
    }

    pub fn genus(&self, tax_id: u32) -> Option<&Taxon> {
        self.lookup(tax_id, Rank::Genus)
    }

    pub fn lookup(&self, tax_id: u32, rank: Rank) -> Option<&Taxon> {
        let (mut memo, name) = match rank {
            Rank::Species => (self.species_memo.borrow_mut(), "species"),
            Rank::Genus => (self.genus_memo.borrow_mut(), "genus"),
        };

        let mut seen = HashSet::new();
        let mut current = tax_id;

        loop {
            if !seen.insert(current) {
                memo.insert(tax_id, None);
                return None;
            }

            let taxon = self.nodes.get(&current)?;

            if taxon.rank == name {
                memo.insert(tax_id, Some(taxon.tax_id));
                return Some(taxon);
            }
            if taxon.parent_tax_id == current {
                memo.insert(tax_id, None);
                return None;
            }
            current = taxon.parent_tax_id;
        }
    }

    pub fn descendants(&self, taxid: u32) -> impl Iterator<Item = &[u32]> {
        let mut deque: VecDeque<u32> = VecDeque::from([taxid]);
        let mut seen: HashSet<u32> = HashSet::new();

        std::iter::from_fn(move || {
            while let Some(cur) = deque.pop_front() {
                if seen.insert(cur) {
                    let children = self.children.get(&cur).map(|v| v.as_slice()).unwrap_or(&[]);
                    deque.extend(children);
                    return Some(children);
                }
            }
            None
        })
    }

    fn from_readers(nodes: impl BufRead, names: impl BufRead) -> Result<Self> {
        let mut taxonomy = Self {
            nodes: HashMap::new(),
            species_memo: RefCell::new(HashMap::new()),
            genus_memo: RefCell::new(HashMap::new()),
            children: HashMap::new(),
        };

        for (line_number, line) in nodes.lines().enumerate() {
            let line = line.with_context(|| format!("read nodes.dmp line {}", line_number + 1))?;
            let mut fields = line.split("\t|\t");
            let tax_id = parse_tax_id(fields.next(), "tax_id", "nodes.dmp", line_number)?;
            let parent_tax_id = parse_tax_id(fields.next(), "parent tax_id", "nodes.dmp", line_number)?;
            let rank = fields.next().with_context(|| format!("nodes.dmp line {} has no rank", line_number + 1))?;

            taxonomy.nodes.insert(tax_id, Taxon { tax_id, parent_tax_id, rank: rank.to_owned(), name: String::new() });
        }

        for (line_number, line) in names.lines().enumerate() {
            let line = line.with_context(|| format!("read names.dmp line {}", line_number + 1))?;
            let mut fields = line.split("\t|\t");
            let tax_id = parse_tax_id(fields.next(), "tax_id", "names.dmp", line_number)?;
            let name = fields.next().with_context(|| format!("names.dmp line {} has no name", line_number + 1))?;
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

        // last iteration to build child mappings
        for (k, v) in taxonomy.nodes.iter() {
            if v.tax_id == v.parent_tax_id {
                continue;
            }

            taxonomy.children.entry(v.parent_tax_id).or_default().push(*k);
        }

        Ok(taxonomy)
    }
}

fn parse_tax_id(field: Option<&str>, field_name: &str, file_name: &str, zero_based_line_number: usize) -> Result<u32> {
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
