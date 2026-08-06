use std::{
    collections::VecDeque,
    fs::ReadDir,
    path::{Path, PathBuf},
};

use anyhow::{Context, Error};

pub struct AssemblyDirIterator {
    q: VecDeque<(ReadDir, Option<String>)>,
}

impl AssemblyDirIterator {
    pub fn new<P: AsRef<Path>>(input: P) -> Result<Self, Error> {
        let input = input.as_ref();
        let assembly_name = if is_assembly_dir(input)? {
            file_name(input)?.map(str::to_owned)
        } else {
            None
        };

        let mut q = VecDeque::new();
        q.push_back((
            std::fs::read_dir(input)
                .with_context(|| format!("couldn't read input directory {:?}", input))?,
            assembly_name,
        ));

        Ok(Self { q })
    }

    pub fn next_item(&mut self) -> Result<Option<(String, String)>, Error> {
        loop {
            let (entry, assembly_name) = {
                let Some((dir, assembly_name)) = self.q.front_mut() else {
                    return Ok(None);
                };
                let Some(entry) = dir.next() else {
                    self.q.pop_front();
                    continue;
                };
                (entry?, assembly_name.clone())
            };

            let path = entry.path();
            let file_type = entry.file_type()?;

            if file_type.is_dir() {
                let child_assembly_name = if is_assembly_dir(&path)? {
                    file_name(&path)?.map(str::to_owned)
                } else {
                    assembly_name
                };
                self.q
                    .push_back((std::fs::read_dir(&path)?, child_assembly_name));
                continue;
            }

            if file_type.is_file()
                && file_name(&path)?.is_some_and(|name| name.contains(".fna"))
                && let Some(assembly_name) = assembly_name
            {
                return Ok(Some((assembly_name, path_to_string(path)?)));
            }
        }
    }
}

fn is_assembly_dir(path: &Path) -> Result<bool, Error> {
    Ok(file_name(path)?.is_some_and(|name| name.starts_with("GCF_") || name.starts_with("GCA_")))
}

fn file_name(p: &Path) -> Result<Option<&str>, Error> {
    p.file_name()
        .map(|fname| {
            fname
                .to_str()
                .with_context(|| format!("couldn't get file name for {:?}", fname))
        })
        .transpose()
}

impl Iterator for AssemblyDirIterator {
    type Item = Result<(String, String), Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_item().transpose()
    }
}

fn path_to_string(path: PathBuf) -> Result<String, Error> {
    path.into_os_string()
        .into_string()
        .map_err(|path| anyhow::anyhow!("couldn't convert path {:?} to UTF-8", path))
}

#[cfg(test)]
mod tests {
    use super::AssemblyDirIterator;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn finds_all_fna_files_below_nested_assembly_directories() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!("bamtax-assembly-iterator-{unique}"));
        let assembly = root.join("download").join("GCF_000001.1_example");
        let nested = assembly.join("nested");
        let unrelated = root.join("download").join("other");

        fs::create_dir_all(&nested).unwrap();
        fs::create_dir_all(&unrelated).unwrap();
        fs::write(assembly.join("genomic.fna.gz"), []).unwrap();
        fs::write(nested.join("extra.fna"), []).unwrap();
        fs::write(assembly.join("notes.txt"), []).unwrap();
        fs::write(unrelated.join("ignored.fna"), []).unwrap();

        let mut found = AssemblyDirIterator::new(&root)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        found.sort();

        assert_eq!(
            found,
            vec![
                (
                    "GCF_000001.1_example".to_owned(),
                    assembly.join("genomic.fna.gz").to_str().unwrap().to_owned(),
                ),
                (
                    "GCF_000001.1_example".to_owned(),
                    nested.join("extra.fna").to_str().unwrap().to_owned(),
                ),
            ]
        );

        fs::remove_dir_all(root).unwrap();
    }
}
