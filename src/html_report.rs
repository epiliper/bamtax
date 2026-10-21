use crate::locus_tracker::AlignmentReport;
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Result};
use rust_htslib::bam::HeaderView;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

const REPORT_TEMPLATE: &str = include_str!("html_report/report.html");
const REPORT_STYLE: &str = include_str!("html_report/report.css");
const REPORT_SCRIPT: &str = include_str!("html_report/report.js");
const KRONA_BRIDGE_SCRIPT: &str = include_str!("html_report/krona_bridge.js");

const STYLE_SLOT: &str = "<!-- BAMTAX:STYLE -->";
const DATASET_OPTIONS_SLOT: &str = "<!-- BAMTAX:DATASET_OPTIONS -->";
const TAXON_PANELS_SLOT: &str = "<!-- BAMTAX:TAXON_PANELS -->";
const KRONA_SLOT: &str = "<!-- BAMTAX:KRONA -->";
const SCRIPT_SLOT: &str = "<!-- BAMTAX:SCRIPT -->";

pub struct HtmlSample {
    name: String,
    extract_source: String,
    species: Vec<SpeciesSummary>,
}

struct SpeciesSummary {
    tax_id: u32,
    depth: usize,
    loci: usize,
    references: Vec<String>,
}

#[derive(Clone)]
struct ExtractCall {
    source: String,
    references: String,
}

#[derive(Default)]
struct TreeNode {
    tax_id: Option<u32>,
    name: String,
    rank: String,
    depth: usize,
    loci: usize,
    references: BTreeSet<String>,
    samples: BTreeSet<usize>,
    extract_calls: Vec<ExtractCall>,
    children: BTreeMap<Option<u32>, TreeNode>,
}

impl HtmlSample {
    pub fn from_report(name: &str, source: &str, report: &AlignmentReport, header: &HeaderView) -> Self {
        let mut species = report
            .alns
            .iter()
            .map(|(tax_id, (alignments, tids))| {
                let mut references: Vec<_> =
                    tids.iter().map(|tid| String::from_utf8_lossy(header.tid2name(*tid as u32)).into_owned()).collect();
                references.sort();

                SpeciesSummary {
                    tax_id: *tax_id,
                    depth: alignments.iter().map(|alignment| alignment.depth).sum(),
                    loci: alignments.len(),
                    references,
                }
            })
            .collect::<Vec<_>>();
        species.sort_by_key(|summary| summary.tax_id);

        Self { name: name.to_owned(), extract_source: absolute_source(source), species }
    }

    pub fn empty(name: &str, source: &str) -> Self {
        Self { name: name.to_owned(), extract_source: absolute_source(source), species: Vec::new() }
    }
}

pub fn write_html_report(
    path: impl AsRef<Path>,
    samples: &[HtmlSample],
    taxonomy: &Taxonomy,
    krona_command: &str,
    krona_taxonomy: Option<&Path>,
) -> Result<()> {
    let path = path.as_ref();
    let krona = if samples.iter().any(|sample| !sample.species.is_empty()) {
        generate_krona(path, samples, krona_command, krona_taxonomy)
    } else {
        Err("No Krona chart was generated because no taxa passed the report filters.".to_owned())
    };
    let html = render_report(samples, taxonomy, krona.as_ref().map(String::as_str).map_err(String::as_str));
    fs::write(path, html).with_context(|| format!("write HTML report {}", path.display()))?;

    if let Err(error) = krona {
        eprintln!("Warning: {error}");
    }

    Ok(())
}

fn generate_krona(
    report_path: &Path,
    samples: &[HtmlSample],
    krona_command: &str,
    krona_taxonomy: Option<&Path>,
) -> Result<String, String> {
    let temporary_dir = create_temporary_dir(report_path)?;
    let result = (|| {
        let executable = if Path::new(krona_command).components().count() > 1 {
            fs::canonicalize(krona_command)
                .map_err(|error| format!("could not resolve Krona command {krona_command}: {error}"))?
        } else {
            PathBuf::from(krona_command)
        };
        let mut command = Command::new(executable);
        command.current_dir(&temporary_dir).args(["-q", "0", "-t", "1", "-m", "2", "-o", "chart.html"]);

        if let Some(path) = krona_taxonomy {
            let path = fs::canonicalize(path)
                .map_err(|error| format!("could not resolve Krona taxonomy {}: {error}", path.display()))?;
            command.arg("-tax").arg(path);
        }

        let mut aggregate = BTreeMap::<u32, usize>::new();
        for sample in samples {
            for species in &sample.species {
                *aggregate.entry(species.tax_id).or_default() += species.depth;
            }
        }

        let aggregate_input = aggregate.into_iter().collect::<Vec<_>>();
        let datasets = std::iter::once(("All samples", aggregate_input)).chain(samples.iter().map(|sample| {
            (
                sample.name.as_str(),
                sample.species.iter().map(|species| (species.tax_id, species.depth)).collect::<Vec<_>>(),
            )
        }));
        for (index, (name, species)) in datasets.enumerate() {
            let filename = format!("dataset-{index}.tsv");
            let mut input = String::new();
            for (tax_id, depth) in species {
                let _ = writeln!(input, "{tax_id}\t{depth}");
            }
            fs::write(temporary_dir.join(&filename), input)
                .map_err(|error| format!("could not write temporary Krona input: {error}"))?;
            let dataset_name = name.replace([',', '"', '\n', '\r'], " ");
            command.arg(format!("{filename},\"{dataset_name}\""));
        }

        let output = command.output().map_err(|error| format!("could not run {krona_command}: {error}"))?;
        let stderr = String::from_utf8_lossy(&output.stderr);
        if !output.status.success() {
            let detail = stderr.trim();
            return Err(if detail.is_empty() {
                format!("{krona_command} exited with status {}", output.status)
            } else {
                format!("{krona_command} failed: {detail}")
            });
        }
        if !stderr.trim().is_empty() {
            eprintln!("Warning from {krona_command}: {}", stderr.trim());
        }

        fs::read_to_string(temporary_dir.join("chart.html"))
            .map(inject_krona_bridge)
            .map_err(|error| format!("could not read generated Krona chart: {error}"))
    })();
    if let Err(error) = fs::remove_dir_all(temporary_dir) {
        eprintln!("Warning: could not remove temporary Krona files: {error}");
    }
    result
}

fn inject_krona_bridge(mut html: String) -> String {
    let bridge = format!("<script>\n{KRONA_BRIDGE_SCRIPT}\n</script>");
    if let Some(index) = html.rfind("</body>") {
        html.insert_str(index, &bridge);
    } else {
        html.push_str(&bridge);
    }
    html
}

fn create_temporary_dir(report_path: &Path) -> Result<PathBuf, String> {
    let stem = report_path.file_stem().and_then(|stem| stem.to_str()).unwrap_or("report");
    let stem = stem.chars().filter(|character| character.is_ascii_alphanumeric()).take(24).collect::<String>();
    for attempt in 0..100_u8 {
        let path = std::env::temp_dir().join(format!("bamtax-{stem}-{}-{attempt}", std::process::id()));
        match fs::create_dir(&path) {
            Ok(()) => {
                #[cfg(unix)]
                if let Err(error) = fs::set_permissions(&path, fs::Permissions::from_mode(0o700)) {
                    let _ = fs::remove_dir(&path);
                    return Err(format!("could not secure temporary Krona directory: {error}"));
                }
                return Ok(path);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(format!("could not create temporary Krona directory: {error}")),
        }
    }
    Err("could not create a unique temporary Krona directory".to_owned())
}

fn render_report(samples: &[HtmlSample], taxonomy: &Taxonomy, krona: Result<&str, &str>) -> String {
    let mut dataset_options = "          <option value=\"panel-0\">All samples</option>\n".to_owned();
    for (index, sample) in samples.iter().enumerate() {
        let _ = writeln!(
            dataset_options,
            "          <option value=\"panel-{}\">{}</option>",
            index + 1,
            escape_html(&sample.name)
        );
    }

    let mut taxon_panels = String::new();
    render_panel(&mut taxon_panels, "panel-0", samples.iter(), taxonomy, false);
    for (index, sample) in samples.iter().enumerate() {
        render_panel(&mut taxon_panels, &format!("panel-{}", index + 1), std::iter::once(sample), taxonomy, true);
    }

    let krona = match krona {
        Ok(krona_html) => format!(
            "    <iframe title=\"Krona taxonomy chart\" sandbox=\"allow-scripts allow-popups allow-downloads\" srcdoc=\"{}\"></iframe>",
            escape_html(krona_html)
        ),
        Err(error) => format!("    <div class=\"notice\">{}</div>", escape_html(error)),
    };

    REPORT_TEMPLATE
        .replace(STYLE_SLOT, REPORT_STYLE)
        .replace(DATASET_OPTIONS_SLOT, &dataset_options)
        .replace(TAXON_PANELS_SLOT, &taxon_panels)
        .replace(KRONA_SLOT, &krona)
        .replace(SCRIPT_SLOT, REPORT_SCRIPT)
}

fn render_panel<'a>(
    html: &mut String,
    id: &str,
    samples: impl Iterator<Item = &'a HtmlSample>,
    taxonomy: &Taxonomy,
    hidden: bool,
) {
    let samples = samples.collect::<Vec<_>>();
    let tree = build_tree(&samples, taxonomy);
    let hidden = if hidden { " hidden" } else { "" };
    let _ = writeln!(html, "    <div class=\"panel\" id=\"{}\"{}>", id, hidden);
    html.push_str(
        "      <div class=\"tree-scroll\"><div class=\"tree\">\n        <div class=\"row header\"><span>Taxon</span><span>Rank</span><span>Alignments</span><span>Loci</span><span>Samples</span><span>References</span><span>Action</span></div>\n",
    );
    if tree.children.is_empty() {
        html.push_str("        <div class=\"empty\">No taxa passed the report filters.</div>\n");
    } else {
        let mut nodes = tree.children.values().collect::<Vec<_>>();
        sort_nodes(&mut nodes);
        for node in nodes {
            render_node(html, node, 0);
        }
    }
    html.push_str("      </div></div>\n    </div>\n");
}

fn build_tree(samples: &[&HtmlSample], taxonomy: &Taxonomy) -> TreeNode {
    let mut root = TreeNode::default();
    for (sample_index, sample) in samples.iter().enumerate() {
        for species in &sample.species {
            let species_node = taxonomy.get(species.tax_id);
            let species_name = species_node
                .filter(|taxon| !taxon.name.is_empty())
                .map(|taxon| taxon.name.clone())
                .unwrap_or_else(|| species.tax_id.to_string());
            let extract_call =
                ExtractCall { source: sample.extract_source.clone(), references: species.references.join(";") };
            let lineage = taxonomy.lineage(species.tax_id).unwrap_or_default();
            let path = if let Some(family) = lineage.iter().position(|taxon| taxon.rank == "family") {
                lineage[family..]
                    .iter()
                    .map(|taxon| NodeInfo {
                        tax_id: Some(taxon.tax_id),
                        name: if taxon.name.is_empty() { taxon.tax_id.to_string() } else { taxon.name.clone() },
                        rank: taxon.rank.clone(),
                    })
                    .collect::<Vec<_>>()
            } else {
                vec![
                    NodeInfo { tax_id: None, name: "Unclassified family".to_owned(), rank: "family".to_owned() },
                    NodeInfo {
                        tax_id: Some(species.tax_id),
                        name: species_name,
                        rank: species_node.map(|taxon| taxon.rank.clone()).unwrap_or_else(|| "species".to_owned()),
                    },
                ]
            };

            let mut node = &mut root;
            for taxon in path {
                node = node.children.entry(taxon.tax_id).or_insert_with(|| TreeNode {
                    tax_id: taxon.tax_id,
                    name: taxon.name,
                    rank: taxon.rank,
                    ..TreeNode::default()
                });
                node.depth += species.depth;
                node.loci += species.loci;
                node.references.extend(species.references.iter().cloned());
                node.samples.insert(sample_index);
                node.extract_calls.push(extract_call.clone());
            }
        }
    }
    root
}

struct NodeInfo {
    tax_id: Option<u32>,
    name: String,
    rank: String,
}

fn sort_nodes(nodes: &mut Vec<&TreeNode>) {
    nodes.sort_by(|a, b| b.depth.cmp(&a.depth).then_with(|| a.name.cmp(&b.name)).then_with(|| a.tax_id.cmp(&b.tax_id)));
}

fn render_node(html: &mut String, node: &TreeNode, depth: usize) {
    let species_class = if node.rank == "species" { " species-call" } else { "" };
    if node.children.is_empty() {
        let _ = write!(html, "        <div class=\"row leaf{}\" style=\"--depth:{}\">", species_class, depth);
        render_cells(html, node, false);
        html.push_str("</div>\n");
        return;
    }

    let _ = write!(html, "        <details><summary class=\"row{}\" style=\"--depth:{}\">", species_class, depth);
    render_cells(html, node, true);
    html.push_str("</summary><div class=\"children\">\n");
    let mut children = node.children.values().collect::<Vec<_>>();
    sort_nodes(&mut children);
    for child in children {
        render_node(html, child, depth + 1);
    }
    html.push_str("        </div></details>\n");
}

fn render_cells(html: &mut String, node: &TreeNode, expandable: bool) {
    let tax_id = node.tax_id.map(|id| format!("taxid {id}")).unwrap_or_else(|| "synthetic group".to_owned());
    let references = node.references.iter().cloned().collect::<Vec<_>>().join("\n");
    let command = extract_command(&node.extract_calls);
    let toggle = if expandable { "<span class=\"toggle\" aria-hidden=\"true\">▸</span>" } else { "" };
    let _ = write!(
        html,
        "<span class=\"taxon-cell\">{}<span class=\"taxon\">{}</span><br><span class=\"meta\">{}</span></span><span>{}</span><span>{}</span><span>{}</span><span>{}</span><button type=\"button\" class=\"references-button\" data-taxon=\"{}\" data-references=\"{}\" aria-label=\"Show {} references for {}\">{}</button><button type=\"button\" class=\"extract-button\" data-command=\"{}\" title=\"Copy command for {} species call{}\">Extract</button>",
        toggle,
        escape_html(&node.name),
        tax_id,
        escape_html(&node.rank),
        node.depth,
        node.loci,
        node.samples.len(),
        escape_html(&node.name),
        escape_html(&references),
        node.references.len(),
        escape_html(&node.name),
        node.references.len(),
        escape_html(&command),
        node.extract_calls.len(),
        if node.extract_calls.len() == 1 { "" } else { "s" },
    );
}

fn extract_command(calls: &[ExtractCall]) -> String {
    calls
        .iter()
        .map(|call| {
            let references = call
                .references
                .split(';')
                .map(|reference| format!("keep[\"{}\"] = 1", escape_awk_string(reference)))
                .collect::<Vec<_>>()
                .join("; ");
            let program = format!("BEGIN {{ {references} }} $3 in keep");
            format!("samtools view {} | awk -F '\\t' {}", shell_quote(&call.source), shell_quote(&program))
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn escape_awk_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '\\' => escaped.push_str("\\\\"),
            '"' => escaped.push_str("\\\""),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            _ => escaped.push(character),
        }
    }
    escaped
}

fn shell_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

fn absolute_source(source: &str) -> String {
    let path = Path::new(source);
    if path.is_absolute() {
        source.to_owned()
    } else {
        std::env::current_dir()
            .map(|directory| directory.join(path).to_string_lossy().into_owned())
            .unwrap_or_else(|_| source.to_owned())
    }
}

fn escape_html(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn taxonomy() -> Taxonomy {
        let nodes = concat!(
            "1\t|\t1\t|\tno rank\t|\tignored\t|\n",
            "2\t|\t1\t|\tsuperkingdom\t|\tignored\t|\n",
            "10\t|\t2\t|\tfamily\t|\tignored\t|\n",
            "11\t|\t10\t|\tno rank\t|\tignored\t|\n",
            "20\t|\t11\t|\tgenus\t|\tignored\t|\n",
            "30\t|\t20\t|\tspecies\t|\tignored\t|\n",
            "40\t|\t2\t|\tfamily\t|\tignored\t|\n",
            "50\t|\t40\t|\tspecies\t|\tignored\t|\n",
            "60\t|\t2\t|\tspecies\t|\tignored\t|\n",
        );
        let names = concat!(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n",
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n",
            "10\t|\tAlpha family\t|\t\t|\tscientific name\t|\n",
            "11\t|\tAlpha group\t|\t\t|\tscientific name\t|\n",
            "20\t|\tAlpha genus\t|\t\t|\tscientific name\t|\n",
            "30\t|\tAlpha <species>\t|\t\t|\tscientific name\t|\n",
            "40\t|\tBeta family\t|\t\t|\tscientific name\t|\n",
            "50\t|\tBeta species\t|\t\t|\tscientific name\t|\n",
            "60\t|\tOrphan species\t|\t\t|\tscientific name\t|\n",
        );
        Taxonomy::from_readers(Cursor::new(nodes), Cursor::new(names)).unwrap()
    }

    fn sample(name: &str, species: Vec<SpeciesSummary>) -> HtmlSample {
        HtmlSample { name: name.to_owned(), extract_source: name.to_owned(), species }
    }

    fn species(tax_id: u32, depth: usize, loci: usize, reference: &str) -> SpeciesSummary {
        SpeciesSummary { tax_id, depth, loci, references: vec![reference.to_owned()] }
    }

    #[test]
    fn builds_and_aggregates_the_induced_family_tree() {
        let samples = [
            sample("one", vec![species(30, 3, 2, "ref-a")]),
            sample("one", vec![species(30, 4, 1, "ref-b"), species(50, 2, 1, "ref-c")]),
        ];
        let tree = build_tree(&samples.iter().collect::<Vec<_>>(), &taxonomy());

        let alpha = &tree.children[&Some(10)];
        assert_eq!(alpha.depth, 7);
        assert_eq!(alpha.loci, 3);
        assert_eq!(alpha.samples.len(), 2);
        assert_eq!(alpha.extract_calls.len(), 2);
        assert_eq!(alpha.references.iter().cloned().collect::<Vec<_>>(), ["ref-a", "ref-b"]);
        assert!(alpha.children[&Some(11)].children[&Some(20)].children.contains_key(&Some(30)));
        assert_eq!(tree.children[&Some(40)].depth, 2);
    }

    #[test]
    fn groups_species_without_a_family() {
        let sample = sample("one", vec![species(60, 1, 1, "ref")]);
        let tree = build_tree(&[&sample], &taxonomy());

        let unclassified = &tree.children[&None];
        assert_eq!(unclassified.name, "Unclassified family");
        assert!(unclassified.children.contains_key(&Some(60)));
    }

    #[test]
    fn renders_expandable_escaped_multi_sample_html() {
        let samples = [
            sample("one & only", vec![species(30, 3, 2, "ref\"<a>")]),
            sample("two", vec![species(50, 2, 1, "ref-b")]),
        ];
        let html = render_report(&samples, &taxonomy(), Err("Krona <unavailable>"));

        assert!(html.contains("<details><summary class=\"row\" style=\"--depth:0\">"));
        assert!(html.contains("class=\"row leaf species-call\" style=\"--depth:3\""));
        assert!(html.contains("class=\"taxon-cell\""));
        assert!(html.contains("aria-hidden=\"true\">▸</span>"));
        assert!(html.contains("class=\"row leaf species-call\""));
        assert!(html.contains("class=\"extract-button\""));
        assert!(html.contains("class=\"references-button\""));
        assert!(html.contains("class=\"report-layout\""));
        assert!(html.contains("class=\"toolbar report-toolbar\""));
        assert!(html.contains("class=\"card taxa-card\">\n      <div class=\"toolbar report-toolbar\""));
        assert!(html.contains("class=\"card taxa-card\""));
        assert!(html.contains("class=\"card krona-card\""));
        assert!(html.contains("@media (min-aspect-ratio: 3 / 2)"));
        assert!(!html.contains("@media (max-width:"));
        assert!(html.contains("grid-template-columns: minmax(24ch, 2.5fr)"));
        assert!(!html.contains("grid-template-columns: subgrid"));
        assert!(html.contains("position: sticky"));
        assert!(!html.contains("--toolbar-height"));
        assert!(html.contains("id=\"extract-dialog\""));
        assert!(html.contains("id=\"references-dialog\""));
        assert!(!html.contains("<details open"));
        assert!(!html.contains("Georgia"));
        assert!(!html.contains("Alignment-supported species"));
        assert!(html.contains("one &amp; only"));
        assert!(html.contains("Alpha &lt;species&gt;"));
        assert!(html.contains("ref&quot;&lt;a&gt;"));
        assert!(html.contains("Krona &lt;unavailable&gt;"));
        assert!(html.contains("id=\"panel-2\" hidden"));
        assert!(html.contains("type: 'bamtax-dataset'"));
        assert!(html.contains("message.type === 'bamtax-dataset-selected'"));
        assert!(html.contains("selectDataset(message.index, false)"));
        assert!(html.contains("--species: #fff0d8"));
        assert!(!html.contains("<!-- BAMTAX:"));
    }

    #[test]
    fn renders_sorted_reference_count_buttons() {
        let samples = [sample(
            "one",
            vec![SpeciesSummary {
                tax_id: 30,
                depth: 1,
                loci: 1,
                references: vec!["ref-b".to_owned(), "ref-<a>".to_owned()],
            }],
        )];
        let html = render_report(&samples, &taxonomy(), Err("no Krona"));

        assert!(html.contains("data-references=\"ref-&lt;a&gt;\nref-b\""));
        assert!(html.contains("aria-label=\"Show 2 references for Alpha &lt;species&gt;\">2</button>"));
    }

    #[test]
    fn creates_streaming_awk_commands_for_parent_nodes() {
        let samples =
            [sample("one", vec![species(30, 3, 2, "ref-a")]), sample("two", vec![species(30, 4, 1, "ref-b")])];
        let tree = build_tree(&samples.iter().collect::<Vec<_>>(), &taxonomy());
        let command = extract_command(&tree.children[&Some(10)].extract_calls);

        assert_eq!(
            command,
            "samtools view 'one' | awk -F '\\t' 'BEGIN { keep[\"ref-a\"] = 1 } $3 in keep'; \
             samtools view 'two' | awk -F '\\t' 'BEGIN { keep[\"ref-b\"] = 1 } $3 in keep'"
        );
    }

    #[test]
    fn quotes_paths_and_escapes_awk_strings() {
        let command = extract_command(&[ExtractCall {
            source: "one's reads".to_owned(),
            references: "ref\"1;ref\\2".to_owned(),
        }]);

        assert_eq!(
            command,
            "samtools view 'one'\"'\"'s reads' | awk -F '\\t' 'BEGIN { keep[\"ref\\\"1\"] = 1; keep[\"ref\\\\2\"] = 1 } $3 in keep'"
        );
    }

    #[test]
    fn isolates_embedded_krona_document_in_srcdoc() {
        let html = render_report(&[], &taxonomy(), Ok("<!doctype html><script>window.x = \"yes\";</script>"));

        assert!(html.contains(
            "<iframe title=\"Krona taxonomy chart\" sandbox=\"allow-scripts allow-popups allow-downloads\" srcdoc=\"&lt;!doctype html&gt;"
        ));
        assert!(html.contains("&lt;script&gt;window.x = &quot;yes&quot;;&lt;/script&gt;"));
        assert!(!html.contains("src=\"style.css\""));
    }

    #[test]
    fn injects_dataset_bridge_into_krona() {
        let html = inject_krona_bridge("<html><body>chart</body></html>".to_owned());

        assert!(html.contains("message.type === 'bamtax-dataset'"));
        assert!(html.contains("window.selectDataset(message.index)"));
        assert!(html.contains("type: 'bamtax-dataset-selected'"));
        assert!(html.contains("bamtaxSelectDataset(index)"));
        assert!(html.contains("hideKronaDatasetControls"));
        assert!(html.contains("controls.hidden = true"));
        assert!(html.find("bamtax-dataset").unwrap() < html.find("</body>").unwrap());
    }

    #[test]
    fn reports_a_missing_krona_executable_without_panicking() {
        let sample = sample("one", vec![species(30, 1, 1, "ref")]);
        let error = generate_krona(Path::new("report.html"), &[sample], "bamtax-command-that-does-not-exist", None)
            .unwrap_err();

        assert!(error.contains("could not run bamtax-command-that-does-not-exist"));
    }
}
