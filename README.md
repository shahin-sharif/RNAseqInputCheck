# RNAseqInputCheck

[![Tests](https://github.com/shahin-sharif/RNAseqInputCheck/actions/workflows/tests.yml/badge.svg)](https://github.com/shahin-sharif/RNAseqInputCheck/actions/workflows/tests.yml)

**Check RNA-seq inputs before starting a long analysis.** RNAseqInputCheck validates sample sheets, GTF annotations, Salmon quantifications, coordinate-sorted BAM files, and experimental design matrices. It produces an offline HTML report, machine-readable JSON, and tab-separated audit files.

This tool grew out of practical problems in differential expression, differential junction usage and isoform-switch workflows: incompatible transcript annotation sets, chromosome naming differences, and accidentally replacing stable gene IDs with shared gene symbols.

It reads your inputs without modifying them. It does not run differential expression, repair annotation, install software, or upload data. It can be used before [Transcriptome_Analysis](https://github.com/shahin-sharif/Transcriptome_Analysis) or another RNA-seq analysis workflow.

## What it checks

| Input | Checks |
|---|---|
| Sample sheet | Required columns; unique sample IDs; existing files; reused input paths; group sizes |
| GTF | Nine-column structure and valid coordinates; exon gene/transcript IDs; transcript mapping consistency; stable gene IDs on multiple sequences; shared gene symbols |
| Salmon | Required `quant.sf` fields; finite nonnegative estimates; duplicate transcript IDs; exact transcript-ID coverage by GTF exons; common ID sets, order and lengths across samples; TPM totals |
| Salmon metadata | Expected library type, when specified; index sequence/name hash consistency when present; mapping percentage recorded |
| BAM | `samtools quickcheck`; declared coordinate sort order; sequence dictionary consistency; index presence and an indexed query |
| BAM and GTF | Exact sequence-name overlap; missing human-style primary chromosome names; annotation coordinates beyond reference lengths; absent reference sequences |
| BAM and FASTA index | Optional exact sequence-name and length comparison against a `.fai` file |
| Design | Restricted additive/interaction formula; categorical reference levels; numeric covariates; matrix rank and residual degrees of freedom |
| R | Optional package loading/version check and presence of an IsoformSwitchAnalyzeR import argument |

**A passing report means the requested structural checks passed. It does not mean a complete downstream analysis has been validated.** See [limitations](#limitations) before interpreting a report.

## 1. Download the tool

Run these commands in your server's Bash terminal, not inside R:

```bash
# Create a directory for software and enter it.
mkdir -p "$HOME/software"
cd "$HOME/software"

# Download the repository.
git clone https://github.com/shahin-sharif/RNAseqInputCheck.git
cd RNAseqInputCheck

# Confirm Python is available. Version 3.9 or newer is required.
python3 --version
python3 RNAseqInputCheck.py --version
```

No Python libraries need to be installed. Running `python3 RNAseqInputCheck.py` tells Python to execute this file; you do not need to open an interactive Python session.

If checking BAM files, install **samtools** or use an existing installation:

```bash
command -v samtools
samtools --version
```

An optional separate Conda environment can provide both programs:

```bash
conda create --name rnaseq_input_check \
  --override-channels -c conda-forge -c bioconda \
  python=3.12 samtools -y
conda activate rnaseq_input_check
```

R is **not required** unless you request the optional R package check. For that check, activate the same environment that your analysis uses. This ensures `Rscript` probes the intended library.

Optional Python package installation, instead of running the repository script:

```bash
python3 -m pip install .
rnaseq-input-check --help
```

## 2. Run the included small example

```bash
python3 RNAseqInputCheck.py \
  --config examples/tiny/config.json \
  --out results/first_check
```

The example contains four synthetic Salmon quantifications and a two-transcript GTF. It is a software demonstration, not a biological benchmark. BAM and R checks are skipped because these inputs were not requested. Expected result: **28 PASS, 2 INFO, 2 SKIP, no errors or warnings**.

The output folder must be new. To repeat the example, use a different name, such as `results/second_check`. Existing reports are never overwritten.

Open `results/first_check/report.html` in a browser. On a Mac that is logged into a remote server, copy the report from a separate **local Mac terminal**:

```bash
# Replace USER and SERVER with your own SSH login and hostname.
scp USER@SERVER:~/software/RNAseqInputCheck/results/first_check/report.html \
  "$HOME/Downloads/RNAseqInputCheck_report.html"
open "$HOME/Downloads/RNAseqInputCheck_report.html"
```

The HTML file is self-contained and requires no internet connection. Keep the other output files for the detailed audits.

## 3. Create your sample sheet

A sample sheet is a plain-text table separated by **tabs**, with one sample per row. Required columns are `sample_id`, `condition`, and at least one of `bam` or `quant`.

- `sample_id`: unique label for the sample.
- `condition`: experimental group.
- `bam`: path to the coordinate-sorted BAM. Its index should be beside it.
- `quant`: path to a Salmon output directory or directly to its `quant.sf` file.
- Additional columns, such as `batch` or `time`, can be used in the design.

Omit the whole `bam` column for a Salmon-only check, or the whole `quant` column for a BAM-only check. An included input column must be populated for every sample.

Create a working folder and a sheet using Python so that the tabs are written correctly:

```bash
mkdir -p "$HOME/RNAseq_preflight"

python3 - <<'PY'
import csv
from pathlib import Path

# Replace these example input directories with your actual directories.
bam_dir = Path.home() / 'my_experiment' / 'bam' / 'sorted_for_analysis'
quant_dir = Path.home() / 'my_experiment' / 'SalmonQuant'
out = Path.home() / 'RNAseq_preflight' / 'samples.tsv'

# Edit these sample labels and group assignments to match your experiment.
samples = [('WT_1', 'WT'), ('WT_2', 'WT'), ('WT_3', 'WT'),
           ('KO_1', 'KO'), ('KO_2', 'KO'), ('KO_3', 'KO')]

with out.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['sample_id', 'condition', 'bam', 'quant'])
    for sample, condition in samples:
        writer.writerow([sample, condition,
                         str(bam_dir / (sample + '.sorted.bam')),
                         str(quant_dir / sample)])
print(out)
PY

cat "$HOME/RNAseq_preflight/samples.tsv"
```

Absolute paths are easiest to understand. Relative BAM/Salmon paths are resolved relative to **the sample sheet's directory**, not the terminal's working directory. `~` is expanded to your home directory. Shell variables such as `$HOME` are not expanded inside JSON or TSV values; use `~` or an absolute path there.

## 4. Create the configuration

Create and edit a plain JSON file:

```bash
nano "$HOME/RNAseq_preflight/config.json"
```

Paste the following, replacing the annotation path:

```json
{
  "samples": "samples.tsv",
  "gtf": "~/reference/gencode.v48.annotation.gtf",
  "design": "~ condition",
  "reference_levels": {"condition": "WT"},
  "expected_library_type": "ISR"
}
```

Save in nano with **Ctrl+O**, press **Enter**, and exit with **Ctrl+X**.

`ISR` is Salmon's inward-facing, paired-end, reverse-stranded library type. Use it only when appropriate for your libraries. Remove `expected_library_type` if you do not know it; the tool will still report available metadata but will not establish an expected-type match. It does not infer strandedness from alignments.

Use the annotation corresponding to the transcript sequences actually used for your Salmon index. Even files from the same GENCODE release can cover different sequence sets. The checker compares exact transcript IDs, including version suffixes. GTF transcripts absent from Salmon are allowed; Salmon transcripts absent from GTF exons are errors.

Configuration file paths are relative to **the configuration file's directory**. This is why `"samples": "samples.tsv"` works when the two files are together.

Validate the JSON syntax and run:

```bash
python3 -m json.tool "$HOME/RNAseq_preflight/config.json"

python3 RNAseqInputCheck.py \
  --config "$HOME/RNAseq_preflight/config.json" \
  --out "$HOME/RNAseq_preflight/check01"
```

For a simple check without a configuration file:

```bash
python3 RNAseqInputCheck.py \
  --samples "$HOME/RNAseq_preflight/samples.tsv" \
  --gtf "$HOME/reference/genes.gtf.gz" \
  --out "$HOME/RNAseq_preflight/check02"
```

Do not combine `--config` with `--samples` or `--gtf`. The configuration is specific to RNAseqInputCheck; do not pass another pipeline's JSON file directly.

## 5. Interpret the report

| Status | Meaning |
|---|---|
| ERROR | Resolve before proceeding with the intended analysis. |
| WARN | Review the finding and document whether it is acceptable for your analysis. |
| PASS | The named check passed, within its stated scope. |
| SKIP | This check was not performed. |
| INFO | Context or a limitation, not a successful test. |

Exit codes: `0` means no errors; `1` means validation errors (or warnings with `--fail-on-warning`); `2` means a configuration, command-line or report-writing problem. A zero exit code can include warnings and skipped checks.

```bash
# Optional strict mode for an automated workflow:
python3 RNAseqInputCheck.py \
  --config "$HOME/RNAseq_preflight/config.json" \
  --out "$HOME/RNAseq_preflight/check03" \
  --fail-on-warning

# Run immediately after the checker to see its exit code.
echo $?
```

To save console output while preserving the checker's exit code in Bash:

```bash
set -o pipefail
python3 RNAseqInputCheck.py \
  --config "$HOME/RNAseq_preflight/config.json" \
  --out "$HOME/RNAseq_preflight/check04" \
  2>&1 | tee "$HOME/RNAseq_preflight/check04.log"
```

### Output files

- `report.html`: offline report with findings and metrics.
- `report.json`: the same findings in machine-readable form.
- `checks.tsv`: one row per check.
- `input_manifest.tsv`: inspected paths, sizes and modification times; optional SHA-256 hashes.
- `resolved_config.json`: configuration with top-level input paths resolved.
- `design_matrix.tsv`: design coefficients for each sample when construction succeeds.
- Additional JSON audits when relevant: missing transcript IDs, shared gene symbols, stable gene IDs on multiple sequences, or GTF sequences absent from BAM references.

Reports contain sample labels and local paths. Review them before posting publicly. Inputs and credentials are never uploaded by the tool.

## Advanced configuration

All supported keys are listed below. Unknown keys produce an error to catch spelling mistakes.

| Key | Default / purpose |
|---|---|
| `samples`, `gtf` | Required paths. GTF can be `.gz`. |
| `fai` | Optional genome FASTA index for BAM sequence-name and length comparison. |
| `design` | `~ condition`. Includes an intercept. |
| `reference_levels` | Object mapping categorical variables to baseline levels. Otherwise lexical order determines the baseline. |
| `numeric_covariates` | List of design columns to treat as numbers; all others are categorical. |
| `expected_library_type` | Optional explicit Salmon type, e.g. `ISR`, `ISF`, `IU`, `SR`, `SF`, `U`. Not `A`. |
| `strip_pipe` | `false`. If true, take text before the first `\|` in Salmon names. Version suffixes are still retained; collisions are errors. |
| `samtools`, `rscript` | Executable names on PATH or absolute executable paths. Defaults: `samtools`, `Rscript`. |
| `r_packages` | Optional list of R packages to load and report. No installation is attempted. |
| `command_timeout` | 120 seconds per external command. Increase for slow network storage or R loading. |
| `hash_inputs` | `false`. If true, SHA-256 every recorded file, including BAMs; this can be very slow. |

Example with batch and a continuous variable:

```json
{
  "samples": "samples.tsv",
  "gtf": "genes.gtf.gz",
  "fai": "genome.fa.fai",
  "design": "~ batch + condition * time",
  "numeric_covariates": ["time"],
  "reference_levels": {"condition": "WT", "batch": "A"},
  "expected_library_type": "ISR",
  "r_packages": ["DESeq2", "edgeR", "Rsubread", "IsoformSwitchAnalyzeR"],
  "command_timeout": 300
}
```

Every design variable must have a sample-sheet column. Supported syntax is `+`, `:`, and `*` with bare variable names. The tool rejects unsupported R formula syntax instead of claiming to validate it. Interactions are expanded into products of dummy-coded or numeric columns. Rank is estimated by elimination after column scaling, with tolerance `1e-10`; this is not an assessment of numerical conditioning for a fitted model.

A single batch does not need a batch term. If batch and condition are identical assignments, they cannot be estimated separately. For two KO clones and one WT line with replicate cultures, a three-level condition can check separate clone comparisons; culture labels do not establish independent knockout events. Biological interpretation still needs experimental context.

### Common findings and their meaning

**BAM uses `1`, annotation uses `chr1`:** these names are not exact matches. A few matching alternate sequences must not conceal missing primary chromosomes. Obtain compatible reference files or prepare a separately verified alias conversion; do not blindly add `chr` to every name.

**GTF contains alternate/patch sequences absent from BAM:** those loci cannot be quantified from these BAMs. A complete transcript annotation may still be needed for Salmon. Decide the scope explicitly; the checker does not silently delete loci.

**Two stable gene IDs share a gene symbol:** a warning is appropriate. Keep the stable IDs separate. For curated GENCODE annotation imported into IsoformSwitchAnalyzeR, inspect `fixStringTieAnnotationProblem` and explicitly disable its gene-ID reassignment when needed. Loading the package alone does not validate an actual import.

**Missing BAM index or unsorted BAM:** prepare a new coordinate-sorted BAM and index outside this tool. Do not overwrite your original file. The checker does not decide whether deduplication is biologically appropriate.

**Salmon IDs are missing from the GTF:** check annotation release, primary versus full annotation sets, and the exact transcript FASTA used to build the index. Removing ID versions is not an automatic solution.

## Limitations

- Version 0.1.0 performs structural preflight, not full read-level QC, genome-sequence verification, read pairing checks, or empirical strandedness inference.
- `samtools quickcheck` checks the header and EOF, not all alignments. A declared sort order and successful indexed query do not prove complete sorting, index freshness or absence of internal corruption.
- Matching reference names and lengths do not prove identical sequences or assemblies. No FASTA sequence checksums are compared to BAM or Salmon index sequences.
- The primary-chromosome naming heuristic covers human-style `1–22`, `X`, `Y`, `M`/`MT`, optionally prefixed by `chr`. Other species still receive exact overlap/coordinate checks, but missing species-specific primary chromosomes are not classified automatically.
- Mapping percentages are reported, not judged by a universal cutoff. Metadata are trusted as metadata; they do not independently verify library preparation.
- GTF checks focus on exon identifiers, coordinates and mapping consistency. They do not establish transcript-model biology or prove exon structures match the indexed transcript sequences.
- A design rank check cannot detect mislabeled samples, hidden confounders, pseudoreplication or all unsuitable statistical models. Contrasts are not tested.
- R probes do not run a full pipeline or certify compatibility across every installed package.
- Do not change input files during a check. A manifest without hashes records file identity hints, not immutable provenance.

## Tests and development

```bash
python3 -m unittest discover -s tests -v
```

The tests cover ID mismatches, gene-symbol collisions, confounded designs, invalid numeric inputs, duplicate sample paths, report escaping and safe output handling. With samtools installed, an additional test creates real tiny BAM files, verifies a primary-chromosome alias failure despite an alternate-sequence match, and checks missing-index detection. GitHub Actions installs samtools and runs the suite on Python 3.9 and 3.12.

The bundled fixtures contain no real experimental data. Please report reproducible issues with small synthetic examples and software versions.

## Reference documentation

File-format and command semantics follow [GENCODE GTF documentation](https://www.gencodegenes.org/pages/data_format.html), [Salmon output documentation](https://salmon.readthedocs.io/en/stable/file_formats.html), and [samtools documentation](https://www.htslib.org/doc/samtools.html).

Maintained by Shahin Behrouz Sharif. MIT license. See [CHANGELOG.md](CHANGELOG.md) for releases.
