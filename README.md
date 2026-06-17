# meta_marker_count

`meta_marker_count` counts marker-derived reads from paired-end clean reads. It uses BBDuk to extract candidate 16S/18S/ITS reads, minimap2 to align them to marker references, and a taxonomy table to produce marker RPM tables.

The default marker set is `16S,ITS`. Add 18S explicitly with `--markers 16S,18S,ITS`.

## Installation

Clone the repository and install the command wrappers:

```bash
git clone <repo-url>
cd meta_marker_count

mamba env create -f environment.yml
mamba activate meta-marker-count

make install PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

Installed commands:

- `meta_marker_count`: main pipeline.
- `meta_marker_build_refs`: converts SILVA/UNITE FASTA files into the default reference layout and minimap2 indexes.

Update an existing installation:

```bash
cd meta_marker_count
git pull --ff-only
make install PREFIX="$HOME/.local"
```

## Default Reference Directory

The software has built-in reference paths. By default, reference files are stored in:

```text
~/.local/share/meta_marker_count/ref
```

Override this directory with either:

```bash
export META_MARKER_COUNT_REF_DIR=/path/to/meta_marker_count_ref
```

or at runtime:

```bash
meta_marker_count --ref-dir /path/to/meta_marker_count_ref ...
```

After references are built in the default directory, normal pipeline runs do not need `--ref-16s`, `--index-16s`, `--ref-its`, `--taxonomy`, or related options.

Default files expected by `meta_marker_count`:

```text
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/UNITE_public_19.02.2025.shortid.fasta
~/.local/share/meta_marker_count/ref/UNITE_public_19.02.2025.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/ref_taxonomy.tsv
```

`meta_marker_count --help` also prints these default paths.

## Download Raw Reference Files

Download SILVA SSU directly:

```bash
curl -L \
  -o SILVA_138.2_SSURef_tax_silva.fasta.gz \
  "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz"
```

Download UNITE ITS from the repository page:

```text
https://unite.ut.ee/repository.php
```

Select the public FASTA gzip release and save it as:

```text
UNITE_public_19.02.2025.fasta.gz
```

## Build And Configure References

Run this once after downloading the raw reference files:

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --force
```

By default this writes to `~/.local/share/meta_marker_count/ref` and builds `.mmi` indexes next to the generated FASTA files. The main pipeline will automatically use these paths.

If minimap2 indexes should be created later, use:

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --skip-index \
  --force
```

Then build indexes manually:

```bash
REF_DIR="${META_MARKER_COUNT_REF_DIR:-$HOME/.local/share/meta_marker_count/ref}"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta"

minimap2 -d "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta.mmi" \
  "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta"
```

Generated reference files:

```text
UNITE_public_19.02.2025.shortid.fasta
UNITE_public_19.02.2025.taxonomy.tsv
SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta
SILVA_138.2_SSURef_tax_silva.dna.arc_bac.taxonomy.tsv
SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta
SILVA_138.2_SSURef_tax_silva.dna.euk.taxonomy.tsv
ref_taxonomy.tsv
```

Combined taxonomy format:

```tsv
ref_id	marker	domain	phylum	class	order	family	genus	species
AB000393.1.1510	16S	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Vibrionaceae	Vibrio	Vibrio_halioticoli
```

Reference parsing rules:

- UNITE ITS headers such as `UDB016649|k__Fungi;...|SH1281904.10FU` use `UDB016649` as `ref_id` and marker `ITS`.
- SILVA `Bacteria` and `Archaea` records are written to the arc/bac FASTA and marked as `16S`.
- SILVA `Eukaryota` records are written to the euk FASTA and marked as `18S`.
- SILVA eukaryotic lineages can have more than seven levels; the last two levels are retained as `genus` and `species`.
- SILVA records with an exact lineage item `Mitochondria` or `Chloroplast` are skipped by default. Use `--keep-organelles` to keep them.

## Sample Input

Multi-sample input is a TSV file with a header. Columns are matched by name, not by order.

Required columns:

- `sample_id`
- `r1_path`
- `r2_path`

Optional columns:

- `year`
- `month`
- `depth`
- any other sample metadata column

Example:

```tsv
sample_id	year	month	depth	site_type	r1_path	r2_path
202311_MF1_00-10	2023	11	00-10	MF	/data/202311_MF1_00-10_R1.fq.gz	/data/202311_MF1_00-10_R2.fq.gz
202311_MF1_10-20	2023	11	10-20	MF	/data/202311_MF1_10-20_R1.fq.gz	/data/202311_MF1_10-20_R2.fq.gz
```

## Run The Pipeline

Default 16S + ITS run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

16S + 18S + ITS run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,18S,ITS \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

Single-sample run:

```bash
meta_marker_count \
  --sample-id 202311_MF1_00-10 \
  --r1 /data/202311_MF1_00-10_R1.fq.gz \
  --r2 /data/202311_MF1_00-10_R2.fq.gz \
  --outdir marker_count_out
```

Dry run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --dry-run
```

## Key Options

- `--steps all`: run all steps. Subsets are allowed: `reads_count,extract,align,abundance`.
- `--markers 16S,ITS`: marker list. Allowed markers are `16S`, `18S`, `ITS`.
- `--rank genus`: output rank. Allowed ranks are `domain,phylum,class,order,family,genus,species,all`.
- `--ref-dir DIR`: use a non-default reference directory.
- `--jobs INT`: samples processed in parallel.
- `--threads-per-sample INT`: BBDuk/minimap2 threads per sample.
- `--bbmap-mem 64G`: Java heap for BBDuk.
- `--force`: overwrite and rerun existing outputs.
- `--no-resume`: ignore checkpoints.
- `--clean-tmp`: remove temporary directories after a successful run.

## Output

Default output directory:

```text
marker_count_out/
  00_manifest/samples.tsv
  00_logs/marker_count.YYYYmmdd_HHMMSS.log
  00_logs/run_config.tsv
  01_reads_count/clean_reads.tsv
  02_marker_reads/
  03_align/
  04_abundance/
    marker_rpm.<rank>.long.tsv
    marker_rpm.<rank>.matrix.tsv
    marker_rpm.<rank>.domain_total.tsv
    marker_rpm.<rank>.assignment_stats.tsv
  99_checkpoints/
```

Important result tables:

- `marker_rpm.<rank>.long.tsv`: one row per sample/marker/domain/taxon.
- `marker_rpm.<rank>.matrix.tsv`: sample by feature matrix.
- `marker_rpm.<rank>.domain_total.tsv`: total marker reads/RPM per domain and marker.
- `marker_rpm.<rank>.assignment_stats.tsv`: alignment filtering and ambiguity statistics.

## Checks

Run syntax checks:

```bash
make check
```

Show command defaults:

```bash
meta_marker_count --print-defaults
meta_marker_count --help
```
