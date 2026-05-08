# scPT-Natural-Aging

This repository contains a Snakemake pipeline for processing scPaired-Tag data.

## Input

The pipeline uses two config files.

### 1. `config/samples.tsv`

This is a tab-separated file describing per-sample FASTQ inputs.

Required columns include:

- `R1`: path to read 1 FASTQ
- `R2`: path to read 2 FASTQ
- `sample`: sample name
- `library`: library type used by the Paired-Tag pipeline
- `species`: genome, for example `mm10`
- `seq_type`: sequencing type, for example `T7` or `nonT7`

If the same sample appears in multiple rows, the pipeline will merge all FASTQ files from that sample before running scPT.

Example:

```tsv
R1	R2	library	sample	tissue	depth	exp_prefix	exp_type	species	seq_type
/path/to/sample_L001_R1.fastq.gz	/path/to/sample_L001_R2.fastq.gz	RNA	NTY623-3	skin	shallow	NTY623-3	paired_tag_sc	mm10	nonT7
/path/to/sample_L002_R1.fastq.gz	/path/to/sample_L002_R2.fastq.gz	RNA	NTY623-3	skin	deep	NTY623-3	paired_tag_sc	mm10	nonT7
```

### 2. `config/experiments.csv`

This file lists experiment IDs in the `proj` column.

Example:

```csv
proj
20240620-scPT-Natural-Aging-liver-H3K27ac-H3K27me3
20240627-scPT-Natural-Aging-colon-H3K27ac-H3K27me3
```

For each experiment in `experiments.csv`, the pipeline expects a matching config file:

```bash
config/{exp}.sh
```

These shell config files define which sample BAMs belong to the same experiment and should be merged together.

## Workflow

The pipeline has two main steps.

### Step 1. Merge FASTQs and run scPT per sample

For each sample in `config/samples.tsv`:

1. Merge all `R1` FASTQs from the same sample.
2. Merge all `R2` FASTQs from the same sample.
3. Run the Paired-Tag pipeline on the merged FASTQs.
4. Generate one BAM file per sample.

### Step 2. Merge BAMs per experiment

For each experiment listed in `config/experiments.csv`:

1. Read `config/{exp}.sh`.
2. Parse the `dna_rna` pairs from that config file.
3. Collect all sample BAMs belonging to that experiment.
4. Run the merge script to generate merged experiment-level BAMs.

## Output

### Per-sample outputs

Merged FASTQs:

```bash
{pipeline_base}/{sample}/merge/fastq/{sample}_R1.fastq.gz
{pipeline_base}/{sample}/merge/fastq/{sample}_R2.fastq.gz
```

Per-sample BAM:

```bash
{pipeline_base}/{sample}/merge/{sample}/bam/{sample}_UsefulReads.bam
```

Set `pipeline_base` in `config/config.yaml`.

### Per-experiment outputs

Merged RNA BAM:

```bash
{merge_base}/{exp}/rna_bam/RNA_final.bam
```

Set `merge_base` in `config/config.yaml`.

## Run

Run the full workflow with:

```bash
snakemake --use-conda -p -c 10 -s workflow/Snakefile
```

## Main Files

- `workflow/Snakefile`: main workflow entry
- `workflow/rules/common.smk`: shared config and path helpers
- `workflow/rules/merge_fastq.smk`: merge FASTQs for the same sample
- `workflow/rules/run_paired_tag.smk`: run the Paired-Tag pipeline per sample
- `workflow/rules/merge_bam.smk`: merge BAMs for each experiment
- `scripts/generate_bam.sh`: experiment-level BAM merge script
