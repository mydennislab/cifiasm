# CiFi Assembly Workflow (active development)

A Snakemake pipeline for generating chromosome-scale, phased de novo assemblies using HiFi long reads and CiFi (Long-reads Chromatin Conformation Capture) data, with support for manual curation via Juicebox Assembly Tools (JBAT).

## Overview

This pipeline takes HiFi reads (BAM or pre-built FASTA) and CiFi reads (BAM, FASTQ, or FASTA) as input and produces haplotype-resolved, scaffold-level assemblies with QC metrics, contact maps, and editable `.hic` files for manual curation in Juicebox.


## Workflow Diagram

![CiFi Assembly Workflow](assets/flowchart.png)


## Prerequisites

### Dependencies

Managed via conda (`environment.yaml`). Key tools: snakemake, samtools, hifiasm, gfatools, yahs, k8, nextflow, cifi (PyPI), biopython, numpy, pairix.

### 3D-DNA (required for JBAT)

Included as a git submodule. After cloning this repo:

```bash
git submodule update --init --recursive
```

### JuicerTools JAR (required for JBAT)

Download the JuicerTools JAR file:

```bash
wget https://github.com/aidenlab/JuicerTools/releases/download/v3.0.0/juicer_tools.jar -O juicer_tools.3.0.0.jar
```


## Configuration

Copy `config.example.yaml` to `config.yaml` and edit:

```bash
cp config.example.yaml config.yaml
```

### Samples
```yaml
samples:
  my_sample:
    hifi_bam: /path/to/hifi_reads.bam      # PacBio HiFi reads (BAM)
    # OR, instead of hifi_bam:
    # hifi_fasta: /path/to/hifi.fa.gz      # pre-built HiFi FASTA (.fa / .fa.gz);
    #                                      # skips BAM->FASTA conversion and is
    #                                      # incompatible with hifi.downsample.
    cifi: /path/to/cifi_reads.bam          # CiFi reads (BAM, FASTQ, or FASTA; .gz OK)
    enzyme: HindIII                        # Restriction enzyme (HindIII, DpnII, NlaIII, ...)
    # Optional: custom restriction site instead of a named enzyme.
    # site: "GANTC"                         # IUPAC recognition sequence
    # cut_pos: 1                            # 0-based cut position within the site
```
Exactly one of `hifi_bam:` or `hifi_fasta:` must be set per sample.

### Output directory (optional, default `results`)
```yaml
output_dir: results
```
All pipeline outputs land under this directory in their existing subdir
structure. Set to a different path to write elsewhere.

### Downsampling (optional)

Three independent ways to downsample reads. Each scenario shows up as a
single `{label}` value in output paths.

**CiFi-only sweep (existing, unchanged):**
```yaml
dilution:
  enabled: false                       # use 100% CiFi reads
  percentages: [20, 40, 60, 80, 100]   # or sweep across coverage levels
```
Labels are the numeric percentages: `20`, `40`, `100`, ...

**HiFi-only sweep (new):**
```yaml
hifi:
  downsample:
    enabled: true
    mode: depth                  # "depth" | "fraction"
    genome_size: 2500000000      # bp (haploid). Required for depth mode.
    depths:      [5, 10, 30]     # target X coverage  (depth mode)
    percentages: [20, 50]        # % of total reads   (fraction mode)
```
Depth-mode labels look like `h5X`, `h10X`, `h30X`. Fraction-mode labels look
like `h20pct`, `h50pct`. The pipeline computes the per-sample fraction from
each sample's actual HiFi base count via `samtools stats` (cached on disk).

**HiFi + CiFi together (zip):** when both `hifi.downsample.enabled` and
`dilution.enabled` are `true`, the lists are paired index-wise (must be the
same length). Labels combine both: `h5X_c20`, `h10X_c50`, ...

**Pre-downsampled CiFi (external):** when you already have one or more
downsampled CiFi inputs and don't want the pipeline to re-sample, declare
them per-sample under `cifi_external`:
```yaml
samples:
  my_sample:
    hifi_bam: /data/hifi.bam
    cifi: /data/cifi.bam          # still required, used only for cifi_qc
    cifi_external:
      c10:  /data/cifi.10pct.bam
      c25:  /data/cifi.25pct.fa.gz
      c100: /data/cifi.bam
    enzyme: HindIII
```
Each entry may be BAM, FASTQ, or FASTA (gzipped or not). BAMs are symlinked
in; FASTQ/FASTA are wrapped into an unmapped BAM via `samtools import`. No
re-sampling happens. Each key becomes a `{label}` in the output paths (must
match `[A-Za-z0-9._-]+`). Cannot be combined with `dilution.enabled` or
`hifi.downsample.enabled`. When multiple samples are defined, every sample
must declare the same label set.

### Tool Paths
```yaml
tools:
  singularity_cache: /path/to/cache         # For Nextflow containers
  threed_dna: "./3d-dna"                    # Path to 3D-DNA repository
  juicer_tools_jar: "./juicer_tools.3.0.0.jar"  # Path to JuicerTools JAR
```

### CiFi toolkit options (optional)
```yaml
cifi:
  qc:
    num_reads: 0         # reads to sample for QC (0 = all)
    min_sites: 1         # min enzyme sites to count a read as usable
  digest:
    min_fragments: 3     # min fragments per read to keep (-m)
    min_frag_len: 20     # min fragment length bp (-l)
    strip_overhang: true # default: strip enzyme overhang from R2
    gzip: false          # gzip-compress R1/R2 output
    fast: false          # streaming stats (lower memory)
```
See the [cifi toolkit](https://pypi.org/project/cifi/) docs for details. To use a custom restriction site instead of a named enzyme, set `site` + `cut_pos` under the sample entry.

### SLURM (optional)
```yaml
slurm:
  partition: "low"       # SLURM partition name
  account: "publicgrp"   # SLURM account name
```

## Running

```bash
conda activate cifiasm

# Dry run (validate DAG)
snakemake --dry-run

# Run locally
snakemake --cores 32

# Run on SLURM cluster
snakemake --cores 32 --slurm

# Run specific target
snakemake results/stats/my_sample/100/summary.tsv         # contig stats only
snakemake results/stats/my_sample/100/yahs_summary.tsv    # scaffold stats only
```


## Output Files

All pipeline outputs land under a single directory (default: `results/`).
Override with `output_dir: /some/other/path` in `config.yaml`. The label
`{label}` encodes the downsampling scenario (e.g. `100`, `h10X`, `h10X_c20`).

```
results/
├── qc_cifi/{sample}/
│   └── qc.pdf                            # CiFi QC report
├── hifi/
│   ├── {sample}.{label}.hifi.bam         # (down)sampled HiFi BAM (only when hifi_bam is set)
│   └── {sample}.{label}.hifi.fa          # HiFi FASTA fed to hifiasm
├── cifi/
│   └── {sample}.{label}.{bam,bam.bai}    # per-label CiFi BAM
├── cifi2pe/
│   └── {sample}.{label}_R{1,2}.fastq     # Hi-C-like paired reads (from cifi digest)
├── asm/{sample}/{label}/
│   ├── *.hic.hap{1,2}.p_ctg.gfa          # hifiasm contigs (GFA)
│   └── *.hap{1,2}.fa                     # Contigs (FASTA)
├── porec/{sample}/{label}/hap{1,2}/
│   ├── bed/*.bed                         # Pore-C contacts
│   ├── pairs/*.pairs.gz                  # Contact pairs
│   ├── pairs/*.mcool                     # Multi-resolution contact matrix
│   ├── hi-c/*.hic                        # Hi-C contact map
│   └── wf-pore-c-report.html             # wf-pore-c HTML report
├── yahs/{sample}/{label}/
│   ├── *_scaffolds_final.fa              # Final scaffolds
│   └── *_scaffolds_final.agp             # Scaffold AGP
├── qc_porec/{sample}/{label}/hap{1,2}/
│   ├── bed/*.bed
│   ├── pairs/{*.pairs.gz,*.mcool}
│   ├── hi-c/*.hic                        # QC contact map
│   ├── bams/*.cs.bam                     # Aligned CiFi reads (coordinate-sorted)
│   ├── paired_end/*.ns.bam               # Name-sorted paired BAM (feeds JBAT)
│   └── wf-pore-c-report.html
├── jbat/{sample}/{label}/hap{1,2}/
│   ├── *.hic                             # Editable Hi-C map for Juicebox
│   ├── *.assembly                        # Assembly file for JBAT
│   └── merged_nodups.txt                 # Contact data (3D-DNA format)
└── stats/{sample}/{label}/
    ├── summary.tsv                       # Contig assembly stats
    └── yahs_summary.tsv                  # Scaffold stats
```

## Manual Curation with JBAT

After the pipeline completes, use Juicebox Assembly Tools for manual curation:

1. Open `results/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic` in [Juicebox](https://github.com/aidenlab/Juicebox)
2. Load the `.assembly` file via **Assembly > Import Map Assembly**
3. Review and correct scaffold joins/orientations
4. Export the corrected assembly as `{sample}.{label}.hap{hap}.review.assembly`
5. Run the post-review rule to generate the final FASTA:
   ```bash
   snakemake results/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.FINAL.fa
   ```


## Workflow Rules

| Rule | Description |
|------|-------------|
| `cifi_qc` | QC report on raw CiFi input (via `cifi qc`) |
| `cifi_fastq_to_bam` | Wrap canonical CiFi FASTQ/FASTA into an unmapped BAM (when `cifi:` is not a BAM) |
| `downsample_hifi_bam` | (Down)sample HiFi BAM to a label-specific copy (skipped when `hifi_fasta:` is set) |
| `hifi_fasta` | Produce per-label HiFi FASTA: `samtools fasta` from BAM, or symlink/gunzip from `hifi_fasta:` |
| `downsample_cifi_bam` | Produce per-label CiFi BAM — BAM input sampled or symlinked, FASTQ/FASTA imported via `samtools import` |
| `cifi_fastq_from_downsampled_bam` | Extract FASTQ from per-label CiFi BAM |
| `cifi2pe_split` | Digest CiFi reads into Hi-C-like PE reads (via `cifi digest`) |
| `hifiasm_dual_scaf` | Assemble with hifiasm --dual-scaf |
| `gfa2fa` | Convert GFA to FASTA |
| `caln50` | Calculate N50 and other stats |
| `summarize_assembly` | Compile contig assembly statistics |
| `porec_nextflow` | Run wf-pore-c for contact mapping |
| `index_fa` | Index FASTA with samtools |
| `yahs_scaffold` | Scaffold with YAHS |
| `yahs_caln50` | Calculate scaffold statistics |
| `summarize_yahs` | Compile scaffold statistics |
| `yahs_index_scaffolds_fa` | Index scaffold FASTA |
| `qc_porec_nextflow` | QC by mapping CiFi to scaffolds |
| `cifi_filter_bam` | Filter paired-end BAM for JBAT |
| `generate_genome_file` | Generate chromosome sizes file |
| `bam2pairs` | Convert BAM to pairs format |
| `pairs_to_mnd` | Convert pairs to merged_nodups.txt |
| `generate_assembly_file` | Generate .assembly from scaffold FASTA |
| `jbat_hic` | Generate editable .hic for Juicebox |
| `jbat_post_review` | Convert curated .assembly back to FASTA |


## Scripts & external tools

This workflow uses the following external scripts and tools:

- [`cifi`](https://pypi.org/project/cifi/) (PyPI) - CiFi QC (`cifi qc`) and in-silico
  restriction digestion to Hi-C-like paired-end reads (`cifi digest`). Options are
  configurable under the `cifi:` key in `config.yaml`.

- `scripts/calN50.js` - Calculates N50 and assembly statistics (requires k8)
  Source: [calN50](https://github.com/lh3/calN50) by Heng Li.

- [3D-DNA](https://github.com/aidenlab/3d-dna) - Assembly visualization and post-review tools
  by Aiden Lab.

- [JuicerTools](https://github.com/aidenlab/JuicerTools) - Hi-C file generation
  by Aiden Lab.
