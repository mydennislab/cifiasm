# CiFi Assembly Workflow

A Snakemake pipeline for generating chromosome-scale, phased de novo assemblies using HiFi long reads and CiFi (Long-reads Chromatin Conformation Capture) data, with support for manual curation via Juicebox Assembly Tools (JBAT).

## Overview

This pipeline takes HiFi reads and CiFi reads (each as BAM, FASTQ, or FASTA, one or more files per sample) and produces haplotype-resolved, scaffold-level assemblies with QC metrics, contact maps, and editable `.hic` files for manual curation in Juicebox.


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
    hifi: /path/to/hifi_reads.bam          # one path, or a list of files
    cifi: /path/to/cifi_reads.bam          # one path, or a list of files
    enzyme: HindIII                        # Restriction enzyme (HindIII, DpnII, NlaIII, ...)
    # Optional: custom restriction site instead of a named enzyme.
    # site: "GANTC"                         # IUPAC recognition sequence
    # cut_pos: 1                            # 0-based cut position within the site
```

`hifi:` and `cifi:` each take a single file or a list of files (e.g. several SMRT
cells). Both accept BAM, FASTQ, and FASTA (`.gz` is fine), detected by extension.
HiFi files go straight to hifiasm; CiFi files are merged before digestion.

```yaml
  multi_cell_sample:
    hifi:
      - /path/to/hifi.cell1.bam
      - /path/to/hifi.cell2.fastq.gz
    cifi:
      - /path/to/cifi.cell1.bam
      - /path/to/cifi.cell2.fastq.gz
    enzyme: HindIII
```

The legacy keys `hifi_bam:`, `hifi_fasta:`, and `cifi_bam:` still work as aliases.

### Output directory (optional, default `results`)
```yaml
output_dir: results
```
All pipeline outputs land under this directory in their existing subdir
structure. Set to a different path to write elsewhere.

### CiFi downsampling (optional)

Only CiFi is downsampled; HiFi always goes to hifiasm at full depth. Each
scenario becomes one `{label}` in the output paths.

To assemble at several CiFi coverage levels, enable the dilution sweep:
```yaml
dilution:
  enabled: false                       # use 100% CiFi reads
  percentages: [20, 40, 60, 80, 100]   # or sweep across these percentages
```
The label is the percentage: `20`, `40`, `100`, ...

If you already have downsampled CiFi files and don't want the pipeline to
re-sample, list them per-sample under `cifi_external`:
```yaml
samples:
  my_sample:
    hifi: /data/hifi.bam
    cifi: /data/cifi.bam          # still required, used for cifi_qc
    cifi_external:
      c10:  /data/cifi.10pct.bam
      c25:  /data/cifi.25pct.fa.gz
      c100: /data/cifi.bam
    enzyme: HindIII
```
Each entry may be BAM, FASTQ, or FASTA (`.gz` OK). BAMs are linked in;
FASTQ/FASTA are imported to BAM. Each key becomes a `{label}` (must match
`[A-Za-z0-9._-]+`). This cannot be combined with `dilution.enabled`, and when
multiple samples are defined they must share the same set of labels.

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
Override with `output_dir: /some/other/path` in `config.yaml`. The `{label}`
encodes the CiFi scenario (e.g. `100`, or a dilution percentage like `20`).

```
results/
├── qc_cifi/{sample}/
│   └── qc.pdf                            # CiFi QC report
├── hifi/{sample}/
│   └── cell{n}.fastq                     # HiFi BAM cells converted to FASTQ (FASTQ/FASTA bypass this)
├── cifi/
│   ├── merged/{sample}.cifi.bam          # all CiFi cells merged
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
| `merge_cifi` | Merge all CiFi cells into one BAM (`samtools import` + `samtools cat`) |
| `cifi_qc` | QC report on the merged CiFi BAM (via `cifi qc`) |
| `hifi_bam_to_fastq` | Convert a HiFi BAM cell to FASTQ (FASTQ/FASTA cells skip this) |
| `downsample_cifi_bam` | Produce per-label CiFi BAM — sampled with `samtools view -s`, or symlinked at 100% |
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


## Citation

If you use this tool, please cite:

> Abuelanin M, Kaya G, Lake JA, Lambert C, Wu MV, Berendzen KM, Krasheninnikova K, Wood JMD, Solomon NG, Donaldson ZR, Bales KL, Howe K, Korlach J, Manoli D, Tollkuhn J, Dennis MY. Single-library chromosome-scale diploid assemblies of vole genomes resolve a species-specific duplication implicated in pair bonding. bioRxiv 2026.03.13.711624. https://doi.org/10.64898/2026.03.13.711624


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
