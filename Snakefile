import os
import re
import subprocess
import sys

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/.]+",
    label=r"[A-Za-z0-9._-]+"

_LABEL_RE = re.compile(r"^[A-Za-z0-9._-]+$")

def get_account_for_jobs(wildcards):
    return SLURM_ACCOUNT

def seeddotfrac_from_fraction(frac: float, seed: int = 100) -> str:
    """Render a fraction (0 < f <= 1) as samtools view -s 'SEED.FRAC'."""
    frac_digits = f"{frac:.10f}".split(".")[1].rstrip("0") or "0"
    return f"{seed}.{frac_digits}"

# Build sample list from config
SAMPLES = []
for sample_name, sample_data in config["samples"].items():
    # HiFi: accept either a BAM (`hifi_bam`) OR a pre-built FASTA (`hifi_fasta`).
    # Exactly one is required. FASTA skips the BAM -> samtools fasta conversion.
    hifi_bam_path   = sample_data.get("hifi_bam", "")
    hifi_fasta_path = sample_data.get("hifi_fasta", "")
    if bool(hifi_bam_path) == bool(hifi_fasta_path):
        raise ValueError(
            f"Sample '{sample_name}': exactly one of `hifi_bam` or `hifi_fasta` "
            f"must be set (got hifi_bam={hifi_bam_path!r}, hifi_fasta={hifi_fasta_path!r})."
        )
    # CiFi: accept either "cifi_bam" or "cifi" key (BAM, FASTQ, or FASTA)
    cifi_path = sample_data.get("cifi", sample_data.get("cifi_bam", ""))
    SAMPLES.append(dict(
        sample=sample_name,
        hifi_bam=hifi_bam_path,
        hifi_fasta=hifi_fasta_path,
        cifi_input=cifi_path,
        enzyme=sample_data.get("enzyme", ""),
        site=sample_data.get("site", ""),
        cut_pos=sample_data.get("cut_pos", ""),
        cifi_external=sample_data.get("cifi_external", {}) or {},
    ))

def samples_list():
    return [s["sample"] for s in SAMPLES]

def get_sample_data(sample_name):
    return next(s for s in SAMPLES if s["sample"] == sample_name)

def get_enzyme(wildcards):
    return get_sample_data(wildcards.sample)["enzyme"]

def get_enzyme_args(wildcards):
    """Return enzyme CLI args: either '-e NAME' or '--site SEQ --cut-pos N'."""
    data = get_sample_data(wildcards.sample)
    if data.get("site"):
        return f"--site {data['site']} --cut-pos {data['cut_pos']}"
    return f"-e {data['enzyme']}"

def get_digest_extra_flags():
    """Build extra CLI flags for cifi digest from config."""
    flags = []
    if not CIFI_DIGEST_OPTS.get("strip_overhang", True):
        flags.append("--revcomp-r2")
    if CIFI_DIGEST_OPTS.get("gzip", False):
        flags.append("--gzip")
    if CIFI_DIGEST_OPTS.get("fast", False):
        flags.append("--fast")
    return " ".join(flags)

_FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
_FASTA_EXTS = (".fasta", ".fa", ".fasta.gz", ".fa.gz")

def _path_is_bam(path):
    return str(path).lower().endswith(".bam")

def _path_is_fastx(path):
    """True for FASTQ/FASTA (gzipped or not) — anything `samtools import` accepts."""
    p = str(path).lower()
    return p.endswith(_FASTQ_EXTS) or p.endswith(_FASTA_EXTS)

def cifi_is_fastq(sample_name):
    """True if the canonical CiFi input needs conversion to BAM (FASTQ or FASTA)."""
    path = get_sample_data(sample_name)["cifi_input"]
    return _path_is_fastx(path)

def hifi_input_is_fasta(sample_name):
    """True if the sample uses `hifi_fasta:` (skip BAM->FASTA conversion)."""
    return bool(get_sample_data(sample_name)["hifi_fasta"])

def _canonical_cifi_bam_path(sample_name: str) -> str:
    """Upstream CiFi BAM: the original if BAM, else the fastq->bam conversion target."""
    if cifi_is_fastq(sample_name):
        return OUTDIR + f"/cifi/{sample_name}.from_fastq.bam"
    return get_sample_data(sample_name)["cifi_input"]

def get_cifi_bam(wildcards):
    """Return the upstream CiFi BAM for `downsample_cifi_bam` at this (sample, label).

    For external mode, each label maps to its own pre-downsampled BAM
    (populated in PER_SAMPLE_CIFI_SRC by `_scenarios_cifi_external`).
    Otherwise fall back to the canonical per-sample CiFi input.
    """
    key = (wildcards.label, wildcards.sample)
    if key in PER_SAMPLE_CIFI_SRC:
        return PER_SAMPLE_CIFI_SRC[key]
    return _canonical_cifi_bam_path(wildcards.sample)

# Script paths (relative to workflow)
CALN50_JS = "scripts/calN50.js"

# CiFi toolkit options from config (with defaults matching cifi CLI defaults)
CIFI_QC_OPTS = config.get("cifi", {}).get("qc", {})
CIFI_DIGEST_OPTS = config.get("cifi", {}).get("digest", {})

# Output directory (everything the pipeline writes lives under this)
OUTDIR = config.get("output_dir", "results")

# External tool paths from config
SINGULARITY_CACHE = config["tools"]["singularity_cache"]
THREED_DNA = os.path.abspath(config["tools"]["threed_dna"])
JUICER_TOOLS_JAR = os.path.abspath(config["tools"]["juicer_tools_jar"])

# SLURM settings from config (with sensible defaults)
SLURM_PARTITION = config.get("slurm", {}).get("partition", "low")
SLURM_ACCOUNT = config.get("slurm", {}).get("account", "publicgrp")

# ============================================================================
# BAM stats helpers (cached to disk so parse-time stays cheap)
# ============================================================================

def _bam_cache_path(bam_path: str, suffix: str) -> str:
    return f"{bam_path}.{suffix}"

def _cache_is_fresh(cache: str, source: str) -> bool:
    return (os.path.exists(cache)
            and os.path.exists(source)
            and os.path.getmtime(cache) >= os.path.getmtime(source))

def get_bam_total_bases(bam_path: str) -> int:
    """Return total sequence length (bp) of all reads in a BAM, with caching."""
    cache = _bam_cache_path(bam_path, "stats.bases")
    if _cache_is_fresh(cache, bam_path):
        return int(open(cache).read().strip())
    out = subprocess.check_output(["samtools", "stats", bam_path], text=True)
    total = 0
    for line in out.splitlines():
        if line.startswith("SN\ttotal length:"):
            total = int(line.split("\t")[2])
            break
    if total == 0:
        raise RuntimeError(f"samtools stats produced no 'total length' for {bam_path}")
    with open(cache, "w") as fh:
        fh.write(str(total))
    return total

# ============================================================================
# Scenario computation: HiFi/CiFi downsampling sweep
#
# A "scenario" is a (label, hifi_frac, cifi_frac) triple. The whole pipeline
# runs once per scenario. SCENARIOS is built once at parse time from config.
# Per-sample fractions are stored in PER_SAMPLE_FRACS keyed by (label, sample).
# ============================================================================

# (label, sample) -> (hifi_frac, cifi_frac)
PER_SAMPLE_FRACS: dict = {}

# (label, sample) -> pre-downsampled CiFi BAM path (only populated in external mode).
# When a key is present here, `get_cifi_bam` returns that path instead of the
# canonical per-sample CiFi input, and `downsample_cifi_bam` sees fraction=1.0
# (which triggers the symlink shortcut — no re-sampling).
PER_SAMPLE_CIFI_SRC: dict = {}

def _fmt_num(x) -> str:
    """Render an int or float without trailing zeros (e.g. 5 -> '5', 5.5 -> '5.5')."""
    if isinstance(x, int):
        return str(x)
    if float(x).is_integer():
        return str(int(x))
    return str(x).rstrip("0").rstrip(".") or "0"

def _depth_to_fraction(depth_x: float, genome_size_bp: int, total_bases: int) -> float:
    target = depth_x * genome_size_bp
    if target >= total_bases:
        sys.stderr.write(
            f"[cifiasm] WARNING: requested depth {depth_x}X needs {target} bp but BAM "
            f"has {total_bases} bp; clamping HiFi fraction to 1.0\n"
        )
        return 1.0
    return target / total_bases

def _hifi_value_label(hifi_cfg: dict, val) -> str:
    if hifi_cfg.get("mode") == "depth":
        return f"h{_fmt_num(val)}X"
    return f"h{_fmt_num(val)}pct"

def _hifi_value_to_frac(hifi_cfg: dict, val, sample_dict) -> float:
    if hifi_cfg.get("mode") == "depth":
        gs = hifi_cfg.get("genome_size")
        if gs is None:
            raise ValueError(
                "hifi.downsample.mode == 'depth' requires hifi.downsample.genome_size (bp)"
            )
        total = get_bam_total_bases(sample_dict["hifi_bam"])
        return _depth_to_fraction(float(val), int(gs), int(total))
    # fraction mode
    return float(val) / 100.0

def _hifi_values(hifi_cfg: dict):
    if hifi_cfg.get("mode") == "depth":
        return list(hifi_cfg.get("depths", []))
    return list(hifi_cfg.get("percentages", []))

def _scenarios_default():
    return [{"label": "100", "hifi_frac": 1.0, "cifi_frac": 1.0}]

def _scenarios_cifi_only(dil_cfg: dict, samples: list):
    out = []
    for pct in dil_cfg["percentages"]:
        label = str(pct)
        cfrac = float(pct) / 100.0
        for s in samples:
            PER_SAMPLE_FRACS[(label, s["sample"])] = (1.0, cfrac)
        out.append({"label": label, "hifi_frac": 1.0, "cifi_frac": cfrac})
    return out

def _scenarios_hifi_only(hifi_cfg: dict, samples: list):
    out = []
    for v in _hifi_values(hifi_cfg):
        label = _hifi_value_label(hifi_cfg, v)
        # store per-sample fracs (depth mode is sample-specific)
        last_frac = 1.0
        for s in samples:
            hf = _hifi_value_to_frac(hifi_cfg, v, s)
            PER_SAMPLE_FRACS[(label, s["sample"])] = (hf, 1.0)
            last_frac = hf
        out.append({"label": label, "hifi_frac": last_frac, "cifi_frac": 1.0})
    return out

def _scenarios_cifi_external(samples: list):
    """User-supplied pre-downsampled CiFi inputs, one per label per sample.

    Each file may be BAM, FASTQ, or FASTA (gzipped or not). `downsample_cifi_bam`
    converts FASTQ/FASTA to BAM via `samtools import` before the rest of the
    pipeline consumes it; BAMs are symlinked.

    All samples must declare the same label set.
    """
    missing = [s["sample"] for s in samples if not s["cifi_external"]]
    if missing:
        raise ValueError(
            f"cifi_external is set for some samples but missing for: {missing}. "
            f"All samples must declare cifi_external with the same labels when "
            f"using pre-downsampled CiFi inputs."
        )

    first_labels = tuple(samples[0]["cifi_external"].keys())
    for s in samples[1:]:
        this_labels = tuple(s["cifi_external"].keys())
        if set(this_labels) != set(first_labels):
            raise ValueError(
                f"cifi_external labels differ across samples: "
                f"'{samples[0]['sample']}' has {sorted(first_labels)}, "
                f"'{s['sample']}' has {sorted(this_labels)}."
            )

    out = []
    for label in first_labels:
        if not _LABEL_RE.match(label):
            raise ValueError(
                f"cifi_external label '{label}' must match {_LABEL_RE.pattern} "
                f"(labels appear in output paths)."
            )
        for s in samples:
            path = s["cifi_external"][label]
            if not (_path_is_bam(path) or _path_is_fastx(path)):
                raise ValueError(
                    f"cifi_external['{label}'] for sample '{s['sample']}' must be "
                    f"a BAM, FASTQ, or FASTA file; got '{path}'."
                )
            PER_SAMPLE_CIFI_SRC[(label, s["sample"])] = path
            PER_SAMPLE_FRACS[(label, s["sample"])] = (1.0, 1.0)
        out.append({"label": label, "hifi_frac": 1.0, "cifi_frac": 1.0})
    return out

def _scenarios_zip(hifi_cfg: dict, dil_cfg: dict, samples: list):
    h_values = _hifi_values(hifi_cfg)
    c_values = list(dil_cfg["percentages"])
    if len(h_values) != len(c_values):
        raise ValueError(
            f"zip mode requires equal-length lists, got "
            f"hifi.downsample ({len(h_values)}) and dilution.percentages ({len(c_values)})"
        )
    out = []
    for hv, cv in zip(h_values, c_values):
        label = f"{_hifi_value_label(hifi_cfg, hv)}_c{_fmt_num(cv)}"
        cfrac = float(cv) / 100.0
        last_hfrac = 1.0
        for s in samples:
            hfrac = _hifi_value_to_frac(hifi_cfg, hv, s)
            PER_SAMPLE_FRACS[(label, s["sample"])] = (hfrac, cfrac)
            last_hfrac = hfrac
        out.append({"label": label, "hifi_frac": last_hfrac, "cifi_frac": cfrac})
    return out

def build_scenarios(config_dict, samples):
    """Return the list of scenario dicts and populate PER_SAMPLE_FRACS."""
    PER_SAMPLE_FRACS.clear()
    PER_SAMPLE_CIFI_SRC.clear()
    hifi_cfg  = config_dict.get("hifi", {}).get("downsample", {})
    dil_cfg   = config_dict.get("dilution", {})

    hifi_on = hifi_cfg.get("enabled", False)
    cifi_on = dil_cfg.get("enabled", False)
    external_on = any(s["cifi_external"] for s in samples)

    if hifi_on:
        fasta_samples = [s["sample"] for s in samples if s["hifi_fasta"]]
        if fasta_samples:
            raise ValueError(
                f"hifi.downsample.enabled cannot be combined with `hifi_fasta` "
                f"(FASTA is not sub-sampleable by samtools view). Samples using "
                f"hifi_fasta: {fasta_samples}. Provide `hifi_bam` instead, or "
                f"disable hifi.downsample."
            )

    if external_on:
        if hifi_on or cifi_on:
            raise ValueError(
                "cifi_external (pre-downsampled CiFi) cannot be combined with "
                "hifi.downsample.enabled or dilution.enabled. Disable the "
                "fractional downsample modes when supplying pre-downsampled inputs."
            )
        return _scenarios_cifi_external(samples)

    if hifi_on and cifi_on:
        return _scenarios_zip(hifi_cfg, dil_cfg, samples)
    if hifi_on:
        return _scenarios_hifi_only(hifi_cfg, samples)
    if cifi_on:
        return _scenarios_cifi_only(dil_cfg, samples)

    # nothing enabled: backward-compat single scenario "100"
    for s in samples:
        PER_SAMPLE_FRACS[("100", s["sample"])] = (1.0, 1.0)
    return _scenarios_default()

SCENARIOS = build_scenarios(config, SAMPLES)
FRAC_LABELS = [s["label"] for s in SCENARIOS]

def get_hifi_frac_for(label: str, sample: str) -> float:
    return PER_SAMPLE_FRACS[(label, sample)][0]

def get_cifi_frac_for(label: str, sample: str) -> float:
    return PER_SAMPLE_FRACS[(label, sample)][1]

# ---------------- Targets ----------------
rule all:
    input:
        expand(OUTDIR + "/qc_cifi/{sample}/qc.pdf", sample=samples_list()),
        expand(OUTDIR + "/cifi/{sample}.{label}.bam",
               sample=samples_list(), label=FRAC_LABELS),
        expand(OUTDIR + "/stats/{sample}/{label}/summary.tsv",
               sample=samples_list(), label=FRAC_LABELS),
        expand(OUTDIR + "/stats/{sample}/{label}/yahs_summary.tsv",
               sample=samples_list(), label=FRAC_LABELS),
        expand(OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/bams/{sample}.{label}.hap{hap}.cs.bam",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        expand(OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        expand(OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),


# ---------------- Core steps ----------------

rule cifi_qc:
    """Run CiFi QC on raw CiFi input (BAM or FASTQ)"""
    priority: 10
    input:
        cifi=lambda w: get_sample_data(w.sample)["cifi_input"]
    output:
        pdf=OUTDIR + "/qc_cifi/{sample}/qc.pdf"
    params:
        outdir=OUTDIR + "/qc_cifi/{sample}",
        enzyme_args=get_enzyme_args,
        num_reads=CIFI_QC_OPTS.get("num_reads", 0),
        min_sites=CIFI_QC_OPTS.get("min_sites", 1),
    threads: 1
    resources:
        mem_mb=8000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "cifi qc {input.cifi} -o {params.outdir} {params.enzyme_args} "
        "-n {params.num_reads} --min-sites {params.min_sites}"

rule cifi_fastq_to_bam:
    """Convert canonical CiFi FASTQ or FASTA (gzipped or not) to an unmapped BAM."""
    priority: 50
    input:
        fq=lambda w: get_sample_data(w.sample)["cifi_input"]
    output:
        bam=OUTDIR + "/cifi/{sample}.from_fastq.bam"
    threads: 8
    resources:
        mem_mb=16000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools import -@ {threads} -0 {input.fq} -o {output.bam}"

rule downsample_hifi_bam:
    """Create label-specific HiFi BAM (fraction == 1.0 → symlink)."""
    priority: 150
    input:
        bam=lambda w: get_sample_data(w.sample)["hifi_bam"]
    output:
        bam=OUTDIR + "/hifi/{sample}.{label}.hifi.bam"
    params:
        frac=lambda w: get_hifi_frac_for(w.label, w.sample),
        sarg=lambda w: seeddotfrac_from_fraction(get_hifi_frac_for(w.label, w.sample), 100),
    threads: 4
    resources:
        mem_mb=16000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        if awk "BEGIN {{ exit !({params.frac} >= 0.9999999999) }}"; then
            ln -sf $(readlink -f {input.bam}) {output.bam}
        else
            samtools view -@ {threads} -b -s {params.sarg} -o {output.bam} {input.bam}
        fi
        '''

def _hifi_fasta_src(w):
    """Input resolver for `rule hifi_fasta`: user FASTA (when provided) or the
    labeled HiFi BAM produced by `downsample_hifi_bam`."""
    if hifi_input_is_fasta(w.sample):
        return get_sample_data(w.sample)["hifi_fasta"]
    return OUTDIR + f"/hifi/{w.sample}.{w.label}.hifi.bam"

rule hifi_fasta:
    """Produce per-label HiFi FASTA for hifiasm.

    BAM input → `samtools fasta`.
    FASTA input → symlink (or gunzip if .gz, since downstream tools want plain FASTA).
    """
    priority: 100
    input:
        src=_hifi_fasta_src
    output:
        fa=OUTDIR + "/hifi/{sample}.{label}.hifi.fa"
    threads: 8
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.fa})
        SRC="{input.src}"
        case "$SRC" in
          *.bam)
            samtools fasta -@ {threads} "$SRC" > {output.fa}
            ;;
          *.fa.gz|*.fasta.gz)
            gunzip -c "$SRC" > {output.fa}
            ;;
          *.fa|*.fasta)
            ln -sf "$(readlink -f "$SRC")" {output.fa}
            ;;
          *)
            echo "ERROR: unsupported HiFi input extension: $SRC" >&2
            exit 1
            ;;
        esac
        '''


rule downsample_cifi_bam:
    """Create label-specific CiFi BAM (for porec_nextflow and downstream FASTQ).

    Input path comes from `get_cifi_bam`:
      - canonical mode: sample's full CiFi BAM (or the from_fastq-converted BAM).
      - external mode: user's pre-downsampled BAM / FASTQ / FASTA for this label.

    BAM input: symlink if fraction=1.0, else `samtools view -s`.
    FASTQ/FASTA input: `samtools import` to BAM. Fractional downsampling is
    only supported for BAM sources.
    """
    priority: 350
    input:
        src=get_cifi_bam
    output:
        bam=OUTDIR + "/cifi/{sample}.{label}.bam",
        bai=OUTDIR + "/cifi/{sample}.{label}.bam.bai"
    params:
        frac=lambda w: get_cifi_frac_for(w.label, w.sample),
        sarg=lambda w: seeddotfrac_from_fraction(get_cifi_frac_for(w.label, w.sample), 100),
    threads: 4
    resources:
        mem_mb=16000, runtime=2 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        SRC="{input.src}"
        case "$SRC" in
          *.bam)
            if awk "BEGIN {{ exit !({params.frac} >= 0.9999999999) }}"; then
                ln -sf "$(readlink -f "$SRC")" {output.bam}
            else
                samtools view -@ {threads} -b -s {params.sarg} -o {output.bam} "$SRC"
            fi
            ;;
          *.fq|*.fq.gz|*.fastq|*.fastq.gz|*.fa|*.fa.gz|*.fasta|*.fasta.gz)
            if ! awk "BEGIN {{ exit !({params.frac} >= 0.9999999999) }}"; then
                echo "ERROR: fractional downsampling is only supported for BAM CiFi sources (got $SRC with frac={params.frac})" >&2
                exit 1
            fi
            samtools import -@ {threads} -0 "$SRC" -o {output.bam}
            ;;
          *)
            echo "ERROR: unsupported CiFi input extension: $SRC" >&2
            exit 1
            ;;
        esac
        samtools index -@ {threads} {output.bam}
        '''


rule cifi_fastq_from_downsampled_bam:
    """Extract FASTQ from downsampled CiFi BAM"""
    priority: 200
    input:
        bam=OUTDIR + "/cifi/{sample}.{label}.bam"
    output:
        fq=OUTDIR + "/cifi/{sample}.{label}.fastq"
    threads: 4
    resources:
        mem_mb=4*1024, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools collate -O -u -@ {threads} {input.bam} | "
        "samtools fastq   -@ {threads} - > {output.fq}"


rule cifi2pe_split:
    """CiFi single FASTQ -> HiC-like PE (R1/R2) with restriction enzyme"""
    priority: 400
    input:  OUTDIR + "/cifi/{sample}.{label}.fastq"
    output:
        r1=OUTDIR + "/cifi2pe/{sample}.{label}_R1.fastq",
        r2=OUTDIR + "/cifi2pe/{sample}.{label}_R2.fastq"
    params:
        out=OUTDIR + "/cifi2pe/{sample}.{label}",
        enzyme_args=get_enzyme_args,
        min_frags=CIFI_DIGEST_OPTS.get("min_fragments", 3),
        min_frag_len=CIFI_DIGEST_OPTS.get("min_frag_len", 20),
        extra=get_digest_extra_flags(),
    threads: 1
    resources:
        mem_mb=16000, runtime=12 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "cifi digest {input} {params.enzyme_args} -o {params.out} "
        "-m {params.min_frags} -l {params.min_frag_len} {params.extra}"

rule hifiasm_dual_scaf:
    """Assemble with hifiasm --dual-scaf (produces hap1/2 ctg GFAs)"""
    priority: 500
    input:
        r1=OUTDIR + "/cifi2pe/{sample}.{label}_R1.fastq",
        r2=OUTDIR + "/cifi2pe/{sample}.{label}_R2.fastq",
        hifi=OUTDIR + "/hifi/{sample}.{label}.hifi.fa"
    output:
        hap1_gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap1.p_ctg.gfa",
        hap2_gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap2.p_ctg.gfa"
    params:
        pref=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm"
    threads: 64
    resources:
        mem_mb=300*1024, runtime=24 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "hifiasm --dual-scaf --telo-m CCCTAAA -t {threads} -o {params.pref} "
        "--h1 {input.r1} --h2 {input.r2} {input.hifi}"

# ---------------- QC: GFA -> FASTA -> calN50 -> summary ----------------

rule gfa2fa:
    """Convert GFA contigs to FASTA (hap1 or hap2)"""
    priority: 600
    input:
        gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap{hap}.p_ctg.gfa"
    output:
        fa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    threads: 4
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "gfatools gfa2fa {input.gfa} > {output.fa}"

rule caln50:
    """Run calN50.js (k8) on each hap FASTA"""
    priority: 700
    input:
        fa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        n50=OUTDIR + "/stats/{sample}/{label}/hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"

rule summarize_assembly:
    """
    Parse k8 calN50.js output from hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 800
    input:
        hap1=OUTDIR + "/stats/{sample}/{label}/hap1.n50.txt",
        hap2=OUTDIR + "/stats/{sample}/{label}/hap2.n50.txt"
    output:
        tsv=OUTDIR + "/stats/{sample}/{label}/summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''


rule porec_nextflow:
    """Run epi2me-labs/wf-pore-c nextflow pipeline for Hi-C contact mapping"""
    priority: 650
    input:
        cifi_bam=OUTDIR + "/cifi/{sample}.{label}.bam",
        ref_fa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        bed=OUTDIR + "/porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed",
        pairs=OUTDIR + "/porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz",
        mcool=OUTDIR + "/porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.mcool",
        hic=OUTDIR + "/porec/{sample}/{label}/hap{hap}/hi-c/{sample}.{label}.hap{hap}.hic",
        report=OUTDIR + "/porec/{sample}/{label}/hap{hap}/wf-pore-c-report.html"
    params:
        outdir=OUTDIR + "/porec/{sample}/{label}/hap{hap}",
        sample_alias="{sample}.{label}.hap{hap}",
        cutter=get_enzyme,
        nxf_cache=SINGULARITY_CACHE,
        cifi_bam_abs=lambda w, input: os.path.abspath(input.cifi_bam),
        ref_fa_abs=lambda w, input: os.path.abspath(input.ref_fa),
        nf_log_abs=lambda w: os.path.abspath(
            f"{OUTDIR}/porec/{w.sample}/{w.label}/logs/{w.sample}.{w.label}.hap{w.hap}.nextflow.log"),
        time_log_abs=lambda w: os.path.abspath(
            f"{OUTDIR}/porec/{w.sample}/{w.label}/logs/{w.sample}.{w.label}.hap{w.hap}.time.txt"),
    log:
        nf=OUTDIR + "/porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.nextflow.log",
        time=OUTDIR + "/porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.time.txt"
    threads: 32
    resources:
        mem_mb= 500*1024,
        runtime=48 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        """
        export NXF_SINGULARITY_CACHEDIR="{params.nxf_cache}"
        export NXF_OPTS='-Xms2g -Xmx4g'
        export NXF_OFFLINE='true'
        export HDF5_USE_FILE_LOCKING='FALSE'
        export TMPDIR="/scratch/tmp/{wildcards.sample}.{wildcards.label}.hap{wildcards.hap}"
        mkdir -p $TMPDIR
        export NXF_TMP="$TMPDIR"
        export NXF_TEMP="$TMPDIR"
        export NXF_EXECUTOR='local'

        # Create and cd to output directory to isolate nextflow run
        mkdir -p {params.outdir}
        cd {params.outdir}
        echo "Running nextflow in $(pwd)"
        echo "Input BAM: {params.cifi_bam_abs}"
        echo "Input REF: {params.ref_fa_abs}"

        /usr/bin/time -v \
            nextflow run epi2me-labs/wf-pore-c \
                -r v1.3.0 \
                -profile singularity -resume \
                --bam {params.cifi_bam_abs} \
                --ref {params.ref_fa_abs} \
                --cutter {params.cutter} \
                --out_dir . \
                --threads {threads} \
                --minimap2_settings '-ax map-hifi' \
                --paired_end \
                --bed --pairs --hi_c --mcool --coverage \
                --sample "{params.sample_alias}" \
                -with-report   pipeline_report.html \
                -with-timeline pipeline_timeline.html \
                -with-trace    pipeline_trace.txt \
                -with-dag      pipeline_dag.svg \
                1> {params.nf_log_abs} 2> {params.time_log_abs}
        """

rule index_fa:
    """Index FASTA files for yahs"""
    priority: 625
    input:
        fa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        fai=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai"
    threads: 1
    resources:
        mem_mb=8 * 1024, runtime=60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools faidx {input.fa}"


rule yahs_scaffold:
    """Scaffold haplotype assemblies with yahs"""
    priority: 750
    input:
        asm_fa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa",
        asm_fai=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai",
        porec_bed=OUTDIR + "/porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed"
    output:
        scaffolds=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa",
        agp=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.agp",
        bin=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}.bin"
    params:
        prefix=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}"
    log:
        OUTDIR + "/yahs/{sample}/{label}/logs/{sample}.{label}.hap{hap}.yahs.log"
    threads: 32
    resources:
        mem_mb=250*1024, runtime=48 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        """
        # removed -q 0 --no-contig-ec
        yahs \
            -o {params.prefix} \
            -v 1 \
            {input.asm_fa} \
            {input.porec_bed} \
            2>&1 | tee {log}
        """

rule yahs_caln50:
    """Run calN50.js (k8) on each yahs scaffold FASTA"""
    priority: 850
    input:
        fa=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        n50=OUTDIR + "/stats/{sample}/{label}/yahs_hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"


rule summarize_yahs:
    """
    Parse k8 calN50.js output from yahs scaffolds hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 900
    input:
        hap1=OUTDIR + "/stats/{sample}/{label}/yahs_hap1.n50.txt",
        hap2=OUTDIR + "/stats/{sample}/{label}/yahs_hap2.n50.txt"
    output:
        tsv=OUTDIR + "/stats/{sample}/{label}/yahs_summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''

rule yahs_index_scaffolds_fa:
    """Index YAHS scaffold FASTA so we can derive chrom.sizes."""
    priority: 740
    input:
        fa=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        fai=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa.fai"
    threads: 1
    resources:
        mem_mb=8*1024, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools faidx {input.fa}"

rule qc_porec_nextflow:
    """Run epi2me-labs/wf-pore-c nextflow pipeline for Hi-C contact mapping on scaffolds"""
    priority: 650
    input:
        cifi_bam=OUTDIR + "/cifi/{sample}.{label}.bam",
        ref_fa=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        bed=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed",
        pairs=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz",
        mcool=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.mcool",
        hic=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/hi-c/{sample}.{label}.hap{hap}.hic",
        report=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/wf-pore-c-report.html",
        bam=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/bams/{sample}.{label}.hap{hap}.cs.bam",
        ns_bam=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/paired_end/{sample}.{label}.hap{hap}.ns.bam",
    params:
        outdir=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}",
        sample_alias="{sample}.{label}.hap{hap}",
        cutter=get_enzyme,
        nxf_cache=SINGULARITY_CACHE,
        cifi_bam_abs=lambda w, input: os.path.abspath(input.cifi_bam),
        ref_fa_abs=lambda w, input: os.path.abspath(input.ref_fa),
        nf_log_abs=lambda w: os.path.abspath(
            f"{OUTDIR}/qc_porec/{w.sample}/{w.label}/logs/{w.sample}.{w.label}.hap{w.hap}.nextflow.log"),
        time_log_abs=lambda w: os.path.abspath(
            f"{OUTDIR}/qc_porec/{w.sample}/{w.label}/logs/{w.sample}.{w.label}.hap{w.hap}.time.txt"),
    log:
        nf=OUTDIR + "/qc_porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.nextflow.log",
        time=OUTDIR + "/qc_porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.time.txt"
    threads: 32
    resources:
        mem_mb= 400*1024,
        runtime=48 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        """
        export NXF_SINGULARITY_CACHEDIR="{params.nxf_cache}"
        export NXF_OPTS='-Xms2g -Xmx4g'
        export NXF_OFFLINE='true'
        export HDF5_USE_FILE_LOCKING='FALSE'
        export TMPDIR="/scratch/tmp/{wildcards.sample}.{wildcards.label}.hap{wildcards.hap}"
        mkdir -p $TMPDIR
        export NXF_TMP="$TMPDIR"
        export NXF_TEMP="$TMPDIR"
        export NXF_EXECUTOR='local'

        # Create and cd to output directory to isolate nextflow run
        mkdir -p {params.outdir}
        cd {params.outdir}
        echo "Running nextflow in $(pwd)"
        echo "Input BAM: {params.cifi_bam_abs}"
        echo "Input REF: {params.ref_fa_abs}"

        /usr/bin/time -v \
            nextflow run epi2me-labs/wf-pore-c \
                -r v1.3.0 \
                -profile singularity -resume \
                --bam {params.cifi_bam_abs} \
                --ref {params.ref_fa_abs} \
                --cutter {params.cutter} \
                --out_dir . \
                --threads {threads} \
                --minimap2_settings '-ax map-hifi' \
                --paired_end \
                --bed --pairs --hi_c --mcool --coverage \
                --sample "{params.sample_alias}" \
                -with-report   pipeline_report.html \
                -with-timeline pipeline_timeline.html \
                -with-trace    pipeline_trace.txt \
                -with-dag      pipeline_dag.svg \
                1> {params.nf_log_abs} 2> {params.time_log_abs}
        """


# ============================================================================
# JBAT (Juicebox Assembly Tools) File Preparation
# ============================================================================
# Prerequisites:
#   git clone https://github.com/aidenlab/3d-dna.git
#   wget https://github.com/aidenlab/JuicerTools/releases/download/v3.0.0/juicer_tools.jar
#
# Add to config.yaml:
#   tools:
#     threed_dna: "/path/to/3d-dna"
#     juicer_tools_jar: "/path/to/juicer_tools.jar"
# ============================================================================

# ----------------------------------------------------------------------------
# Pre-processing: Filter BAM and convert to pairs
# ----------------------------------------------------------------------------

rule cifi_filter_bam:
    """
    Filter the paired-end BAM from qc_porec to remove unmapped reads etc.
    """
    input:
        bam=OUTDIR + "/qc_porec/{sample}/{label}/hap{hap}/paired_end/{sample}.{label}.hap{hap}.ns.bam"
    output:
        bam=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.filtered.bam"
    params:
        outdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}"
    threads: 8
    resources:
        mem_mb=32*1024, runtime=4*60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p {params.outdir}
        cifi filter {input.bam} -t {threads} -o {output.bam}
        '''


rule generate_genome_file:
    """
    Generate .genome file (chrom sizes) from yahs scaffold FASTA index.
    """
    input:
        fai=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa.fai"
    output:
        genome=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.genome"
    params:
        outdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}"
    threads: 1
    resources:
        mem_mb=1024, runtime=10, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p {params.outdir}
        cut -f1,2 {input.fai} > {output.genome}
        '''


rule bam2pairs:
    """
    Convert filtered BAM to pairs format using pairix bam2pairs.
    """
    input:
        bam=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.filtered.bam",
        genome=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.genome"
    output:
        pairs=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.bam2pairs.bsorted.pairs.gz"
    params:
        prefix=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.bam2pairs"
    threads: 8
    resources:
        mem_mb=32*1024, runtime=4*60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mamba run -n pairix bam2pairs -l -c {input.genome} {input.bam} {params.prefix}
        '''


rule pairs_to_mnd:
    """
    Convert .pairs.gz from bam2pairs to merged_nodups.txt format for 3D-DNA.

    Short format (11 columns):
    str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 mapq2

    Only include pairs where both mates are mapped (filter out '!' entries).
    """
    input:
        pairs=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.bam2pairs.bsorted.pairs.gz"
    output:
        mnd=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    params:
        outdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}"
    threads: 4
    resources:
        mem_mb=32*1024, runtime=4*60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p {params.outdir}
        # Convert pairs format to merged_nodups.txt long format (16 columns)
        # Input pairs columns: readID chr1 pos1 chr2 pos2 strand1 strand2 ...
        # Output mnd: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 cigar2 seq2 readname1 readname2
        # Filter: only include pairs where BOTH chr1 and chr2 are mapped (not "!" or ".")
        zcat {input.pairs} | awk 'BEGIN{{OFS=" "}}
            !/^#/ && $2 != "!" && $2 != "." && $4 != "!" && $4 != "." {{
                str1 = ($6 == "+") ? 0 : 16
                str2 = ($7 == "+") ? 0 : 16
                print str1, $2, $3, 0, str2, $4, $5, 1, 60, "-", "-", 60, "-", "-", $1, $1
            }}' > {output.mnd}
        '''


rule generate_assembly_file:
    """
    Generate .assembly file from yahs scaffold FASTA using 3D-DNA utility.
    """
    input:
        fa=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        assembly=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly"
    params:
        awk_script=f"{THREED_DNA}/utils/generate-assembly-file-from-fasta.awk",
        outdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}"
    threads: 1
    resources:
        mem_mb=8*1024, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p {params.outdir}
        awk -f {params.awk_script} {input.fa} > {output.assembly}
        '''


rule jbat_hic:
    """
    Generate JBAT-compatible .hic file using 3D-DNA visualizer.
    This creates the 'assembly' pseudo-chromosome format required by JBAT.
    """
    input:
        assembly=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly",
        mnd=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    output:
        hic=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic"
    params:
        threed_dna=THREED_DNA,
        juicer_jar=JUICER_TOOLS_JAR,
        workdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}",
        prefix="{sample}.{label}.hap{hap}"
    log:
        OUTDIR + "/jbat/{sample}/{label}/hap{hap}/logs/visualize.log"
    threads: 16
    resources:
        mem_mb=128*1024, runtime=12*60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        mkdir -p {params.workdir}/logs
        cd {params.workdir}

        # Set juicer_tools path for 3D-DNA scripts
        export JUICER_TOOLS_JAR="{params.juicer_jar}"

        # Run 3D-DNA visualizer (disable GNU Parallel with -p false)
        # Note: run-assembly-visualizer.sh exits with code 1 due to its
        # cleanup guard: [ "$clean_up" == "true" ] && rm ... returns 1
        # when clean_up is false (default). The .hic file is valid, so we
        # capture the exit code and only fail on real errors.
        set +o pipefail
        bash {params.threed_dna}/visualize/run-assembly-visualizer.sh \
            -p false \
            {params.prefix}.assembly \
            merged_nodups.txt \
            2>&1 | tee logs/visualize.log
        set -o pipefail

        # Verify the .hic file was actually created
        if [ ! -s {params.prefix}.hic ]; then
            echo "ERROR: .hic file was not created" >&2
            exit 1
        fi
        '''


rule jbat_post_review:
    """
    After manual curation in JBAT, convert reviewed .assembly to FASTA.

    Usage:
      1. Open results/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic in Juicebox
      2. Load the .assembly file via Assembly -> Import Map Assembly
      3. Make corrections, then Export Assembly as {sample}.{label}.hap{hap}.review.assembly
      4. Run: snakemake results/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.FINAL.fa
    """
    input:
        review=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.review.assembly",
        orig_fa=OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa",
        mnd=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    output:
        final_fa=OUTDIR + "/jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.FINAL.fa"
    params:
        threed_dna=THREED_DNA,
        juicer_jar=JUICER_TOOLS_JAR,
        workdir=OUTDIR + "/jbat/{sample}/{label}/hap{hap}",
        orig_fa_abs=lambda w, input: os.path.abspath(input.orig_fa),
    log:
        OUTDIR + "/jbat/{sample}/{label}/hap{hap}/logs/post_review.log"
    threads: 16
    resources:
        mem_mb=128*1024, runtime=12*60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        cd {params.workdir}
        export JUICER_TOOLS_JAR="{params.juicer_jar}"

        bash {params.threed_dna}/run-asm-pipeline-post-review.sh \
            -r {wildcards.sample}.{wildcards.label}.hap{wildcards.hap}.review.assembly \
            {params.orig_fa_abs} \
            merged_nodups.txt \
            2>&1 | tee logs/post_review.log

        # 3D-DNA outputs .FINAL.fasta, rename to .FINAL.fa for consistency
        if [ -f *.FINAL.fasta ]; then
            mv *.FINAL.fasta {wildcards.sample}.{wildcards.label}.hap{wildcards.hap}.FINAL.fa
        fi
        '''
