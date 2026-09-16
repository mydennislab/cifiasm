import os
import re
import sys

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/.]+",
    label=r"[A-Za-z0-9._-]+",
    hap=r"[12]"

_LABEL_RE = re.compile(r"^[A-Za-z0-9._-]+$")

def get_account_for_jobs(wildcards):
    return SLURM_ACCOUNT

def seeddotfrac_from_fraction(frac: float, seed: int = 100) -> str:
    """Format a seeded CiFi sampling fraction."""
    frac_digits = f"{frac:.10f}".split(".")[1].rstrip("0") or "0"
    return f"{seed}.{frac_digits}"

_FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
_FASTA_EXTS = (".fasta", ".fa", ".fasta.gz", ".fa.gz")

def _as_list(x):
    """Normalize sample paths to a list."""
    if not x:
        return []
    return [str(p) for p in x] if isinstance(x, (list, tuple)) else [str(x)]

def _classify_input(path):
    """Identify the input format from its extension."""
    p = str(path).lower()
    if p.endswith(".bam"):
        return "bam"
    if p.endswith(_FASTQ_EXTS):
        return "fastq"
    if p.endswith(_FASTA_EXTS):
        return "fasta"
    raise ValueError(f"Unrecognized input extension (need BAM/FASTQ/FASTA): {path}")

# Sample inputs
SAMPLES = []
for sample_name, sample_data in config["samples"].items():
    hifi_files = (_as_list(sample_data.get("hifi"))
                  + _as_list(sample_data.get("hifi_bam"))
                  + _as_list(sample_data.get("hifi_fasta")))
    if not hifi_files:
        raise ValueError(
            f"Sample '{sample_name}': no HiFi input found. Set `hifi:` "
            f"(a path or list of BAM/FASTQ/FASTA files)."
        )
    for p in hifi_files:
        _classify_input(p)

    cifi_files = _as_list(sample_data.get("cifi") or sample_data.get("cifi_bam"))
    if not cifi_files:
        raise ValueError(
            f"Sample '{sample_name}': no CiFi input found. Set `cifi:` "
            f"(a path or list of BAM/FASTQ/FASTA files)."
        )
    for p in cifi_files:
        _classify_input(p)

    SAMPLES.append(dict(
        sample=sample_name,
        hifi_files=hifi_files,
        cifi_files=cifi_files,
        enzyme=sample_data.get("enzyme", ""),
        site=sample_data.get("site", ""),
        cut_pos=sample_data.get("cut_pos", ""),
        # Normalize numeric YAML labels for output paths.
        cifi_external={str(k): v for k, v in (sample_data.get("cifi_external") or {}).items()},
    ))

def samples_list():
    return [s["sample"] for s in SAMPLES]

def get_sample_data(sample_name):
    return next(s for s in SAMPLES if s["sample"] == sample_name)

def get_enzyme(wildcards):
    return get_sample_data(wildcards.sample)["enzyme"]

def get_enzyme_args(wildcards):
    """Use a custom site when configured; otherwise use the named enzyme."""
    data = get_sample_data(wildcards.sample)
    if data.get("site"):
        return f"--site {data['site']} --cut-pos {data['cut_pos']}"
    return f"-e {data['enzyme']}"

def get_digest_extra_flags():
    """Read optional digestion flags from the configuration."""
    flags = []
    if not CIFI_DIGEST_OPTS.get("strip_overhang", True):
        flags.append("--no-strip-overhang")
    if CIFI_DIGEST_OPTS.get("revcomp_r2", False):
        flags.append("--revcomp-r2")
    if CIFI_DIGEST_OPTS.get("gzip", False):
        flags.append("--gzip")
    if CIFI_DIGEST_OPTS.get("fast", False):
        flags.append("--fast")
    return " ".join(flags)

def _canonical_cifi_bam_path(sample_name: str) -> str:
    """Return the merged CiFi BAM path."""
    return OUTDIR + f"/cifi/merged/{sample_name}.cifi.bam"

def get_cifi_bam(wildcards):
    """Select the CiFi input for this sample and label."""
    key = (wildcards.label, wildcards.sample)
    if key in PER_SAMPLE_CIFI_SRC:
        return PER_SAMPLE_CIFI_SRC[key]
    return _canonical_cifi_bam_path(wildcards.sample)

CALN50_JS = "scripts/calN50.js"

# Tool settings
HIFIASM_OPTS = config.get("hifiasm", {}) or {}
HIFIASM_TELOMERE_MOTIF = HIFIASM_OPTS.get("telomere_motif", "CCCTAAA")
if HIFIASM_TELOMERE_MOTIF is not None and not isinstance(HIFIASM_TELOMERE_MOTIF, str):
    raise ValueError("hifiasm.telomere_motif must be a string or null")

CIFI_QC_OPTS = config.get("cifi", {}).get("qc", {})
CIFI_DIGEST_OPTS = config.get("cifi", {}).get("digest", {})

# Current digest keys take precedence over legacy aliases.
DIGEST_MIN_SEGMENTS = CIFI_DIGEST_OPTS.get(
    "min_segments", CIFI_DIGEST_OPTS.get("min_fragments", 3))
DIGEST_MIN_SEGMENT_LEN = CIFI_DIGEST_OPTS.get(
    "min_segment_len", CIFI_DIGEST_OPTS.get("min_frag_len", 20))
DIGEST_PE_EXT = ".fastq.gz" if CIFI_DIGEST_OPTS.get("gzip", False) else ".fastq"

MAPPING_OPTS = config.get("mapping", {})
MINIMAP2_PRESET = MAPPING_OPTS.get("minimap2_preset", "map-hifi")
CONTACTS_MAPQ = int(MAPPING_OPTS.get("contacts_mapq", 1))
# Sorting shares the mapping job's thread allocation.
MAP_SORT_THREADS = 4

# YaHS and contact maps use the same additional MAPQ threshold.
YAHS_OPTS = config.get("yahs", {})
YAHS_MAPQ = int(YAHS_OPTS.get("mapq", 0))

# Contig error correction is disabled by default for CiFi scaffolding.
YAHS_CONTIG_EC = bool(YAHS_OPTS.get("contig_ec", False))
YAHS_ARGS = f"-q {YAHS_MAPQ}" + ("" if YAHS_CONTIG_EC else " --no-contig-ec")

OUTDIR = config.get("output_dir", "results")

JUICER_TOOLS_JAR = os.path.abspath(config["tools"]["juicer_tools_jar"])

OBSOLETE_TOOL_KEYS = ("singularity_cache", "threed_dna")
for _key in OBSOLETE_TOOL_KEYS:
    if _key in config.get("tools", {}):
        print(f"warning: tools.{_key} is no longer used and is ignored", file=sys.stderr)

# Optional outputs for the default target.
CONTACT_MAPS = config.get("contact_maps", {}) or {}
CM_JUICEBOX = bool(CONTACT_MAPS.get("juicebox", True))
CM_PRETEXT = bool(CONTACT_MAPS.get("pretext", True))
CM_PAIRS = bool(CONTACT_MAPS.get("pairs", False))
CM_SNAPSHOT = bool(CONTACT_MAPS.get("snapshot", False))
# Sorting shares the pairs job's thread allocation.
PAIRS_SORT_THREADS = 2

SLURM_PARTITION = config.get("slurm", {}).get("partition", "low")
SLURM_ACCOUNT = config.get("slurm", {}).get("account", "publicgrp")


# CiFi scenarios: (label, sample) -> sampling fraction.
PER_SAMPLE_FRACS: dict = {}

# External CiFi inputs: (label, sample) -> input path.
PER_SAMPLE_CIFI_SRC: dict = {}

def _scenarios_default():
    return [{"label": "100", "cifi_frac": 1.0}]

def _scenarios_cifi_only(dil_cfg: dict, samples: list):
    out = []
    for pct in dil_cfg["percentages"]:
        label = str(pct)
        cfrac = float(pct) / 100.0
        for s in samples:
            PER_SAMPLE_FRACS[(label, s["sample"])] = cfrac
        out.append({"label": label, "cifi_frac": cfrac})
    return out

def _scenarios_cifi_external(samples: list):
    """Validate external CiFi inputs and their shared labels."""
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
            try:
                _classify_input(path)
            except ValueError:
                raise ValueError(
                    f"cifi_external['{label}'] for sample '{s['sample']}' must be "
                    f"a BAM, FASTQ, or FASTA file; got '{path}'."
                )
            PER_SAMPLE_CIFI_SRC[(label, s["sample"])] = path
            PER_SAMPLE_FRACS[(label, s["sample"])] = 1.0
        out.append({"label": label, "cifi_frac": 1.0})
    return out

def build_scenarios(config_dict, samples):
    """Build CiFi scenarios and record each sample's sampling fraction."""
    PER_SAMPLE_FRACS.clear()
    PER_SAMPLE_CIFI_SRC.clear()
    dil_cfg = config_dict.get("dilution", {})

    cifi_on = dil_cfg.get("enabled", False)
    external_on = any(s["cifi_external"] for s in samples)

    if external_on:
        if cifi_on:
            raise ValueError(
                "cifi_external (pre-downsampled CiFi) cannot be combined with "
                "dilution.enabled. Disable the dilution sweep when supplying "
                "pre-downsampled inputs."
            )
        return _scenarios_cifi_external(samples)

    if cifi_on:
        return _scenarios_cifi_only(dil_cfg, samples)

    for s in samples:
        PER_SAMPLE_FRACS[("100", s["sample"])] = 1.0
    return _scenarios_default()

SCENARIOS = build_scenarios(config, SAMPLES)
FRAC_LABELS = [s["label"] for s in SCENARIOS]

def get_cifi_frac_for(label: str, sample: str) -> float:
    return PER_SAMPLE_FRACS[(label, sample)]

# Shared output paths
ASM_FA = OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
YAHS_PREFIX = OUTDIR + "/yahs/{sample}/{label}/{sample}.{label}.hap{hap}"
JBAT_DIR = OUTDIR + "/jbat/{sample}/{label}/hap{hap}"
JBAT_PREFIX = JBAT_DIR + "/{sample}.{label}.hap{hap}"
CM_DIR = OUTDIR + "/contact_maps/{sample}/{label}/hap{hap}"
CM_PREFIX = CM_DIR + "/{sample}.{label}.hap{hap}"

def contact_map_targets():
    """Collect the enabled contact-map outputs."""
    per_hap = dict(sample=samples_list(), label=FRAC_LABELS, hap=[1, 2])
    out = []
    if CM_JUICEBOX:
        out += expand(JBAT_PREFIX + ".assembly", **per_hap)
        out += expand(JBAT_PREFIX + ".hic", **per_hap)
    if CM_PRETEXT or CM_PAIRS:
        out += expand(CM_PREFIX + ".scaffolds.pairs.gz", **per_hap)
    if CM_PRETEXT:
        out += expand(CM_PREFIX + ".pretext", **per_hap)
    if CM_SNAPSHOT:
        out += expand(CM_PREFIX + ".pretext.png", **per_hap)
    return out

rule all:
    input:
        expand(OUTDIR + "/qc_cifi/{sample}/qc.html", sample=samples_list()),
        expand(OUTDIR + "/cifi/{sample}.{label}.bam",
               sample=samples_list(), label=FRAC_LABELS),
        expand(OUTDIR + "/stats/{sample}/{label}/summary.tsv",
               sample=samples_list(), label=FRAC_LABELS),
        expand(OUTDIR + "/stats/{sample}/{label}/yahs_summary.tsv",
               sample=samples_list(), label=FRAC_LABELS),
        expand(YAHS_PREFIX + "_scaffolds_final.fa",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        expand(YAHS_PREFIX + "_scaffolds_final.agp",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        expand(YAHS_PREFIX + ".bin",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        contact_map_targets(),


include: "rules/inputs.smk"

include: "rules/assembly.smk"


include: "rules/scaffolding.smk"


include: "rules/juicebox.smk"


include: "rules/contact_maps.smk"
