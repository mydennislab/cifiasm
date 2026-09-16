rule cifi_qc:
    """Generate per-sample CiFi QC reports."""
    priority: 10
    input:
        cifi=OUTDIR + "/cifi/merged/{sample}.cifi.bam"
    output:
        html=OUTDIR + "/qc_cifi/{sample}/qc.html",
        json=OUTDIR + "/qc_cifi/{sample}/qc.json"
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

rule merge_cifi:
    """Combine each sample's CiFi inputs into an unmapped BAM."""
    priority: 50
    input:
        files=lambda w: get_sample_data(w.sample)["cifi_files"]
    output:
        bam=OUTDIR + "/cifi/merged/{sample}.cifi.bam"
    threads: 8
    resources:
        mem_mb=16000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        files=( {input.files} )
        if [ "${{#files[@]}}" -eq 1 ]; then
            case "${{files[0]}}" in
              *.bam) ln -sf "$(readlink -f "${{files[0]}}")" {output.bam} ;;
              *)     samtools import -@ {threads} -0 "${{files[0]}}" -o {output.bam} ;;
            esac
        else
            TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
            bams=(); i=0
            for f in "${{files[@]}}"; do
              case "$f" in
                *.bam) bams+=("$f") ;;
                *)     imp="$TMP/imp_$i.bam"; samtools import -@ {threads} -0 "$f" -o "$imp"; bams+=("$imp"); i=$((i+1)) ;;
              esac
            done
            samtools cat -@ {threads} -o {output.bam} "${{bams[@]}}"
        fi
        '''

def get_hifi_inputs(sample_name):
    """Resolve HiFi inputs, using FASTQ conversions for BAM files."""
    out = []
    for idx, path in enumerate(get_sample_data(sample_name)["hifi_files"]):
        if _classify_input(path) == "bam":
            out.append(OUTDIR + f"/hifi/{sample_name}/cell{idx}.fastq")
        else:
            out.append(path)
    return out

def _hifi_bam_for_idx(w):
    """Select the HiFi BAM for a conversion job."""
    return get_sample_data(w.sample)["hifi_files"][int(w.idx)]

rule hifi_bam_to_fastq:
    """Convert a HiFi BAM input to FASTQ."""
    priority: 150
    wildcard_constraints:
        idx=r"[0-9]+"
    input:
        bam=_hifi_bam_for_idx
    output:
        fq=OUTDIR + "/hifi/{sample}/cell{idx}.fastq"
    threads: 4
    resources:
        mem_mb=16000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools fastq -@ {threads} {input.bam} > {output.fq}"


rule downsample_cifi_bam:
    """Prepare the CiFi BAM for one sample and label."""
    priority: 350
    input:
        src=get_cifi_bam
    output:
        bam=OUTDIR + "/cifi/{sample}.{label}.bam"
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
        '''


rule cifi_digest:
    """Generate CiFi read pairs, segments, and digestion reports."""
    priority: 400
    input:
        bam=OUTDIR + "/cifi/{sample}.{label}.bam"
    output:
        r1=OUTDIR + "/cifi2pe/{sample}.{label}_R1" + DIGEST_PE_EXT,
        r2=OUTDIR + "/cifi2pe/{sample}.{label}_R2" + DIGEST_PE_EXT,
        segments=OUTDIR + "/cifi2pe/{sample}.{label}.segments.fastq.gz",
        stats=OUTDIR + "/cifi2pe/{sample}.{label}_stats.json",
        report=OUTDIR + "/cifi2pe/{sample}.{label}_digestion_report.html"
    params:
        out=OUTDIR + "/cifi2pe/{sample}.{label}",
        enzyme_args=get_enzyme_args,
        min_segments=DIGEST_MIN_SEGMENTS,
        min_segment_len=DIGEST_MIN_SEGMENT_LEN,
        extra=get_digest_extra_flags(),
    threads: 1
    resources:
        mem_mb=32000, runtime=12 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "cifi digest {input.bam} {params.enzyme_args} -o {params.out} "
        "-m {params.min_segments} -l {params.min_segment_len} "
        "--segments-out {output.segments} {params.extra}"
