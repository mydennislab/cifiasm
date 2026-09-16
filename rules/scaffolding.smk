rule map_segments:
    """Align CiFi segments to each haplotype for contact calling."""
    priority: 650
    input:
        fa=ASM_FA,
        segments=OUTDIR + "/cifi2pe/{sample}.{label}.segments.fastq.gz"
    output:
        bam=OUTDIR + "/mapping/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.segments.ns.bam"
    # Name sorting groups segments for contact calling.
    params:
        preset=MINIMAP2_PRESET,
        mm2_threads=lambda w, threads: max(1, threads - MAP_SORT_THREADS),
        sort_threads=MAP_SORT_THREADS,
        sort_mem="2G",
        tmp_prefix=OUTDIR + "/mapping/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.sort_tmp"
    log:
        OUTDIR + "/mapping/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.minimap2.log"
    threads: 16 + MAP_SORT_THREADS
    resources:
        mem_mb=48*1024, runtime=6 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        minimap2 -t {params.mm2_threads} -ax {params.preset} --secondary=no --no-hash-name {input.fa} {input.segments} 2> {log} \
          | samtools sort -n -@ {params.sort_threads} -m {params.sort_mem} -T {params.tmp_prefix} -o {output.bam}
        '''

rule cifi_contacts:
    """Generate scaffolding contacts and reports from segment alignments."""
    priority: 675
    input:
        bam=OUTDIR + "/mapping/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.segments.ns.bam"
    output:
        bed=OUTDIR + "/contacts/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.bed",
        stats=OUTDIR + "/contacts/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}_contacts_stats.json",
        report=OUTDIR + "/contacts/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}_contacts_report.html"
    params:
        mapq=CONTACTS_MAPQ
    threads: 4
    resources:
        mem_mb=8*1024, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "cifi contacts {input.bam} --format bed -q {params.mapq} -o {output.bed} -t {threads}"


rule yahs_scaffold:
    """Scaffold each haplotype using CiFi contacts."""
    priority: 750
    input:
        asm_fa=ASM_FA,
        asm_fai=ASM_FA + ".fai",
        bed=OUTDIR + "/contacts/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.bed"
    output:
        scaffolds=YAHS_PREFIX + "_scaffolds_final.fa",
        agp=YAHS_PREFIX + "_scaffolds_final.agp",
        bin=YAHS_PREFIX + ".bin"
    params:
        prefix=YAHS_PREFIX,
        args=YAHS_ARGS
    log:
        OUTDIR + "/yahs/{sample}/{label}/logs/{sample}.{label}.hap{hap}.yahs.log"
    threads: 1
    resources:
        mem_mb=32768, runtime=24 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        yahs {params.args} -o {params.prefix} -v 1 {input.asm_fa} {input.bed} 2>&1 | tee {log}
        '''

rule yahs_caln50:
    """Calculate scaffold assembly statistics."""
    priority: 850
    input:
        fa=YAHS_PREFIX + "_scaffolds_final.fa"
    output:
        n50=OUTDIR + "/stats/{sample}/{label}/yahs_hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"


rule summarize_yahs:
    """Combine haplotype scaffold statistics into one table."""
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
