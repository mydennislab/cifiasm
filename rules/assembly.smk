rule hifiasm_dual_scaf:
    """Assemble phased contigs from HiFi reads and CiFi pairs."""
    priority: 500
    input:
        r1=OUTDIR + "/cifi2pe/{sample}.{label}_R1" + DIGEST_PE_EXT,
        r2=OUTDIR + "/cifi2pe/{sample}.{label}_R2" + DIGEST_PE_EXT,
        hifi=lambda w: get_hifi_inputs(w.sample)
    output:
        hap1_gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap1.p_ctg.gfa",
        hap2_gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap2.p_ctg.gfa"
    benchmark:
        OUTDIR + "/benchmarks/hifiasm_dual_scaf/{sample}/{label}.tsv"
    params:
        pref=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm",
        telomere_args=["--telo-m", HIFIASM_TELOMERE_MOTIF] if HIFIASM_TELOMERE_MOTIF else []
    threads: get_threads("hifiasm", 64)
    resources:
        mem_mb=300*1024, runtime=24 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "hifiasm --dual-scaf {params.telomere_args:q} -t {threads} -o {params.pref} "
        "--h1 {input.r1} --h2 {input.r2} {input.hifi}"


rule gfa2fa:
    """Export each haplotype's contigs as FASTA."""
    priority: 600
    input:
        gfa=OUTDIR + "/asm/{sample}/{label}/{sample}.{label}.asm.hic.hap{hap}.p_ctg.gfa"
    output:
        fa=ASM_FA
    threads: 4
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "gfatools gfa2fa {input.gfa} > {output.fa}"

rule index_fa:
    """Index the contig FASTA."""
    priority: 625
    input:
        fa=ASM_FA
    output:
        fai=ASM_FA + ".fai"
    threads: 1
    resources:
        mem_mb=8 * 1024, runtime=60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools faidx {input.fa}"

rule caln50:
    """Calculate contig assembly statistics."""
    priority: 700
    input:
        fa=ASM_FA
    output:
        n50=OUTDIR + "/stats/{sample}/{label}/hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"

rule summarize_assembly:
    """Combine haplotype contig statistics into one table."""
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
