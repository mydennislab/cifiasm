rule jbat_pre:
    """Prepare each haplotype for Juicebox curation."""
    priority: 760
    input:
        bin=YAHS_PREFIX + ".bin",
        agp=YAHS_PREFIX + "_scaffolds_final.agp",
        fai=ASM_FA + ".fai"
    output:
        txt=JBAT_PREFIX + ".txt",
        assembly=JBAT_PREFIX + ".assembly",
        liftover=JBAT_PREFIX + ".liftover.agp",
        assembly_agp=JBAT_PREFIX + ".assembly.agp",
        chrom_sizes=JBAT_PREFIX + ".chrom.sizes",
        scale_factor=JBAT_PREFIX + ".scale_factor.txt"
    params:
        mapq=YAHS_MAPQ,
        prefix=JBAT_PREFIX,
        outdir=JBAT_DIR
    log:
        JBAT_PREFIX + ".log"
    threads: 2
    resources:
        mem_mb=16*1024, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        mkdir -p {params.outdir}
        juicer pre -a -q {params.mapq} -o {params.prefix} {input.bin} {input.agp} {input.fai} > {log} 2>&1
        awk 'BEGIN{{OFS="\t"}} /^PRE_C_SIZE:/ {{print $2, $3}}' {log} > {output.chrom_sizes}
        if [ ! -s {output.chrom_sizes} ]; then
            echo "ERROR: no PRE_C_SIZE line in {log}" >&2
            exit 1
        fi
        sf=$(sed -n 's/.*scale factor: \([0-9][0-9]*\).*/\1/p' {log} | tail -n 1)
        echo "${{sf:-1}}" > {output.scale_factor}
        '''

rule jbat_hic:
    """Generate the Juicebox contact map."""
    priority: 770
    input:
        txt=JBAT_PREFIX + ".txt",
        chrom_sizes=JBAT_PREFIX + ".chrom.sizes"
    output:
        hic=JBAT_PREFIX + ".hic"
    params:
        jar=JUICER_TOOLS_JAR,
        heap_mb=lambda w, resources: int(resources.mem_mb * 0.9)
    log:
        JBAT_PREFIX + ".hic.log"
    threads: 8
    resources:
        mem_mb=64*1024, runtime=12 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        rm -f {output.hic}.part
        java -Xmx{params.heap_mb}m -jar {params.jar} pre {input.txt} {output.hic}.part {input.chrom_sizes} > {log} 2>&1 \
          && mv {output.hic}.part {output.hic}
        '''

rule jbat_post:
    """Build the curated FASTA and AGP from a reviewed assembly."""
    priority: 780
    input:
        review=JBAT_PREFIX + ".review.assembly",
        liftover=JBAT_PREFIX + ".liftover.agp",
        fa=ASM_FA
    output:
        final_fa=JBAT_PREFIX + ".FINAL.fa",
        final_agp=JBAT_PREFIX + ".FINAL.agp"
    params:
        prefix=JBAT_PREFIX
    log:
        JBAT_PREFIX + ".post.log"
    threads: 2
    resources:
        mem_mb=16*1024, runtime=4 * 60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        juicer post -o {params.prefix} {input.review} {input.liftover} {input.fa} > {log} 2>&1
        '''
