rule index_scaffolds_fa:
    """Index the scaffold FASTA."""
    priority: 755
    input:
        fa=YAHS_PREFIX + "_scaffolds_final.fa"
    output:
        fai=YAHS_PREFIX + "_scaffolds_final.fa.fai"
    threads: 1
    resources:
        mem_mb=4 * 1024, runtime=60, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        "samtools faidx {input.fa}"

rule scaffold_pairs:
    """Write scaffold-coordinate contacts in 4DN pairs format."""
    # Export 1-based positions and reject scaled coordinates.
    priority: 765
    input:
        bin=YAHS_PREFIX + ".bin",
        agp=YAHS_PREFIX + "_scaffolds_final.agp",
        ctg_fai=ASM_FA + ".fai",
        scaf_fai=YAHS_PREFIX + "_scaffolds_final.fa.fai"
    output:
        pairs=CM_PREFIX + ".scaffolds.pairs.gz"
    params:
        script=os.path.join(workflow.basedir, "scripts", "scaffold_pairs.sh"),
        mapq=YAHS_MAPQ,
        outdir=CM_DIR,
        sort_tmp=CM_PREFIX + ".sort_tmp",
        sort_threads=PAIRS_SORT_THREADS,
        sort_mem=lambda w, resources: f"{int(resources.mem_mb * 0.5)}M"
    log:
        CM_PREFIX + ".scaffolds.pairs.log"
    threads: 2 + PAIRS_SORT_THREADS
    resources:
        mem_mb=2 * 1024, runtime=30, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        bash {params.script:q} \
            --bin {input.bin:q} --agp {input.agp:q} \
            --ctg-fai {input.ctg_fai:q} --scaf-fai {input.scaf_fai:q} \
            --output {output.pairs:q} --log {log:q} \
            --mapq {params.mapq} --outdir {params.outdir:q} \
            --sort-tmp {params.sort_tmp:q} --sort-mem {params.sort_mem:q} \
            --sort-threads {params.sort_threads} --threads {threads}
        '''

rule pretext_map:
    """Generate the Pretext contact map."""
    # Preserve scaffold order; MAPQ filtering is applied upstream.
    priority: 775
    input:
        pairs=CM_PREFIX + ".scaffolds.pairs.gz"
    output:
        pretext=CM_PREFIX + ".pretext"
    log:
        CM_PREFIX + ".pretext.log"
    threads: 2
    resources:
        mem_mb=4 * 1024, runtime=15, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        rm -f {output.pretext}.part
        trap 'rm -f {output.pretext}.part' EXIT
        zcat {input.pairs} | PretextMap -o {output.pretext}.part --sortby nosort > {log} 2>&1
        mv {output.pretext}.part {output.pretext}
        '''

rule pretext_snapshot:
    """Export a PNG preview of the Pretext map."""
    priority: 785
    input:
        pretext=CM_PREFIX + ".pretext"
    output:
        png=CM_PREFIX + ".pretext.png"
    params:
        outdir=CM_DIR,
        prefix="{sample}.{label}.hap{hap}.",
        resolution=1080
    log:
        CM_PREFIX + ".pretext.png.log"
    threads: 1
    resources:
        mem_mb=1024, runtime=10, slurm_partition=SLURM_PARTITION, slurm_account=SLURM_ACCOUNT
    shell:
        r'''
        set -euo pipefail
        PretextSnapshot -m {input.pretext} -f png -r {params.resolution} --sequences "=full" \
          -o {params.outdir} --prefix {params.prefix} > {log} 2>&1
        mv {params.outdir}/{params.prefix}FullMap.png {output.png}
        '''
