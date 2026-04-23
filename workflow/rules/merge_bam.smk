rule merge_bam:
    input:
        config=lambda wc: str(experiment_config_path(wc.exp)),
        bams=lambda wc: experiment_bam_inputs(wc.exp)
    output:
        bam=EXPERIMENT_OUTPUT_BAM
    params:
        outdir=EXPERIMENT_OUTPUT_DIR
    conda:
        "scDNA"
    log:
        EXPERIMENT_MERGE_BAM_LOG
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/rna_bam
        mkdir -p $(dirname {log})
        bash scripts/generate_bam.sh {input.config} 2>&1 | tee {log}
        """
