rule run_paired_tag:
    input:
        r1=lambda wc: sample_fastq_r1_path(wc.sample),
        r2=lambda wc: sample_fastq_r2_path(wc.sample)
    output:
        bam=SAMPLE_USEFUL_BAM
    params:
        library=lambda wc: get_sample_meta(wc.sample, "library"),
        species=lambda wc: get_sample_meta(wc.sample, "species"),
        seq_type=lambda wc: get_sample_meta(wc.sample, "seq_type"),
        workdir=lambda wc: sample_run_dir(wc.sample),
        parent_dir=lambda wc: sample_merge_dir(wc.sample)
    conda:
        "paired_tag"
    log:
        SAMPLE_BAM_LOG
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.parent_dir}
        mkdir -p $(dirname {log})
        cd {params.workdir}

        bash ~/github/Paired-Tag-Pipeline/pipelines/pipeline_PT.sh \
                -g {params.species} -l {params.library} -p 4 \
                2>&1 | tee {log}

        # Remove large intermediate files
        rm -f bw/*.bw
        rm -f process/*.bam
        rm -f process/*.gz
        rm -f process/*.fa
        rm -f process/*ReadType.txt
        rm -f bam/*Assign2BC.bam
        rm -f bam/*Barcode.bam
        """
