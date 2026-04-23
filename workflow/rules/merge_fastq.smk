rule merge_fastq:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample],
        r2=lambda wc: SAMPLE_TO_R2[wc.sample]
    output:
        r1=SAMPLE_FASTQ_R1,
        r2=SAMPLE_FASTQ_R2,
        cmd=SAMPLE_MERGE_CMD
    log:
        SAMPLE_MERGE_LOG
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {log})

        # clean only target outputs
        rm -f {output.r1} {output.r2} {output.cmd}

        # save commands for reproducibility
        echo "zcat {input.r1} | gzip > {output.r1}" >  {output.cmd}
        echo "zcat {input.r2} | gzip > {output.r2}" >> {output.cmd}

        # merge (space-separated expansion handles multiple input files)
        zcat {input.r1} | gzip > {output.r1}
        zcat {input.r2} | gzip > {output.r2}
        """
