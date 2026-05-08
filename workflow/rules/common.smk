from pathlib import Path
import re

import pandas as pd
from snakemake.exceptions import WorkflowError


WORKFLOW_DIR = Path(workflow.basedir).resolve()
PROJECT_ROOT = WORKFLOW_DIR.parent

PIPELINE_BASE = config.get(
    "pipeline_base",
    "/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/data/xx_pipeline",
)
MERGE_BASE = config.get(
    "merge_base",
    "/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/data/merge",
)


def project_path(path_str):
    path = Path(path_str)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def sample_fastq_r1_path(sample):
    return f"{PIPELINE_BASE}/{sample}/merge/fastq/{sample}_R1.fastq.gz"


def sample_fastq_r2_path(sample):
    return f"{PIPELINE_BASE}/{sample}/merge/fastq/{sample}_R2.fastq.gz"


def sample_merge_cmd_path(sample):
    return f"{PIPELINE_BASE}/{sample}/merge/fastq/merge_cmd.txt"


def sample_merge_dir(sample):
    return f"{PIPELINE_BASE}/{sample}/merge"


def sample_run_dir(sample):
    return f"{sample_merge_dir(sample)}"


def sample_bam_path(sample):
    return f"{sample_run_dir(sample)}/{sample}/bam/{sample}_UsefulReads.bam"


def experiment_output_bam_path(exp):
    return f"{MERGE_BASE}/{exp}/rna_bam/RNA_final.bam"


def experiment_output_dir(exp):
    return f"{MERGE_BASE}/{exp}"


def experiment_config_path(exp):
    return project_path(f"config/{exp}.sh")


def log_path(filename):
    return str(project_path(f"logs/{filename}"))


samples = pd.read_csv(project_path(config["samples"]), sep="\t").fillna("")
experiments = pd.read_csv(project_path(config["experiments"]), encoding="utf-8-sig")

SAMPLE_TO_R1 = samples.groupby("sample")["R1"].apply(list).to_dict()
SAMPLE_TO_R2 = samples.groupby("sample")["R2"].apply(list).to_dict()
SAMPLE_TO_META = samples.groupby("sample").first().to_dict("index")
SAMPLES = sorted(SAMPLE_TO_R1)
EXPS = [
    str(exp).strip()
    for exp in experiments["proj"].dropna().tolist()
    if str(exp).strip()
]

EXP_TO_SAMPLES = {}


def parse_experiment_samples(config_script):
    config_text = Path(config_script).read_text()
    match = re.search(r"dna_rna=\((.*?)\)", config_text, re.S)
    if not match:
        return []

    pairs = re.findall(r'"([^"]+)"', match.group(1))
    samples = []
    for pair in pairs:
        samples.extend(part for part in pair.split("_") if part)
    return sorted(set(samples))


def get_sample_meta(sample, field):
    return SAMPLE_TO_META[sample][field]


def require_ready_experiment(exp):
    config_script = experiment_config_path(exp)
    if not config_script.exists():
        raise WorkflowError(f"Missing config script for experiment '{exp}': {config_script}")

    if exp not in EXP_TO_SAMPLES:
        EXP_TO_SAMPLES[exp] = parse_experiment_samples(config_script)

    missing_samples = [sample for sample in EXP_TO_SAMPLES[exp] if sample not in SAMPLE_TO_R1]
    if missing_samples:
        missing_text = ", ".join(missing_samples)
        raise WorkflowError(
            f"Experiment '{exp}' has samples missing from {config['samples']}: {missing_text}"
        )


def experiment_bam_inputs(exp):
    require_ready_experiment(exp)
    return [sample_bam_path(sample) for sample in EXP_TO_SAMPLES[exp]]


SAMPLE_FASTQ_R1 = sample_fastq_r1_path("{sample}")
SAMPLE_FASTQ_R2 = sample_fastq_r2_path("{sample}")
SAMPLE_MERGE_CMD = sample_merge_cmd_path("{sample}")
SAMPLE_USEFUL_BAM = sample_bam_path("{sample}")
EXPERIMENT_OUTPUT_BAM = experiment_output_bam_path("{exp}")
EXPERIMENT_OUTPUT_DIR = experiment_output_dir("{exp}")
SAMPLE_MERGE_LOG = log_path("{sample}_merge.log")
SAMPLE_BAM_LOG = log_path("{sample}_bam.log")
EXPERIMENT_MERGE_BAM_LOG = log_path("{exp}_merge_bam.log")
