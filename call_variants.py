### PYTHON SCRIPT 2 ###

import os
import yaml
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from datetime import datetime
import subprocess

### Determine the full path of which minigene_splicing_assay dir was copied to ###
minigene_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(minigene_dir, "minigene_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from minigene_config.yaml ###
root_output_dir = config["root_output_dir"]
reference_fasta = config["reference_fasta"]
deepvariant_sif = config["deepvariant_sif"]
model_type = config["model_type"]
max_workers = config["max_workers"]

### Define directories ###
base_output_dir = os.path.join(root_output_dir, "variant_barcode_results")
bam_files_dir = os.path.join(base_output_dir, "bam_files")
NUM_SHARDS = 1

def run_deepvariant(group_id):
    """Run DeepVariant for all BAM files belonging to a specific group."""
    bam_files = list(Path(bam_files_dir).glob(f"sorted_*_{group_id}.bam"))

    if not bam_files:
        print(f"No BAM files found for group {group_id}")
        return

    for bam_file in bam_files:
        barcode = bam_file.stem.split('_')[1]
        dv_output_dir = Path(base_output_dir) / "deepvariant" / f"{barcode}_{group_id}"
        dv_output_dir.mkdir(parents=True, exist_ok=True)

        output_vcf = dv_output_dir / f"{barcode}.vcf.gz"

        ### DeepVariant command ###
        command = [
            "singularity", "run",
            "-B", "/usr/lib/locale/:/usr/lib/locale/",
            "-B", f"{bam_files_dir}:{bam_files_dir}",
            "-B", f"{dv_output_dir}:{dv_output_dir}",
            "-B", f"{Path(reference_fasta).parent}:{Path(reference_fasta).parent}",
            deepvariant_sif,
            "/opt/deepvariant/bin/run_deepvariant",
            f"--model_type={model_type}",
            f"--ref={reference_fasta}",
            f"--reads={bam_file}",
            f"--output_vcf={output_vcf}",
            f"--num_shards={NUM_SHARDS}",
        ]

        print(f"Running DeepVariant for barcode {barcode} using file {bam_file}")

        ### Start timing ###
        job_start_time = datetime.now()

        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {bam_file}: {e}")
            continue

        ### End timing ###
        job_end_time = datetime.now()
        job_duration = (job_end_time - job_start_time).total_seconds()

        ### Log runtime ###
        runtime_log = dv_output_dir / f"runtime_{group_id}.log"
        with open(runtime_log, "a") as log_file:
            log_file.write(f"Runtime for task {group_id}: {job_duration} seconds\n")

        print(f"DeepVariant processing completed for barcode {barcode}")

def main():
    ### Identify unique group numbers from BAM file names ###
    bam_files = list(Path(bam_files_dir).glob("sorted_*.bam"))
    group_ids = {bam_file.stem.split('_')[-1] for bam_file in bam_files}

    if not group_ids:
        print("No groups found in BAM files.")
        return

    ### Run DeepVariant in parallel using ProcessPoolExecutor ###
    with ProcessPoolExecutor(max_workers=len(group_ids)) as executor:
        executor.map(run_deepvariant, group_ids)

if __name__ == "__main__":
    main()