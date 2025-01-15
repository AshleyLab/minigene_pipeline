#!/bin/bash
#SBATCH --job-name=setup
#SBATCH --output=setup.%j.out
#SBATCH --error=setup.%j.err
#SBATCH --time=03:30:00
#SBATCH -p euan
#SBATCH -c 8
#SBATCH --mem=80GB

# Function to convert seconds to hours, minutes, and seconds
convert_seconds() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(( (total_seconds % 3600) / 60 ))
    local seconds=$((total_seconds % 60))
    printf "%02d:%02d:%02d\n" $hours $minutes $seconds
}

# Start time
start_time=$(date +%s)

# Check if the project directory is passed as an argument
if [ -z "$1" ]; then
    echo "Error: No project directory specified"
    exit 1
fi

# Set the project directory
minigene_dir="$1"

# Activate minigene_env
source "${minigene_dir}/minigene_env/bin/activate"

# Optional second argument for the final script to run (default is splice_analysis.py)
final_script="${2:-splice_analysis.py}"

# Optional third argument to disable running the final script (default is true)
run_final_script=${3:-true}

# Change to the project directory
cd "$minigene_dir" || exit 1

# Set the config file path dynamically
config_file="${minigene_dir}/minigene_config.yaml"

# Load configuration variables using the dynamically set config file
eval $(python3 "${minigene_dir}/export_config.py" -c "$config_file")

ROOT_OUTPUT_DIR="${ROOT_OUTPUT_DIR}"
BASE_OUTPUT_DIR="${ROOT_OUTPUT_DIR}/variant_barcode_results"
BAM_FILES_DIR="${BASE_OUTPUT_DIR}/bam_files"
FLAG_FILE="${BASE_OUTPUT_DIR}/initial_setup_done.flag"
REFERENCE_FASTA="${REFERENCE_FASTA}"
ARR_NUMBER="${ARR_NUMBER}"
MODEL_TYPE="${MODEL_TYPE}"

if [ ! -f "$FLAG_FILE" ]; then
    echo "Running setup.py"
    python3 "${minigene_dir}/setup.py"
    echo "Initial setup completed on $(date)" > "$FLAG_FILE"

    # Define the log directory
    log_dir="${minigene_dir}/log"

    # Submit the job array for all files and capture job ID
    array_job_id=$(sbatch --parsable --array=1-${ARR_NUMBER}%${ARR_NUMBER} -p euan --job-name="DeepVariant" --output="${log_dir}/DeepVariant.%j.out" --error="${log_dir}/DeepVariant.%j.err" --mem=3GB --time=07:00:00 -c 1 $0 "$minigene_dir")
    echo "Submitted DeepVariant processing array with Job ID: $array_job_id"

    # Submit the second Python script process_variants.py to run after the DeepVariant array jobs completes
    process_variants_job_id=$(sbatch --parsable --dependency=afterok:$array_job_id -p euan --job-name="process_variants" --output="${log_dir}/process_variants.%j.out" --error="${log_dir}/process_variants.%j.err" --mem=1GB --time=00:30:00 -c 4 --wrap="python3 ${minigene_dir}/process_variants.py ${minigene_dir}")
    echo "Submitted process_variants.py job with Job ID: $process_variants_job_id"

    if [ "$run_final_script" = true ]; then
        # Set resource parameters for the final script based on the script name
        if [ "$final_script" = "barcode_codon_rs.py" ]; then
            final_script_mem=1GB
            final_script_time=00:30:00
            final_script_c=1
        else
            final_script_mem=30GB
            final_script_time=07:00:00
            final_script_c=32
        fi

        # Set job name, output, and error file names based on the final script name
        final_script_name=$(basename "$final_script" .py)
        final_job_name="${final_script_name}"
        final_output="${log_dir}/${final_script_name}.%j.out"
        final_error="${log_dir}/${final_script_name}.%j.err"

        # Submit the user-specified final script with dependency on process_variants.py job
        sbatch --dependency=afterok:$process_variants_job_id -p euan --job-name="$final_job_name" --output="$final_output" --error="$final_error" --mem=$final_script_mem --time=$final_script_time -c $final_script_c --wrap="python3 ${minigene_dir}/${final_script} ${minigene_dir}"
    else
        echo "Skipping final script as requested by the user."
    fi

else
    echo "Running DeepVariant processing for task $SLURM_ARRAY_TASK_ID..."

    singularity_image="/oak/stanford/groups/euan/projects/variant_effect_mapping/tools/DeepVariant/deepvariant_1.6.1.sif"
    reference_fasta="${REFERENCE_FASTA}"
    model_type="${MODEL_TYPE}"
    num_shards=1

    for file in ${BAM_FILES_DIR}/sorted_*_${SLURM_ARRAY_TASK_ID}.bam; do
        barcode=$(basename "$file" .bam | cut -d'_' -f2) # Adjust based on your naming convention
        DV_OUTPUT_DIR="${BASE_OUTPUT_DIR}/deepvariant/${barcode}_${SLURM_ARRAY_TASK_ID}"
        mkdir -p "$DV_OUTPUT_DIR"

        # Construct output file paths
        output_vcf="${DV_OUTPUT_DIR}/${barcode}.vcf.gz"

        # Run DeepVariant command
        echo "Running DeepVariant for barcode ${barcode} using file ${file}"

        # Start timing the job
        JOB_START_TIME=$(date +%s)

        singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${BAM_FILES_DIR}":"${BAM_FILES_DIR}" \
        -B "${DV_OUTPUT_DIR}":"${DV_OUTPUT_DIR}" \
        -B "$(dirname "${reference_fasta}")":"$(dirname "${reference_fasta}")" \
        "${singularity_image}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=${model_type} \
        --ref="${reference_fasta}" \
        --reads="${file}" \
        --output_vcf="${output_vcf}" \
        --num_shards=${num_shards}

        # End timing the job
        JOB_END_TIME=$(date +%s)
        JOB_DURATION=$((JOB_END_TIME-JOB_START_TIME))
        echo "Runtime for task $SLURM_ARRAY_TASK_ID: $JOB_DURATION seconds" >> "${DV_OUTPUT_DIR}/runtime_${SLURM_ARRAY_TASK_ID}.log"

        echo "DeepVariant processing completed for barcode ${barcode}"
    done
fi

# End time
end_time=$(date +%s)

# Calculate the execution time
execution_time=$((end_time - start_time))

# Convert the execution time to hours, minutes, and seconds
formatted_execution_time=$(convert_seconds $execution_time)

echo "Execution Time: $formatted_execution_time"

# Deactivate the Python virtual environment
deactivate
