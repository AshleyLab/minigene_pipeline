import subprocess
import os
import yaml

### Determine the full path of which minigene_splicing_assay dir was copied to ###
minigene_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(minigene_dir, "minigene_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### List of scripts to run in order ###
scripts = [
    os.path.join(minigene_dir, "setup.py"),
    os.path.join(minigene_dir, "call_variants.py"),
    os.path.join(minigene_dir, "process_variants.py"),
    os.path.join(minigene_dir, "splice_analysis.py")
]

### Assign variables from minigene_config.yaml ###
root_output_dir = config["root_output_dir"]

### Create directory for logs ###
log_dir = os.path.join(root_output_dir, "minigene_pipeline_log")
os.makedirs(log_dir, exist_ok=True)

### Run each script sequentially ###
for idx, script in enumerate(scripts, start=1):
    log_file = os.path.join(log_dir, f"{os.path.basename(script)}.log")
    print(f"Running {script}... Output directed to {log_file}")

    ### Run the script and redirect stdout and stderr to the log file ###
    with open(log_file, "w") as log:
        result = subprocess.run(
            ["python3", script],
            stdout=log,
            stderr=subprocess.PIPE,  
            text=True
        )

        ### Write captured stderr to log file ###
        log.write(result.stderr)

    ### Verify successful completion of each script ###
    if result.returncode != 0:
        print(f"Error: {script} failed with return code {result.returncode}. See {log_file} for details.")

        ### Include a note about possible memory issues if return code matches known signals ###
        if result.returncode == -9:
            print(f"Note: Return code -9 usually indicates the process was killed, likely due to insufficient memory.")
            with open(log_file, "a") as log:
                log.write("\nNOTE: Return code -9 usually indicates the process was killed, likely due to insufficient memory.\n")

        ### Exit the pipeline if a script fails ###
        exit(1)
    else:
        print(f"{script} completed successfully.")

print("All scripts ran successfully.")
