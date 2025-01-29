### PYTHON SCRIPT 1 ###

import os
import yaml
import re
import gzip
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, as_completed
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
from kneebow.rotor import Rotor
import time
import numpy as np
from scipy.optimize import curve_fit
import subprocess
from itertools import cycle
import shlex

### Determine the full path of which minigene_pipeline dir was copied to ###
minigene_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(minigene_dir, "minigene_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from minigene_config.yaml ###
fastq_path = config["fastq_path"]
DNA_baseq_threshold = config["DNA_baseq_threshold"]
reference_fasta = config["reference_fasta"]
gene = config["gene"]
root_output_dir = config["root_output_dir"]
fw_barcode_prefix = config['fw_barcode_prefix']
fw_barcode_suffix = config['fw_barcode_suffix']
rev_barcode_prefix = config['rev_barcode_prefix']
rev_barcode_suffix = config['rev_barcode_suffix']
barcode_length = config['barcode_length']
k = config['k']
arr_number = config["arr_number"]
max_workers = config["max_workers"]

mm2_max_workers = max_workers // 4

### Compile regex for finding 3' barcodes in both forward and reverse reads ###
fwd_bc_regex = re.compile(f"{fw_barcode_prefix}([ATCG]{{{barcode_length}}}){fw_barcode_suffix}")
rev_bc_regex = re.compile(f"{rev_barcode_prefix}([ATCG]{{{barcode_length}}}){rev_barcode_suffix}")
wt_fwd_regex = re.compile(f"{fw_barcode_prefix}{fw_barcode_suffix}")
wt_rev_regex = re.compile(f"{rev_barcode_prefix}{rev_barcode_suffix}")

def create_directories():
    """Create necessary directories if they don't exist"""
    if not os.path.exists(root_output_dir):
        os.makedirs(root_output_dir, exist_ok=True)

    global base_output_dir, chopper_minigenes_dir, barcode_rank_dir, demuxed_dir, sam_dir, bam_files_dir

    base_output_dir = os.path.join(root_output_dir, "variant_barcode_results")
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir, exist_ok=True)

    chopper_minigenes_dir = os.path.join(base_output_dir, "choppered_minigene")
    barcode_rank_dir =  os.path.join(base_output_dir, "barcode_rank")
    demuxed_dir = os.path.join(base_output_dir, "demuxed")
    sam_dir = os.path.join(base_output_dir, "sam_files")
    bam_files_dir = os.path.join(base_output_dir, "bam_files")
    for dir_path in [chopper_minigenes_dir, barcode_rank_dir, demuxed_dir, sam_dir, bam_files_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)

################################
### GENERATE REFERENCE INDEX ###
################################

def index_reference(reference_fasta):
    """Checks for *.fa.fai file; if not, generates the index using samtools faidx"""
    fai_file = f"{reference_fasta}.fai"
    if not os.path.exists(fai_file):
        print(f"Index file {fai_file} not found. Generating index...")
        try:
            subprocess.run(['samtools', 'faidx', reference_fasta], check=True)
            print(f"Index file {fai_file} created successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error generating index for {reference_fasta}: {e}")
    else:
        print(f"Index file {fai_file} already exists.")

####################################
### QUALITY CONTROL WITH CHOPPER ###
####################################

def run_chopper(fastq_path, chopper_minigenes_dir):
    """Run Chopper on the input cDNA FASTQ file for quality control"""
    input_fastq_file = shlex.quote(fastq_path)

    ### Generate the path for the Chopper-processed FASTQ file ###
    choppered_fastq_path = os.path.join(chopper_minigenes_dir, "choppered_minigenes.fastq")
    choppered_fastq = shlex.quote(choppered_fastq_path)

    ### Run the modified Chopper command using cat and piping ###
    chopper_cmd = f"cat {input_fastq_file} | chopper --quality {DNA_baseq_threshold} --threads {max_workers} > {choppered_fastq}"
    print(f"Running Chopper: {chopper_cmd}")
    subprocess.run(chopper_cmd, shell=True, check=True)

    return choppered_fastq_path


#########################################################################################                                                                  
### FIND ALL UNIQUE BARCODES AND FILTER OUT LOW-FIDELITY ONES USING BARCODE-RANK PLOT ###                                                                                  
#########################################################################################

### Function to compute the reverse complement of a DNA sequence ###
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def process_fastq(fastq_path, fwd_bc_regex, rev_bc_regex, wt_fwd_regex, wt_rev_regex):
    """
    Processes a fastq file to classify reads based on barcode and wildtype sequences i.e. demultiplexing

    Parameters:
    * fastq_path (str): Path to input fastq file
    * fwd_bc_regex (regex pattern): Compiled regex pattern to match forward barcodes in the forward 5' to 3' read
    * rev_bc_regex (regex pattern): Compiled regex pattern to match reverse barcodes in the reverse 3' to 5' read
    * wt_fwd_regex (regex pattern): Compiled regex pattern to match forward wildtype sequences in the forward read
    * wt_rev_regex (regex pattern): Compiled regex pattern to match reverse wildtype sequences in the reverse read

    Returns: read_storage (dict) containing classified reads
    * 'barcodes': A dictionary where keys are barcodes and values are lists of reads
    * 'wildtype': A list of wildtype reads
    * 'unclassified': A list of unclassified reads
    """
    read_storage = {'barcodes': {}, 'wildtype': [], 'unclassified': []}

    with open(fastq_path, 'r') as fastq_file:
        while True:
            id_line = fastq_file.readline().strip()
            if not id_line: break
            seq_line = fastq_file.readline().strip()
            plus_line = fastq_file.readline().strip()
            qual_line = fastq_file.readline().strip()

            fwd_match = fwd_bc_regex.search(seq_line)
            rev_match = rev_bc_regex.search(seq_line)
            wt_fwd_match = wt_fwd_regex.search(seq_line)
            wt_rev_match = wt_rev_regex.search(seq_line)

            if fwd_match:
                barcode = fwd_match.group(1)
                read_info = f"{id_line}\n{seq_line}\n{plus_line}\n{qual_line}"
                read_storage['barcodes'].setdefault(barcode, []).append(read_info)
            elif rev_match:
                barcode = reverse_complement(rev_match.group(1))
                seq_line = reverse_complement(seq_line)
                qual_line = qual_line[::-1]
                read_info = f"{id_line}\n{seq_line}\n{plus_line}\n{qual_line}"
                read_storage['barcodes'].setdefault(barcode, []).append(read_info)
            elif wt_fwd_match or wt_rev_match:
                if wt_rev_match:
                    seq_line = reverse_complement(seq_line)
                    qual_line = qual_line[::-1]
                read_info = f"{id_line}\n{seq_line}\n{plus_line}\n{qual_line}"
                read_storage['wildtype'].append(read_info)
            else:
                read_storage['unclassified'].append(f"{id_line}\n{seq_line}\n{plus_line}\n{qual_line}")

    return read_storage

def create_barcode_df(read_storage):
    """Creates a df with barcode counts from read_storage"""
    barcodes = [{'Barcode': barcode, 'Count': len(reads)} for barcode, reads in read_storage['barcodes'].items()]
    return pd.DataFrame(barcodes)

def process_barcode_counts(df_barcode_counts):
    """
    Processes barcode counts to rank barcodes and plot a barcode rank plot, 
    using kneebow to identify the kneepoint of plot i.e. threshold that separates high from low-quality barcodes

    Parameters:
    * df_barcode_counts (df): df containing barcode counts with 'Barcode' and 'Count' colums

    Workflow:
    * Rank barcodes based on their frequencies, with rank 1 given to barcode with highest count
    * Sort df by barcode rank in descending order
    * Plot Barcode Rank Plot of log(count of barcodes) vs barcode rank using Plotly
    * Use kneebow to automatically identify barcode knee point 
    * Barcode thresholds, valid barcodes, and the number of barcodes below the threshold stored as global variables for downstream processing
    """
    global barcode_rank_threshold
    global barcode_freq_threshold
    global valid_barcodes
    global num_barcodes_not_written

    df_barcode_counts['Barcode_Rank'] = df_barcode_counts['Count'].rank(method='dense', ascending=False).astype(int)
    df_barcode_counts.sort_values('Barcode_Rank', inplace=True)
    df_barcode = df_barcode_counts.reset_index(drop=True)

    ### Extract x and y values for the curve (log-transformed for better knee detection) ###
    x = np.log10(df_barcode_counts["Barcode_Rank"].values)
    y = np.log10(df_barcode_counts["Count"].values)

    ### Start and end points of the line ###
    x_start, y_start = x[0], y[0]
    x_end, y_end = x[-1], y[-1]

    ### Calculate the line equation coefficients (a, b, -c) for ax + by + c = 0 ###
    a = y_end - y_start
    b = -(x_end - x_start)
    c = (y_start * (x_end - x_start)) - (x_start * (y_end - y_start))

    ### Calculate the perpendicular distance from each point to the line ###
    distances = np.abs(a * x + b * y + c) / np.sqrt(a**2 + b**2)

    ### Find the index of the maximum distance ###
    knee_idx = np.argmax(distances)
    knee_x = 10 ** x[knee_idx]                                                          ### convert from log scale to original scale ###
    knee_y = 10 ** y[knee_idx]                                                          ### convert from log scale to original scale ###

    ### Get barcode rank and count threshold ###
    barcode_rank_threshold = int(round(k * knee_x))
    barcode_freq_threshold = df_barcode_counts.loc[df_barcode_counts["Barcode_Rank"] == barcode_rank_threshold, "Count"].values[0]

    ### Plot barcode rank plot with Plotly ###
    fig = go.Figure()

    ### Add main barcode rank plot as a scatter plot ###
    fig.add_trace(go.Scatter(x=df_barcode_counts["Barcode_Rank"] ,y=df_barcode_counts["Count"],
                             mode='markers', name="Barcode Rank Plot", marker=dict(color='darkblue', size=2)))

    ### Add barcode rank and count threshold as vertical line ###
    fig.add_trace(go.Scatter(x=[knee_x, knee_x], y=[df_barcode_counts["Count"].min()-10, df_barcode_counts["Count"].max()], 
                             mode='lines', name=f"Knee Point at Rank: {int(knee_x)}", line=dict(color='#800020', dash='dash')))

    ### Add knee point as marker dot ###
    fig.add_trace(go.Scatter(x=[knee_x], y=[knee_y], 
                             mode='markers', marker=dict(color='#800020', size=10), name="knee point"))

    ### Update plot layout ###
    fig.update_layout(
        title=f"Barcode Rank Plot of {gene} Minigene Plasmid Library",
        title_font=dict(size=30), 
        xaxis=dict(
            type="log",
            title="Barcode Rank (log)",
            title_font=dict(size=25),  
            tickfont=dict(size=20) 
        ),
        yaxis=dict(
            title="Barcode Frequency (log)",
            type="log",
            title_font=dict(size=25), 
            tickfont=dict(size=20)  
        ),
        plot_bgcolor='#fcf9ff',
        showlegend=False, 
        hoverlabel=dict(
            font_size=19 
        )
    )

    ### Add grid lines ###
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', gridcolor='lightgrey')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', gridcolor='lightgrey')

    ### Create path to save barcode_rank_plot.html ###
    html_file_path = os.path.join(barcode_rank_dir, "barcode_rank_plot.html")
    pio.write_html(fig, html_file_path)
    csv_file_path = os.path.join(barcode_rank_dir, "barcode_rank.csv")
    df_barcode.to_csv(csv_file_path, index=False)

    valid_barcodes = df_barcode_counts[df_barcode_counts['Barcode_Rank'] <= barcode_rank_threshold]['Barcode'].tolist()

    below_threshold_barcodes = df_barcode_counts[df_barcode_counts['Barcode_Rank'] > barcode_rank_threshold]
    num_barcodes_not_written = below_threshold_barcodes.shape[0]

def write_reads_to_file(reads, file_path):
    """Writes a list of read information to the specified FASTQ file"""
    with open(file_path, 'w') as f:
        for read in reads:
            f.write(f"{read}\n")

def write_reads_to_file_parallel(read_storage, valid_barcodes, output_dir):
    """Writes reads to their respective FASTQ files in parallel, including wildtype and unclassified reads."""
    def task(writer_function, reads, file_path):
        writer_function(reads, file_path)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        ### Write high-quality barcode reads as determined by knee plot ###
        for barcode, reads in read_storage['barcodes'].items():
            if barcode in valid_barcodes:
                output_file_path = os.path.join(output_dir, f"{barcode}.fq")
                executor.submit(task, write_reads_to_file, reads, output_file_path)

        ### Write wildtype and unclassified reads if they exist ###
        if read_storage['wildtype']:
            wt_file_path = os.path.join(output_dir, "wildtype.fq")
            executor.submit(task, write_reads_to_file, read_storage['wildtype'], wt_file_path)

        if read_storage['unclassified']:
            unclassified_file_path = os.path.join(output_dir, "unclassified.fq")
            executor.submit(task, write_reads_to_file, read_storage['unclassified'], unclassified_file_path)

######################################################################################                                                             
### ALIGNMENT WITH MINIMAP2 AND CONVERT OUTPUT TO SORTED BAM FILES FOR DEEPVARIANT ###                                                                                  
######################################################################################

def run_minimap2(fq_file, output_dir, reference_fasta):
    """Runs minimap2 for a given FASTQ file and saves the SAM output"""
    sam_file_name = os.path.basename(fq_file).replace('.fq', '.sam')
    output_file_path = os.path.join(output_dir, sam_file_name)
    minimap2_cmd = f"minimap2 -ax map-ont --cs -t 4 {reference_fasta} {fq_file} > {output_file_path}"
    subprocess.run(minimap2_cmd, shell=True, check=True)

def process_fastq_files_with_minimap2(input_dir, output_dir, reference_fasta):
    """Finds all FASTQ files in the input directory, except 'wildtype.fq' and 'unclassified.fq', and processes them with minimap2 in parallel; 
       no need for wildtype.fq to be processed with mm2 since it does not contain barcodes and no variant-barcode relationship can be drawn"""
    excluded_files = {'wildtype.fq', 'unclassified.fq'}
    fastq_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fq') and f not in excluded_files]
    with ProcessPoolExecutor(max_workers=mm2_max_workers) as executor:
        for fq_file in fastq_files:
            executor.submit(run_minimap2, fq_file, output_dir, reference_fasta)

def convert_sort_index_sam(sam_file, bam_files_dir, group_number):
    """Converts a SAM file to a sorted and indexed BAM file, including group number in the filename"""
    basename = os.path.basename(sam_file).replace('.sam', '')
    bam_file_path = os.path.join(bam_files_dir, f"{basename}_{group_number}.bam")
    sorted_bam_file_path = os.path.join(bam_files_dir, f"sorted_{basename}_{group_number}.bam")
    subprocess.run(['samtools', 'view', '-bS', sam_file, '-o', bam_file_path], check=True)
    subprocess.run(['samtools', 'sort', bam_file_path, '-o', sorted_bam_file_path], check=True)
    subprocess.run(['samtools', 'index', sorted_bam_file_path], check=True)

def distribute_and_process_sam_files(input_dir, output_dir):
    sam_files = [f for f in os.listdir(input_dir) if f.endswith('.sam')]
    total_files = len(sam_files)

    ### generates a repeating sequence from 1 to arr_number specified in minigene_config.yaml ###
    group_numbers = [i % arr_number + 1 for i in range(total_files)]
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for sam_file, group_number in zip(sam_files, group_numbers):
            full_sam_path = os.path.join(input_dir, sam_file)
            future = executor.submit(convert_sort_index_sam, full_sam_path, output_dir, group_number)
            futures.append(future)

        ### wait for all submitted jobs to complete ###
        for future in as_completed(futures):
            try:
                ### each future.result() call re-raises any exception caught during execution ###
                future.result()
            except Exception as exc:
                print(f'Generated an exception: {exc}')

    print(f"Processed {total_files} files across {max(group_numbers)} groups in parallel.")

def main():
    create_directories()

    index_reference(reference_fasta)

    chopper_start_time = time.time()
    choppered_fastq_path = run_chopper(fastq_path, chopper_minigenes_dir)
    chopper_end_time = time.time()
    print(f"Time taken to QC fastq with chopper: {chopper_end_time - chopper_start_time:.3f} seconds")

    id_bc_start_time = time.time()
    read_storage = process_fastq(choppered_fastq_path, fwd_bc_regex, rev_bc_regex, wt_fwd_regex, wt_rev_regex)
    id_bc_end_time = time.time()
    print(f"Time taken to ID unique barcodes from fastq: {id_bc_end_time - id_bc_start_time:.3f} seconds")

    filter_bc_start_time = time.time()
    df_barcode_counts = create_barcode_df(read_storage)
    process_barcode_counts(df_barcode_counts)
    filter_bc_end_time = time.time()
    print(f"Time taken to assess barcodes quality and filter: {filter_bc_end_time - filter_bc_start_time:.3f} seconds")

    demux_start_time = time.time()
    write_reads_to_file_parallel(read_storage, valid_barcodes, demuxed_dir)
    print(f"{num_barcodes_not_written} barcodes discarded due to barcode rank < barcode rank threshold of {barcode_rank_threshold} and barcode frequency < {barcode_freq_threshold}")
    print("FASTQ file generation complete.")
    demux_end_time = time.time()
    print(f"Time taken to demultiplex fastq by filtered barcodes: {demux_end_time - demux_start_time:.3f} seconds")

    mm2_aln_start_time = time.time()
    process_fastq_files_with_minimap2(demuxed_dir, sam_dir, reference_fasta)
    mm2_aln_end_time = time.time()
    print(f"Time taken to align fastq with minimap2: {mm2_aln_end_time - mm2_aln_start_time:.3f} seconds")

    generate_bam_start_time = time.time()
    distribute_and_process_sam_files(sam_dir, bam_files_dir)
    print("SAM to BAM conversion, sorting, and indexing complete.")
    generate_bam_end_time = time.time()
    print(f"Time taken to convert SAM to sorted BAM files: {generate_bam_end_time - generate_bam_start_time:.3f} seconds")

if __name__ == "__main__":
    main()


