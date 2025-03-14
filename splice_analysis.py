### PYTHON SCRIPT 3 ###

import os
import yaml
import re
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, as_completed
from collections import defaultdict
import pandas as pd
import glob
import subprocess
import shutil
from tempfile import NamedTemporaryFile
import time
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
import warnings
import shlex
import logging
import multiprocessing
import threading
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots

### Determine the full path of which minigene_splicing_assay dir was copied to ###
minigene_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(minigene_dir, "minigene_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from minigene_config.yaml ###
root_output_dir = config["root_output_dir"]
cDNA_fastq_file = config["cDNA_fastq_file"]
cDNA_baseq_threshold = config["cDNA_baseq_threshold"]
reference_fasta = config["reference_fasta"]
index_sequence_patterns = config["index_sequence_patterns"]
fw_barcode_prefix = config["fw_barcode_prefix"]
fw_barcode_suffix = config["fw_barcode_suffix"]
database_name = config["database_name"]
ref_start_pos = config["ref_start_pos"]
ref_end_pos = config["ref_end_pos"]
breakpoint_freq_threshold = config["breakpoint_freq_threshold"]
very_close_bp = config["very_close_bp"]
breakpoint_padding = config["breakpoint_padding"]
max_workers = config["max_workers"]

### Files needed for second demux ###
barcode_txt_file_path = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "unique_SNV_barcodes.txt")
variant_info_csv =  os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "SNV_barcode_info.csv")

### Path needed for GMAP Alignment ###
gmap_path = os.path.join(minigene_dir, "minigene_env", "bin")

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

def create_directories():
    """Create necessary directories if they don't exist"""
    global base_result_directory
    base_result_directory = os.path.join(root_output_dir, "splice_results")
    if not os.path.exists(base_result_directory):
        os.makedirs(base_result_directory, exist_ok=True)
    global visoqlr_base, breakpoints_base
    visoqlr_base = os.path.join(base_result_directory, "VIsoQLR")
    breakpoints_base = os.path.join(base_result_directory, "Breakpoints_data")
    for base_path in [visoqlr_base, breakpoints_base]:
        if not os.path.exists(base_path):
            os.makedirs(base_path, exist_ok=True)

    ### make these directories accessible globally ###
    global chopper_cDNA_dir, demux_index_dir, sample_barcode_demux_dir, gmap_dir, visoqlr_dir, visoqlr_no_cdna_dir, visoqlr_combined_bc_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir, sample_gff3_dir, sample_breakpoint_freq_plot_dir
    chopper_cDNA_dir = os.path.join(base_result_directory, "cDNA_QC_fastq")
    demux_index_dir = os.path.join(base_result_directory, "Sample_Fastq")
    sample_barcode_demux_dir = os.path.join(base_result_directory, "Sample_Barcode_Demuxed_Fastq")
    gmap_dir = os.path.join(base_result_directory, "GMAP")
    visoqlr_dir = os.path.join(visoqlr_base, "VIsoQLR_with_cDNA_csv")
    visoqlr_no_cdna_dir = os.path.join(visoqlr_base, "VIsoQLR_no_cDNA_csv")
    visoqlr_combined_bc_dir = os.path.join(visoqlr_base, "VIsoQLR_per_sample_variant")
    VIsoQLR_plots_dir = os.path.join(visoqlr_base, "VIsoQLR_plots")
    sample_variant_gff3_dir = os.path.join(breakpoints_base, "breakpoints_per_sample_variant_csv")
    sample_gff3_dir = os.path.join(breakpoints_base, "breakpoints_per_sample_csv")
    sample_breakpoint_freq_plot_dir = os.path.join(breakpoints_base, "breakpoint_frequency_plots_per_sample")
    for dir_path in [chopper_cDNA_dir, demux_index_dir, sample_barcode_demux_dir, gmap_dir, visoqlr_dir, visoqlr_no_cdna_dir, visoqlr_combined_bc_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir, sample_gff3_dir, sample_breakpoint_freq_plot_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)

##################################
### PREPROCESSING WITH CHOPPER ###
##################################

def run_chopper(cDNA_fastq_file, chopper_cDNA_dir):
    """Run Chopper on the input cDNA FASTQ file for quality control"""
    input_fastq_file = shlex.quote(cDNA_fastq_file)

    ### Generate the path for the Chopper-processed FASTQ file ###
    choppered_fastq_path = os.path.join(chopper_cDNA_dir, "choppered_cDNA.fastq")
    choppered_fastq = shlex.quote(choppered_fastq_path)

    ### Run the modified Chopper command using cat and piping ###
    chopper_cmd = f"cat {input_fastq_file} | chopper --quality {cDNA_baseq_threshold} --threads {max_workers} > {choppered_fastq}"
    print(f"Running Chopper: {chopper_cmd}")
    subprocess.run(chopper_cmd, shell=True, check=True)

    return choppered_fastq_path


############################################################################################                                                                               
### FIRST DEMULTIPLEXING AND TRIMMING OF INDEXES/ADAPTERS TO OBTAIN {SAMPLE}.FASTQ FILES ###                                                                                  
############################################################################################

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])

def find_and_trim_sequence(sequence, quality, fwd_front, fwd_end, rev_front, rev_end):
    """
    Searches for specified forward and reverse barcode sequences within each fastq read, extracting region
    between them if the sequences are found, along with its seqID and corresponding quality scores

    Parameters:
    * sequence (str): The DNA sequence to be searched
    * quality (str): The quality scores corresponding to the DNA sequence
    * fwd_front (str): The forward barcode's start sequence
    * fwd_end (str): The forward barcode's end sequence
    * rev_front (str): The reverse barcode's start sequence
    * rev_end (str): The reverse barcode's end sequence

    Returns:
    * extracted_seq (str), extracted_quality (str): The extracted DNA sequence and quality scores if a match is found; otherwise, returns (None, None)

    Workflow:
    1. Contains inner is_close_match function that: 
        * checks if a subsequence matches a given pattern with up to a specified number of allowed mismatches
        * returns True if the number of mismatches is within the allowed limit, False otherwise

    2. Processing Forward Matches:
        * searches for `fwd_front` sequence starting from the beginning of `sequence`, storing index as `fwd_start_index` if match is found
        * searches for `fwd_end` sequence starting from the end of `sequence`, moving backwards and storing index as `fwd_end_index` if a match is found
        * if both `fwd_front` and `fwd_end` are found and `fwd_start_index` < `fwd_end_index`, the function extracts the subsequence and corresponding quality scores between these indices
        * returns the extracted sequence and quality scores

    3. Processing Reverse Match:
        * searches for `rev_front` sequence starting from the beginning of `sequence` if 'fwd_front' is not found, storing index as `rev_start_index` if 'rev_front' match is found
        * searches for `rev_end` sequence starting from the end of `sequence', moving backwards and stores index as `rev_end_index` if `rev_end` match is found
        * if both `rev_front` and `rev_end` are found and `rev_start_index` is < `rev_end_index`, the function extracts the subsequence and corresponding quality scores between these indices
        * the extracted sequence is reverse complemented, and the quality scores are reversed before returning the reverse complemented sequence and reversed quality scores

    4. Fallback:
        * if neither a forward nor a reverse match is found, the function returns `(None, None)` to indicate that no valid match was detected
    """
    def is_close_match(seq, pattern, allowed_mismatches=1):
        mismatches = sum(1 for a, b in zip(seq, pattern) if a != b)
        return mismatches <= allowed_mismatches

    ### Process forward match ###
    fwd_start_index = -1
    fwd_end_index = -1

    ### Look for fwd_front from start of sequence ###
    for i in range(len(sequence) - len(fwd_front) + 1):
        if is_close_match(sequence[i:i+len(fwd_front)], fwd_front):
            fwd_start_index = i
            break

    ### Look for fwd_end starting from the end of sequence and moving backward ###
    for i in range(len(sequence) - len(fwd_end), -1, -1):
        if is_close_match(sequence[i:i+len(fwd_end)], fwd_end):
            fwd_end_index = i + len(fwd_end)
            break

    if fwd_start_index != -1 and fwd_end_index != -1 and fwd_start_index < fwd_end_index:
        extracted_seq = sequence[fwd_start_index:fwd_end_index]
        extracted_quality = quality[fwd_start_index:fwd_end_index]
        return extracted_seq, extracted_quality

    ### Process reverse match ###
    rev_start_index = -1
    rev_end_index = -1

    ### Look for rev_front from start of sequence ###
    for i in range(len(sequence) - len(rev_front) + 1):
        if is_close_match(sequence[i:i+len(rev_front)], rev_front):
            rev_start_index = i
            break

    ### Look for rev_end starting from the end of sequence and moving backward ###
    for i in range(len(sequence) - len(rev_end), -1, -1):
        if is_close_match(sequence[i:i+len(rev_end)], rev_end):
            rev_end_index = i + len(rev_end)
            break

    if rev_start_index != -1 and rev_end_index != -1 and rev_start_index < rev_end_index:
        extracted_seq = sequence[rev_start_index:rev_end_index]
        extracted_quality = quality[rev_start_index:rev_end_index]
        ### Reverse the extracted sequence and quality ###
        return reverse_complement(extracted_seq), extracted_quality[::-1]

    return None, None

def process_fastq_chunk(data_chunk, index_sequence_patterns):
    """
    Processes a chunk of FASTQ data to identify sequences that match specified index patterns,
    trimming and organising matched sequences into separate outputs for each sample, based on their index.

    Parameters:
    * data_chunk (list of str): A list of strings containing a chunk of FASTQ data
    * index_sequence_patterns (dict): A dictionary containing the index patterns for each sample where keys are sample names, and the values are dictionaries with the following keys: 'fwd_front', 'fwd_end', 'rev_front' & 'rev_end'

    Returns:
    * results (dict): A dictionary where the keys are sample names and the values are lists of matched and trimmed FASTQ entries for each sample

    Workflow:
    1. Initialize Output Dictionary:
        * create a dictionary `results` with keys corresponding to sample names from `index_sequence_patterns`, and initalise values as empty lists

    2. Process Each FASTQ Record:
        * convert the `data_chunk` list into an iterator `it` for sequential processing
        * enter a loop to process each FASTQ record within chunk

    3. Barcode Matching:
        * attempt to find a match using the `find_and_trim_sequence` function for each sample in `index_sequence_patterns`
        * use `find_and_trim_sequence` function to search for specified index sequences in the `sequence`, extract the sequence and corresponding quality scores
       
    4. Validation and Trimming:
        * validate that both the extracted sequence and quality scores are non-empty after finding a match
        * further trim 8 bases from the front and end of both `trimmed_seq` and `trimmed_quality` to obtain `final_trimmed_seq` and `final_trimmed_quality`
        * ensure that the final trimmed sequence and quality are still non-empty after trimming

    5. Store Results:
        * append FASTQ record (with `seqID`, `final_trimmed_seq`, `separator`, `final_trimmed_quality`) to list for sample in `results` dictionary if the final trimmed sequence and quality are valid

    6. Return Results:
        * return the `results` dictionary, which contains the demultiplexed and trimmed FASTQ entries for each sample after processing all records in the chunk

    Exception Handling:
    * If the loop reaches the end of the data chunk, a `StopIteration` exception is caught, and the loop exits
    """
    results = {sample: [] for sample in index_sequence_patterns}
    it = iter(data_chunk)
    while True:
        try:
            identifier = next(it).strip()                                                                                       ### read the seqID line ###
            if not identifier:                                                                                                  ### break if empty ###
                break
            sequence = next(it).strip()                                                                                         ### read sequence line ###
            separator = next(it).strip()                                                                                        ### read separator line '+' ###
            quality = next(it).strip()                                                                                          ### read quality line ###
            for sample, patterns in index_sequence_patterns.items():
                trimmed_seq, trimmed_quality = find_and_trim_sequence(
                    sequence, quality,
                    patterns['fwd_front'], patterns['fwd_end'],
                    patterns['rev_front'], patterns['rev_end']
                )
                ### check that both sequence and quality are non-empty ###
                if trimmed_seq and trimmed_quality and len(trimmed_seq) > 0 and len(trimmed_quality) > 0:
                    final_trimmed_seq = trimmed_seq[8:-8]
                    final_trimmed_quality = trimmed_quality[8:-8]
                    if len(final_trimmed_seq) > 0 and len(final_trimmed_quality) > 0:
                        results[sample].append(f"{identifier}\n{final_trimmed_seq}\n{separator}\n{final_trimmed_quality}\n")
        except StopIteration:
            break
    return results

def demux_index(choppered_fastq, demux_index_dir, index_sequence_patterns):
    """
    Demultiplexes a FASTQ file containing cDNA sequences based on specified index sequence patterns, processing the file in chunks, 
    and uses parallel processing to improve efficiency, before finally writing the resulting sequences to separate FASTQ files for each sample i.e. demultiplexing

    Parameters:
    * cDNA_fastq_file (str): Path to the input FASTQ file containing cDNA sequences
    * demux_index_dir (str): Directory where the demultiplexed {sample}.fastq files will be saved
    * index_sequence_patterns (dict): A dictionary where keys are sample names and values are dictionaries containing the forward and reverse sequence patterns for demultiplexing; each value dictionary has the keys: fwd_front, fwd_end, rev_front, rev_end

    Workflow:
    * Reads all lines in fastq file into 'data' list
    * Splits the data list into smaller chunks of size chunk_size, where each chunk is processed independently
    * Creates a temporary file for each sample to store the demultiplexed sequences before they are written to the final output files
    * Uses ProcessPoolExecutor to process each chunk in parallel, where process_fastq_chunk function is submitted as a task for each chunk, and results are written to corresponding sample's temp file as each task completes
    * Close and merge temporary files
        - closes each temp file
        - opens the final output file for each sample in append mode ('ab')
        - reads from the temp file and writes its contents to the final output file
        - deletes the temp file after merging its contents
    """
    with open(choppered_fastq, 'r') as file:
        data = file.readlines()

    chunk_size = 100000
    chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]                                                  ### create chunks of data to distribute to workers ###

    temp_files = {sample: NamedTemporaryFile(mode='w+', delete=False) for sample in index_sequence_patterns}                    ### prepare temporary files for each sample ###

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:                                           ### use a ProcessPoolExecutor to process each chunk in parallel ###
        futures = [executor.submit(process_fastq_chunk, chunk, index_sequence_patterns) for chunk in chunks]

        for future in concurrent.futures.as_completed(futures):                                                                 ### collect results as they complete ###
            results = future.result()
            for sample, lines in results.items():
                temp_files[sample].write(''.join(lines))

    for sample, temp_file in temp_files.items():                                                                                ### close and merge temporary files into final outputs ###
        temp_file.close()
        final_path = os.path.join(demux_index_dir, f'{sample}.fastq')
        with open(final_path, 'ab') as final_file:
            with open(temp_file.name, 'rb') as tf:
                shutil.copyfileobj(tf, final_file)
        os.unlink(temp_file.name)                                                                                               ### delete temp files ###

###############################################################################################                                                                         
### SECOND DEMULTIPLEXING BY KNOWN BARCODES TO OBTAIN {BARCODE}.FASTQ FILES FOR EACH SAMPLE ### 
###############################################################################################

### Define constants ###
LINES_PER_RECORD = 4
CHUNK_SIZE = 25000 * LINES_PER_RECORD

### Define regex pattern for exact matching to wildtype reads ### 
wt_regex = re.compile(f"{fw_barcode_prefix}{fw_barcode_suffix}")

def read_barcodes(filename):
    """
    Reads unique_SNV_barcodes.txt file containing barcodes and generates a dictionary of compiled regular expression patterns based on these barcodes, 
    where each pattern is constructed by embedding the barcode between pre-defined preffix and suffix

    Parameters:
    * filename (str): Path to barcode.txt file, where each barcode is written on a separate line

    Returns:
    * barcode_regexes (dict): A dictionary where keys are barcodes (str) and values are compiled regex patterns (re.compile objects)
    """
    barcode_regexes = {}
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            print(f"Total lines read: {len(lines)}")
            for line in lines:
                barcode = line.strip()
                if barcode:
                    pattern = f"{fw_barcode_prefix}{barcode}{fw_barcode_suffix}"
                    barcode_regexes[barcode] = re.compile(pattern)
                else:
                    print("Empty or whitespace-only line detected.")
    except Exception as e:
        print(f"Error reading barcode file: {e}")
    return barcode_regexes

def process_sample_fastq_chunk(data, barcodes):
    """
    Processes a chunk of FASTQ data to identify sequences that match specified barcode patterns,
    organising the matched sequences into separate FASTQ files for each barcode and wildtype sequences ({barcode}.fastq, wildtype.fastq)

    Parameters:
    * data (str): A string containing FASTQ data, where each set of four lines represents a single sequence entry (seqID, sequence, separator, quality)
    * barcodes (dict): A dictionary where keys are barcodes (str) and values are compiled regex patterns (re.compile objects) to search within the sequences

    Returns:
    * outputs (defaultdict): A dictionary where keys are filenames ("{barcode}.fastq" or "wildtype.fastq") and values are lists of matched FASTQ entries

    Workflow:
    * Creates a defaultdict of lists to store the classified FASTQ records where each key corresponds to a filename, and each value is a list of FASTQ records
    * Splits the input FASTQ data string into individual lines, stored in the lines list
    * Iterates over the lines in chunks of four, extracting the seqID, sequence, '+', and quality scores for each read
    * Match sequences to barcode and wildtype patterns:
        - initialises a matched flag to False
        - iterates over each barcode and its corresponding regex pattern in the barcodes dictionary
        - checks if the sequence matches the barcode regex using regex.search(sequence)
        - appends the FASTQ read to the corresponding output list (keyed by {barcode}.fastq), sets the matched flag to True, and breaks out of loop if barcode match is found
        - checks if the sequence matches the wildtype pattern using wt_regex.search(sequence) if no barcode match was found (matched is False)
        - appends the FASTQ read to the wildtype.fastq output list if wildtype match is found
    * Returns the outputs dictionary containing demultipleded FASTQ reads
    """
    outputs = defaultdict(list)
    lines = data.split('\n')
    for i in range(0, len(lines) - LINES_PER_RECORD + 1, LINES_PER_RECORD):
        header = lines[i].strip()
        sequence = lines[i+1].strip()
        plus = lines[i+2].strip()
        quality = lines[i+3].strip()

        matched = False
        for barcode, regex in barcodes.items():
            if regex.search(sequence):
                outputs[f"{barcode}.fastq"].append(f"{header}\n{sequence}\n{plus}\n{quality}\n")
                matched = True
                break

        if not matched and wt_regex.search(sequence):
            outputs["wildtype.fastq"].append(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    return outputs

def process_sample_fastq_file(fastq_file, barcodes_file, output_dir):
    """
    Processes each {sample}.fastq file in chunks to classify sequences based on pre-defined barcode patterns before writing the results to separate FASTQ files in the specified output directory

    Parameters:
    * fastq_file (str): Path to input {sample}.fastq file containing sequence data
    * barcodes_file (str): Path to unique_SNV_barcodes.txt file containing barcodes that are each written on a separate line
    * output_dir (str): Directory where the demultiplexed {sample}_{barcode}.fastq files will be saved

    Workflow:
    * Calls read_barcodes function to read and compile regex patterns from barcode.txt file and store the resulting dictionary of barcode patterns in "barcodes" variable
    * Creates a defaultdict of lists to store the classified FASTQ reads, where each key corresponds to a filename, and each value is a list of FASTQ reads
    * Read and process FASTQ file in chunks:
        - opens the specified fastq_file in read mode  and reads the file in chunks of CHUNK_SIZE lines, storing each chunk in "data" variable
        - breaks out of the loop if data is empty (end of file)
        - calls the process_sample_fastq_chunk function to process the chunk of data using the barcodes dictionary
        - appends the resulting classified FASTQ reads to the corresponding lists in the final_outputs dictionary
    * Write demultiplexed data to output files:
        - iterates over the final_outputs dictionary
        - constructs the output file path for each filename
        - opens the output file in append mode and writes the demuxed FASTQ reads
        - catches any exceptions that occur during the file reading and processing, printing an error message indicating the nature of the exception and the file being processed
    """
    barcodes = read_barcodes(barcodes_file)

    try:
        with open(fastq_file, 'r') as file:
            while True:
                data = ''.join([file.readline() for _ in range(CHUNK_SIZE)])
                if not data:
                    break
                chunk_outputs = process_sample_fastq_chunk(data, barcodes)
                for filename, content in chunk_outputs.items():
                    output_filepath = os.path.join(output_dir, filename)
                    with open(output_filepath, 'a') as f:
                        f.writelines(content)
        print(f"Data written for {fastq_file}")
    except Exception as e:
        print(f"Error processing file {fastq_file}: {e}")

def create_barcode_variant_mapping(variant_info_csv):
    """
    Reads variant information from variants_info.csv and creates a mapping from barcodes to variants

    Parameters:
    * variant_info_csv (str): Path to variants_info.csv containing variant information, with columns named 'BARCODE' and 'VARIANT'

    Returns:
    * barcode_to_variant (dict): A dictionary where keys are barcodes (str) and values are variants (str)

    Workflow:
    * Reads variants_info.csv into a df
    * Initialises an empty dictionary for mapping barcodes to variants
    * Iterates over each row of the df and splits the 'BARCODE' column into individual barcodes
    * Maps each barcode to the corresponding 'VARIANT' value in the dictionary
    * Returns the barcode-to-variant mapping dictionary
    """
    df = pd.read_csv(variant_info_csv)
    barcode_to_variant = {}
    for _, row in df.iterrows():
        barcodes = row['BARCODE'].split(',')
        for barcode in barcodes:
            barcode_to_variant[barcode] = row['VARIANT']
    return barcode_to_variant

def organise_files_by_variant(output_dir, barcode_to_variant):
    """
    Organises {sample}_{barcode}.fastq files into subdirectories based on the variants they contain

    Parameters:
    * output_dir (str): Directory containing the {sample}_{barcode}.fastq files to be organised
    * barcode_to_variant (dict): A dictionary mapping barcodes (str) to variants (str)

    Workflow:
    * Defines specific directories for 'wildtype' and ensure they exist
    * Iterates over each file in the output directory.
    * Moves 'wildtype.fastq' files to the 'wildtype' subdirectory
    * For other files, extracts the barcode from the filename
    * Maps the barcode to the corresponding variant using the barcode_to_variant dictionary
    * Move the file to a subdirectory named after the variant, creating the subdirectory if it doesn't exist
    """
    wildtype_dir = os.path.join(output_dir, "wildtype")

    if not os.path.exists(wildtype_dir):
        os.makedirs(wildtype_dir)

    for file in os.listdir(output_dir):
        if file == 'wildtype.fastq':
            shutil.move(os.path.join(output_dir, file), os.path.join(wildtype_dir, file))
        else:
            barcode = file.split('.')[0]
            variant = barcode_to_variant.get(barcode)
            if variant:
                variant_dir = os.path.join(output_dir, variant)
                if not os.path.exists(variant_dir):
                    os.makedirs(variant_dir)
                shutil.move(os.path.join(output_dir, file), os.path.join(variant_dir, file))

def process_barcode_fastq_directory(input_dir, barcodes_file, variant_info_csv):
    """
    Processes all {sample}_{barcode}.fastq files using specified barcode patterns with parallel processing and organises the resulting files by variant

    Parameters:
    * input_dir (str): Path to the input directory containing {sample}_{barcode}.fastq files
    * barcodes_file (str): Path to unique_SNV_barcodes.txt file containing barcodes that are each written on a separate line
    * variant_info_csv (str): Path to variants_info.csv containing variant information

    Workflow:
    * Creates a barcode-to-variant mapping using create_barcode_variant_mapping
    * Lists all FASTQ files in the input_dir and prepares tasks for processing each FASTQ file, ensuring that the output directories exist
    * Uses a ProcessPoolExecutor to process each FASTQ file in parallel, submitting tasks to process_sample_fastq_file
    * Organises resulting files by variant using organise_files_by_variant after processing each file
    """
    barcode_to_variant = create_barcode_variant_mapping(variant_info_csv)
    tasks = []
    print("Listing files in directory:", input_dir)
    for file in os.listdir(input_dir):
        if file.endswith(".fastq"):
            fastq_file = os.path.join(input_dir, file)
            input_name = os.path.splitext(file)[0]
            output_dir = os.path.join(sample_barcode_demux_dir, input_name)
            os.makedirs(output_dir, exist_ok=True)
            tasks.append((fastq_file, barcodes_file, output_dir))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {executor.submit(process_sample_fastq_file, *task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                future.result()
                organise_files_by_variant(task[2], barcode_to_variant)
            except Exception as e:
                print(f"Task failed: {task}, error: {e}")


################################################
### ALIGNING {BARCODE}.FASTQ FILES WITH GMAP ###
################################################

def setup_and_build_gmap_database(gmap_path, reference_fasta):
    """Set up and build GMAP database"""
    database_dir = os.path.join(gmap_path, "GMAP_database")

    ### Ensure the database directory exists ###
    if not os.path.exists(database_dir):
        os.makedirs(database_dir, exist_ok=True)

    ### Check if existing database directory name already exists, if it does, remove and rebuild ###
    database_path = os.path.join(database_dir, database_name)

    if os.path.exists(database_path):
        print(f"GMAP database {database_name} already exists. Deleting the existing database...")
        shutil.rmtree(database_path)

    print("Building GMAP database...")
    cmd = [
        'gmap_build',
        '-D', database_dir,
        '-d', database_name,
        reference_fasta
    ]
    subprocess.run(cmd, check=True)
    print("Database built successfully.")

    return database_dir, database_name

def align_fastq_with_gmap(input_fastq, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name):
    """Align FASTQ files using GMAP and save the output as *.sam, *.gff3 and *.err files

    Parameters:
    * input_fastq (str): Path to the input FASTQ file to be aligned
    * gmap_sam_output_dir (str): Directory where the *.sam output files will be saved
    * gmap_gff3_output_dir (str): Directory where *.gff3 output files will be saved
    * gmap_gff3_err_dir (str): Directory where *.err output files will be saved
    * database_dir (str): Path to the directory containing the GMAP database
    * database_name (str): Name of the GMAP database

    Workflow:
    * Construct the output filenames for SAM, GFF3, and error files based on the input FASTQ filename
    * Ensure the output directories exist by creating them if necessary
    * Construct the GMAP command for SAM alignment
    * Execute the SAM alignment command, redirecting output and error logs to the specified files
    * Construct the GMAP command for GFF3 alignment
    * Execute the GFF3 alignment command, appending error logs to the specified file
    * Return a success message if both alignments complete successfully, or an error message if an exception occurs
"""
    base_filename = os.path.basename(input_fastq).replace('.fastq', '')
    output_sam = os.path.join(gmap_sam_output_dir, f"{base_filename}.sam")
    output_gff3 = os.path.join(gmap_gff3_output_dir, f"{base_filename}.gff3")
    output_errors = os.path.join(gmap_gff3_err_dir, f"{base_filename}.err")

    directories = [gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)

    ### Create gmap command for sam output ###
    cmd_sam = [
        'gmap',
        '-n1',                                                                                                                  ### report only the best alignment ###
        '-t', '1',                                                                                                              ### set number of threads to 1 ###
        '-f', 'samse',                                                                                                          ### define output format for SAM (single-end reads) ###
        '-D', database_dir,                                                                                                     ### path to the GMAP database directory ###
        '-d', database_name,                                                                                                    ### database name ###
        input_fastq
    ]

    ### Create gmap command for gff3 and err files output based on original mini-isoqlr documentation ###
    cmd_gff3 = [
        'gmap',
        '-n1',
        '-t', '1',
        '--cross-species',
        '--gff3-add-separators=0',
        '-f', '2',
        '-z', 'auto',
        '-D', database_dir,
        '-d', database_name,
        input_fastq
    ]

    ### Execute the SAM alignment and handle outputs ###
    try:
        with open(output_sam, 'w') as sam_file, open(output_errors, 'w') as err_file:
            subprocess.run(cmd_sam, stdout=sam_file, stderr=err_file, check=True)

        ### Execute the GFF3 alignment and handle outputs ###
        with open(output_gff3, 'w') as gff3_file, open(output_errors, 'a') as err_file:
            subprocess.run(cmd_gff3, stdout=gff3_file, stderr=err_file, check=True)
    except subprocess.CalledProcessError as e:
        return f"Error processing {input_fastq}: {e}"

    return "Alignment completed successfully"

def align_wrapper(args):
    """
    Wrapper function for parallel processing with GMAP alignment.

    Parameters:
    * args (tuple): A tuple containing the following elements:
        - input_fastq (str): Path to the input FASTQ file to be aligned
        - gmap_sam_output_dir (str): Directory where the SAM output files will be saved
        - gmap_gff3_output_dir (str): Directory where the GFF3 output files will be saved
        - gmap_gff3_err_dir (str): Directory where the error files will be saved
        - database_dir (str): Path to the directory containing the GMAP database
        - database_name (str): Name of the GMAP database

    Returns:
    * The result of the align_fastq_with_gmap function call

    Workflow:
    * Unpack the tuple args into individual variables
    * Call the align_fastq_with_gmap function with the unpacked arguments
    * Return the result of the align_fastq_with_gmap function
    """
    input_fastq, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name = args
    return align_fastq_with_gmap(input_fastq, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name)

def process_alignments(sample_barcode_demux_dir, gmap_dir, database_dir, database_name):
    """
    Processes all FASTQ files for alignment using GMAP, except for 'unknown.fastq' files, 
    sets up the necessary output directories and uses parallel processing to align the FASTQ files

    Parameters:
    * sample_barcode_demux_dir (str): Directory containing the demultiplexed sample FASTQ files
    * gmap_dir (str): Base directory where the GMAP output files will be saved
    * database_dir (str): Path to the directory containing the GMAP database
    * database_name (str): Name of the GMAP database

    Workflow:
    * Find all FASTQ files in the sample_barcode_demux_dir directory recursively
    * Initialise a list to store tasks for parallel processing
    * Iterate over each FASTQ file:
        - skip the file if it is named 'unknown.fastq'
        - extract the sample and variant information from the file path
        - construct the output directories for SAM files, GFF3 files, and error files
        - append a tuple of arguments to the tasks list for each FASTQ file
        - Use ProcessPoolExecutor to process each FASTQ file in parallel by mapping the align_wrapper function to the tasks list
        - Collect and store the results of the alignment tasks
    """
    fastq_files = glob.glob(os.path.join(sample_barcode_demux_dir, '**', '*.fastq'), recursive=True)
    tasks = []
    for fastq_file in fastq_files:
        if 'unknown.fastq' in fastq_file:                                                                                       ### skip processing of unknown.fastq file ###
            continue
        parts = fastq_file.split(os.sep)
        sample = parts[-3]
        variant = parts[-2]
        gmap_sam_output_dir = os.path.join(gmap_dir, sample, variant, "sam_files")
        gmap_gff3_output_dir = os.path.join(gmap_dir, sample, variant, "gff3")
        gmap_gff3_err_dir = os.path.join(gmap_dir, sample, variant, "gff3_error_files")
        tasks.append((fastq_file, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(align_wrapper, tasks))

def delete_gmap_database(database_dir, database_name):
    """Delete the GMAP database directory

    Parameters:
    * database_dir (str): Path to the directory containing the GMAP database
    * database_name (str): Name of the GMAP database to delete

    Workflow:
    * Construct the database path and delete it if it exists
    """
    database_path = os.path.join(database_dir, database_name)

    if os.path.exists(database_path):
        print(f"Deleting GMAP database at {database_path}...")
        shutil.rmtree(database_path)
        print("GMAP database deleted successfully.")
    else:
        print(f"GMAP database at {database_path} does not exist or has already been deleted.")


####################################
### (ADAPTED) VIsoQLR PROCESSING ###
####################################

def gff3_data_loading(input_path):
    gff3 = pd.read_csv(input_path, sep='\t', header=None, comment='#', dtype=str)
    gff3.index = range(1, len(gff3) + 1)
    gff3 = gff3[gff3[2] == 'exon'][[0, 3, 4, 8, 5]]
    gff3.columns = ['gene', 'start', 'end', 'id', 'score']
    gff3['score'] = pd.to_numeric(gff3['score'])
    gff3['id'] = gff3['id'].str.extract(r'Name=([^;]+)')
    gff3['start'] = gff3['start'].astype(int)
    gff3['end'] = gff3['end'].astype(int)

    num_reads_post_trimming = gff3['id'].nunique()
    # print(gff3)

    return gff3, num_reads_post_trimming

def break_point_calculation(raw_exons, num_reads_post_trimming):
    start_count = raw_exons['start'].value_counts()
    end_count = raw_exons['end'].value_counts()

    start_count = start_count[start_count > (num_reads_post_trimming * breakpoint_freq_threshold / 100)]
    end_count = end_count[end_count > (num_reads_post_trimming * breakpoint_freq_threshold / 100)]

    start_position = pd.DataFrame({'breakp': start_count.index.astype(int), 'freq': start_count.values})
    end_position = pd.DataFrame({'breakp': end_count.index.astype(int), 'freq': end_count.values})
    start_position = start_position.sort_values(by=['breakp'], ascending=True)
    end_position = end_position.sort_values(by=['breakp'], ascending=True)
    start_position.index = range(1, len(start_position) + 1)
    end_position.index = range(1, len(end_position) + 1)

    if very_close_bp > 0:
        if breakpoint_padding >= very_close_bp:
            start_pos_to_remove = []
            for i in range(len(start_position)):
                condition = (
                    (start_position.loc[start_position.index[i], 'breakp'] + very_close_bp >= start_position['breakp']) &
                    (start_position.loc[start_position.index[i], 'breakp'] - very_close_bp <= start_position['breakp']) &
                    (start_position.loc[start_position.index[i], 'freq'] < start_position['freq'])
                )
                if condition.any():
                    start_pos_to_remove.append(start_position.index[i])
            start_position.drop(start_pos_to_remove, inplace=True)

            end_pos_to_remove = []
            for i in range(len(end_position)):
                condition = (
                    (end_position.loc[end_position.index[i], 'breakp'] + very_close_bp >= end_position['breakp']) &
                    (end_position.loc[end_position.index[i], 'breakp'] - very_close_bp <= end_position['breakp']) &
                    (end_position.loc[end_position.index[i], 'freq'] < end_position['freq'])
                )
                if condition.any():
                    end_pos_to_remove.append(end_position.index[i])
            end_position.drop(end_pos_to_remove, inplace=True)
        else:
            raise ValueError("-v, --very_close_bp cannot be bigger than -p, --padding")

    start_position['left'] = start_position['breakp'] - breakpoint_padding
    start_position['right'] = start_position['breakp'] + breakpoint_padding
    end_position['left'] = end_position['breakp'] - breakpoint_padding
    end_position['right'] = end_position['breakp'] + breakpoint_padding

    for i in range(len(start_position) - 1):
        if start_position.loc[start_position.index[i], 'right'] > start_position.loc[start_position.index[i + 1], 'left']:
            media_pos = (start_position.loc[start_position.index[i], 'right'] + start_position.loc[start_position.index[i + 1], 'left']) / 2
            if media_pos.is_integer():
                start_position.at[start_position.index[i], 'right'] = media_pos - 1
                start_position.at[start_position.index[i + 1], 'left'] = media_pos + 1
            else:
                start_position.at[start_position.index[i], 'right'] = np.floor(media_pos)
                start_position.at[start_position.index[i + 1], 'left'] = np.ceil(media_pos)

    for i in range(len(end_position) - 1):
        if end_position.loc[end_position.index[i], 'right'] > end_position.loc[end_position.index[i + 1], 'left']:
            media_pos = (end_position.loc[end_position.index[i], 'right'] + end_position.loc[end_position.index[i + 1], 'left']) / 2
            if media_pos.is_integer():
                end_position.at[end_position.index[i], 'right'] = media_pos - 1
                end_position.at[end_position.index[i + 1], 'left'] = media_pos + 1
            else:
                end_position.at[end_position.index[i], 'right'] = np.floor(media_pos)
                end_position.at[end_position.index[i + 1], 'left'] = np.ceil(media_pos)

    start_position = start_position[['breakp', 'left', 'right']]
    end_position = end_position[['breakp', 'left', 'right']]

    return start_position, end_position

def break_point_assignment(raw_exons, break_points_list):
    start_position, end_position = break_points_list
    bp_assigned_exons = raw_exons.copy()

    bp_assigned_exons['start_tag'] = np.nan
    for idx in start_position.index:
        condition = (
            (bp_assigned_exons['start'] >= start_position.loc[idx, 'left']) &
            (bp_assigned_exons['start'] <= start_position.loc[idx, 'right'])
        )
        bp_assigned_exons.loc[condition, 'start_tag'] = start_position.loc[idx, 'breakp']

    bp_assigned_exons['end_tag'] = np.nan
    for idx in end_position.index:
        condition = (
            (bp_assigned_exons['end'] >= end_position.loc[idx, 'left']) &
            (bp_assigned_exons['end'] <= end_position.loc[idx, 'right'])
        )
        bp_assigned_exons.loc[condition, 'end_tag'] = end_position.loc[idx, 'breakp']

    bp_assigned_exons['start_tag'] = bp_assigned_exons['start_tag'].fillna(pd.NA).astype('Int64')
    bp_assigned_exons['end_tag'] = bp_assigned_exons['end_tag'].fillna(pd.NA).astype('Int64')
    # print("bp_assigned_exons:", bp_assigned_exons)
    return bp_assigned_exons

def break_point_info_calculator(break_points_list, num_reads_post_trimming, bp_assigned_exons):
    start_position_in = break_points_list[0].copy()
    end_position_in = break_points_list[1].copy()

    start_position_in['type'] = 'start'
    end_position_in['type'] = 'end'

    ### merge start positions ###
    start_counts = bp_assigned_exons['start'].value_counts().reset_index()
    start_counts.columns = ['breakp', 'Freq']
    start_position = start_position_in.merge(start_counts, on='breakp', how='left').fillna(0)

    ### merge end positions ###
    end_counts = bp_assigned_exons['end'].value_counts().reset_index()
    end_counts.columns = ['breakp', 'Freq']
    end_position = end_position_in.merge(end_counts, on='breakp', how='left').fillna(0)

    ### merge start tags ###
    start_tag_counts = bp_assigned_exons['start_tag'].value_counts().reset_index()
    start_tag_counts.columns = ['breakp', 'Freq']
    start_position = start_position.merge(start_tag_counts, on='breakp', how='left', suffixes=('_x', '_y')).fillna(0)

    ### merge end tags ###
    end_tag_counts = bp_assigned_exons['end_tag'].value_counts().reset_index()
    end_tag_counts.columns = ['breakp', 'Freq']
    end_position = end_position.merge(end_tag_counts, on='breakp', how='left', suffixes=('_x', '_y')).fillna(0)

    ### combine start and end positions ###
    start_end_position = pd.concat([start_position, end_position], ignore_index=True)
    start_end_position['relative_freq'] = round(start_end_position['Freq_x'] / num_reads_post_trimming * 100, 2)
    start_end_position['relative_freq_padding'] = round(start_end_position['Freq_y'] / num_reads_post_trimming * 100, 2)
    start_end_position['increase_percentage'] = round((start_end_position['Freq_y'] - start_end_position['Freq_x']) / start_end_position['Freq_x'] * 100, 2)
    start_end_position['gene'] = bp_assigned_exons['gene'].iloc[0]
    start_end_position = start_end_position[['type', 'gene', 'breakp', 'left', 'right', 'Freq_x', 'relative_freq', 'Freq_y', 'relative_freq_padding', 'increase_percentage']]
    start_end_position.columns = ['#type', 'gene', 'breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads_breakpoint', 'percent_of_reads_breakpoint', 'number_of_reads_interval', 'percent_of_reads_interval', 'num_of_reads_increase_percentage_interval']

    ### create formatted read counts ###
    start_end_position['num_of_reads (%)'] = start_end_position['num_of_reads_breakpoint'].astype(str) + ' (' + start_end_position['percent_of_reads_breakpoint'].astype(str) + '%)'
    start_end_position['total_num_of_reads (%)'] = start_end_position['number_of_reads_interval'].astype(str) + ' (' + start_end_position['percent_of_reads_interval'].astype(str) + '%)'

    start_output = start_end_position[start_end_position['#type'] == 'start'][['breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads (%)', 'total_num_of_reads (%)']]
    end_output = start_end_position[start_end_position['#type'] == 'end'][['breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads (%)', 'total_num_of_reads (%)']]

    start_end_position = start_end_position[['#type', 'gene', 'breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads_breakpoint', 'percent_of_reads_breakpoint', 'number_of_reads_interval', 'percent_of_reads_interval', 'num_of_reads_increase_percentage_interval']]

    return start_output, end_output, start_end_position

def exon_definition(bp_assigned_exons):
    exon_id_filter1 = bp_assigned_exons[bp_assigned_exons.isna().any(axis=1)]['id']
    exon_id_filter2 = bp_assigned_exons[bp_assigned_exons['start_tag'] > bp_assigned_exons['end_tag']]['id']
    dup1 = bp_assigned_exons.duplicated(subset=['id', 'start_tag'], keep=False)
    dup2 = bp_assigned_exons.duplicated(subset=['id', 'end_tag'], keep=False)
    exon_id_filter3 = bp_assigned_exons[dup1 | dup2]['id']
    exon_id_filter = pd.concat([exon_id_filter1, exon_id_filter2, exon_id_filter3]).unique()

    defined_exons = bp_assigned_exons[~bp_assigned_exons['id'].isin(exon_id_filter)]

    defined_exons['coordinates'] = defined_exons['start_tag'].astype(str) + '-' + defined_exons['end_tag'].astype(str)

    return defined_exons, exon_id_filter

def exon_info_fun(defined_exons):
    exon_info_df = defined_exons['coordinates'].value_counts().reset_index()
    exon_info_df.columns = ['exon', 'number_of_reads']
    exon_info_df['percent_of_reads'] = round(exon_info_df['number_of_reads'] / defined_exons['id'].nunique() * 100, 1)
    exon_info_df = exon_info_df.sort_values('exon')

    exon_info_df['Size'] = exon_info_df['exon'].str.split('-', expand=True).apply(lambda x: int(float(x[1])) - int(float(x[0])) + 1, axis=1)
    exon_info_df.index = range(1, len(exon_info_df) + 1)
    exon_info_df = exon_info_df[["exon", "Size", "number_of_reads", "percent_of_reads"]]

    return exon_info_df

def isoform_definition(defined_exons):
    read_list = defined_exons.groupby('id')['coordinates'].apply(lambda x: '_'.join(sorted(x)))
    read_isoform = read_list.values

    isoform_frequencies = read_list.value_counts().reset_index()
    isoform_frequencies.columns = ['read_isoform', 'Freq']
    isoform_frequencies['perc'] = round(isoform_frequencies['Freq'] * 100 / isoform_frequencies['Freq'].sum(), 1)

    return isoform_frequencies, read_isoform

def isoform_full_length_filter(isoforms, break_points):
    isoform_info, isoform_read = isoforms
    start = str(break_points[0]['breakp'].min())
    end = str(break_points[1]['breakp'].max())

    positions_split = isoform_info['read_isoform'].str.split('-|_', expand=False)
    keep = []
    discard = []

    for i in range(len(isoform_info)):
        if start in positions_split[i] and end in positions_split[i]:
            keep.append(i)
        else:
            discard.append(i)

    keep_isoform_info = isoform_info.iloc[keep]
    keep_isoform_info['perc'] = round(keep_isoform_info['Freq'] / keep_isoform_info['Freq'].sum() * 100, 2)
    keep_isoform_read = isoform_read[np.isin(isoform_read, keep_isoform_info['read_isoform'])]
    partial_reads = isoform_read[np.isin(isoform_read, isoform_info.iloc[discard]['read_isoform'])]

    return keep_isoform_info, keep_isoform_read, partial_reads

def sort_isoform_segments(isoform):
    segments = isoform.split('_')
    sorted_segments = sorted(segments, key=lambda x: int(x.split('-')[0]))
    return '_'.join(sorted_segments)

def calculate_size(isoform):
    segments = isoform.split('_')
    size = sum([(int(part.split('-')[1]) - int(part.split('-')[0]) + 1) for part in segments])
    return size

def isoform_info_fun(isoform_frequencies, n_inicial, n_post_trim, n_final, n_not_full_length=None):
    if n_not_full_length is None:
        n_vector_reads = n_inicial - n_post_trim
        n_no_consensus_breakpoint_reads = n_post_trim - n_final

        other_groups_df = pd.DataFrame({
            'isoform_id': ["Only vector reads", "No consensus breakpoint reads"],
            'read_isoform': ["-", "-"],
            'Freq': [n_vector_reads, n_no_consensus_breakpoint_reads],
            'perc': ["-", "-"],
            'Size': [0, 0]
        })
    else:
        n_vector_reads = n_inicial - n_post_trim
        n_no_consensus_breakpoint_reads = n_post_trim - n_final - n_not_full_length

        other_groups_df = pd.DataFrame({
            'isoform_id': ["Only vector reads", "No consensus breakpoint reads", "Partial length reads"],
            'read_isoform': ["-", "-", "-"],
            'Freq': [n_vector_reads, n_no_consensus_breakpoint_reads, n_not_full_length],
            'perc': ["-", "-", "-"],
            'Size': [0, 0, 0]
        })

    isoform_frequencies['isoform_id'] = ["Iso" + str(i + 1) for i in range(len(isoform_frequencies))]
    isoform_frequencies['Size'] = isoform_frequencies['read_isoform'].str.split('-|_', expand=False).apply(
        lambda x: sum([-int(y) + 1 for y in x if y.isdigit()])
    )

    all_groups = pd.concat([other_groups_df, isoform_frequencies], ignore_index=True)
    all_groups['perc_total'] = round(all_groups['Freq'] * 100 / all_groups['Freq'].sum(), 2)
    all_groups.columns = ["Isoform_id", "Isoform", "Size", "Number_of_Reads", "Prerc_partial", "Perc_total"]
    all_groups = all_groups[["Isoform_id", "Isoform", "Size", "Number_of_Reads"]]
    all_groups = all_groups.iloc[2:]
    all_groups = all_groups.rename(columns={'Size': 'Number_of_Reads', 'Number_of_Reads': 'Percentage'})
    all_groups['Isoform'] = all_groups['Isoform'].apply(sort_isoform_segments)
    all_groups['Size'] = all_groups['Isoform'].apply(calculate_size)
    all_groups = all_groups[["Isoform_id", "Isoform", "Size", "Number_of_Reads", "Percentage"]]
    all_groups = all_groups.reset_index(drop=True)
    all_groups.index = range(1, len(all_groups)+1)

    return all_groups

def process_gff3_file(gff3_path, sample, variant, barcode):
    try:
        ### load GFF3 data ###
        gff3_data, num_reads_post_trimming = gff3_data_loading(gff3_path)
        if gff3_data.empty:
            logging.warning(f"GFF3 data is empty for {sample}/{variant}/{barcode}")
            return None

        ### calculate break points ###
        break_points_list = break_point_calculation(gff3_data, num_reads_post_trimming)

        ### assign break points to exons ###
        bp_assigned_exons = break_point_assignment(gff3_data, break_points_list)
        if bp_assigned_exons.empty:
            logging.warning(f"No break points assigned for {sample}/{variant}/{barcode}")
            return None

        ### calculate break point information ###
        start_output, end_output, start_end_position = break_point_info_calculator(break_points_list, num_reads_post_trimming, bp_assigned_exons)

        ### define exons ###
        defined_exons, exon_id_filter = exon_definition(bp_assigned_exons)
        if defined_exons.empty:
            logging.warning(f"Defined exons are empty for {sample}/{variant}/{barcode}")
            return None

        ### generate exon info df ###
        exon_info_df = exon_info_fun(defined_exons)

        ### define isoforms ###
        isoforms = isoform_definition(defined_exons)

        ### filter for full-length isoforms ###
        full_length_isoforms = isoform_full_length_filter(isoforms, break_points_list)

        ### prepare ID lists ###
        all_ids = gff3_data['id'].unique()
        no_consensus_ids = defined_exons[defined_exons['start_tag'] > defined_exons['end_tag']]['id'].unique()

        ### create dictionaries for classified IDs ###
        classified_ids = {k: v for k, v in zip(full_length_isoforms[1], range(len(full_length_isoforms[1])))}
        not_full_length = {k: v for k, v in zip(full_length_isoforms[2], range(len(full_length_isoforms[2])))}

        ### generate isoform dictionary ###
        iso_dict = isoform_info_fun(full_length_isoforms[0], num_reads_post_trimming, num_reads_post_trimming, len(classified_ids))
        if iso_dict.empty:
            logging.warning(f"Isoform dictionary is empty for {sample}/{variant}/{barcode}")
            return None

        ### save iso_dict df ###
        output_dir = os.path.join(visoqlr_dir, f"{sample}/{variant}")
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{sample}_{barcode}_iso_freq.csv")
        iso_dict.to_csv(output_path, index=False)

        return (sample, variant, barcode, iso_dict)

    except Exception as e:
        logging.error(f"Error processing {gff3_path}: {e}", exc_info=True)
        return None

def process_and_combine_csv_files(sample, variant, lock):
    with lock:  
        variant_dir = os.path.join(visoqlr_dir, sample, variant)
        csv_files = glob.glob(os.path.join(variant_dir, f"{sample}_*_iso_freq.csv"))
        dataframes = []
        for file in csv_files:
            df = pd.read_csv(file)

            if df.empty:
                print(f"Skipping empty file: {file}")
                continue

            ### check and drop columns if they exist ###
            for column in ["Isoform_id", "Percentage"]:
                if column in df.columns:
                    df = df.drop(columns=[column])

            dataframes.append(df)

        if dataframes:
            combined_df = pd.concat(dataframes, ignore_index=True)
            ### ensure the 'Size' column exists before applying groupby and aggregation ###
            if 'Size' in combined_df.columns:
                result_df = combined_df.groupby("Isoform").agg({"Size": "first", "Number_of_Reads": "sum"}).reset_index()
                result_df = result_df.sort_values(by=['Number_of_Reads'], ascending=False)
                total_reads = result_df['Number_of_Reads'].sum()
                result_df['Percentage'] = (result_df['Number_of_Reads'] / total_reads) * 100
                result_df['Isoform_id'] = ['Iso' + str(i + 1) for i in range(len(result_df))]

                output_path = os.path.join(variant_dir, f"{sample}_{variant}_combined.csv")
                result_df.to_csv(output_path, index=False)
            else:
                print(f"No 'Size' column found in files for {sample}_{variant}.")
        else:
            print(f"No valid CSV files found for {sample}_{variant}.")

def process_all_gff3_files(gmap_dir):
    lock = multiprocessing.Lock()  
    iso_results = []

    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for root, dirs, files in os.walk(gmap_dir):
                for file in files:
                    if file.endswith('.gff3'):
                        parts = root.split('/')
                        try:
                            sample = parts[-3]
                            variant = parts[-2]
                            barcode = file.split('.')[0]
                            gff3_path = os.path.join(root, file)
                            futures.append(executor.submit(process_gff3_file, gff3_path, sample, variant, barcode))
                        except IndexError as e:
                            logging.error(f"Error extracting sample/variant/barcode from path {root}: {e}")

            ### wait for all iso_dicts to be saved before proceeding ###
            for future in as_completed(futures):
                result = future.result()
                if result:
                    iso_results.append(result)

        ### combine CSV files for each variant once, ensuring no race conditions ###
        variants_processed = set((sample, variant) for sample, variant, _, _ in iso_results)
        for sample, variant in variants_processed:
            process_and_combine_csv_files(sample, variant, lock)

    except Exception as e:
        logging.error(f"Error processing GFF3 files in {gmap_dir}: {e}", exc_info=True)

def process_iso_freq_csv_files(visoqlr_dir, visoqlr_no_cdna_dir):
    """
    Processes all {sample}_{barcode}_iso_freq.csv files, removing specified columns and rows,
    recalculating the 'Percentage' column, and adding back an 'Isoform_id' column.
    The processed files are saved in a separate directory with the same filenames.
    """
    for root, dirs, files in os.walk(visoqlr_dir):
        for file in files:
            if file.endswith('_iso_freq.csv'):
                ### extract file path and sample/barcode info ###
                csv_path = os.path.join(root, file)
                try:
                    ### process each {sample}_{barcode}_iso_freq.csv file ###
                    df = pd.read_csv(csv_path)
                    df = df.drop(columns=['Isoform_id', 'Percentage'])                                                          ### drop columns if they exist ###
                    df = df[df['Isoform'] != f'{ref_start_pos}-{ref_end_pos}']                                                  ### remove cDNA ###         

                    ### recalculate 'Percentage' ###
                    total_reads = df['Number_of_Reads'].sum()
                    df['Percentage'] = (df['Number_of_Reads'] / total_reads) * 100
                    df = df.sort_values(by=['Percentage'], ascending=False)

                    ### add back 'Isoform_id' column with values like 'Iso1', 'Iso2', etc. ###
                    df['Isoform_id'] = ['Iso' + str(i + 1) for i in range(len(df))]
                    df = df[['Isoform_id', 'Isoform', 'Size', 'Number_of_Reads', 'Percentage']]

                    ### Define new output path in visoqlr_no_cdna_dir ###
                    relative_path = os.path.relpath(root, visoqlr_dir)
                    output_dir = os.path.join(visoqlr_no_cdna_dir, relative_path)
                    os.makedirs(output_dir, exist_ok=True)
                    output_path = os.path.join(output_dir, file)
                    df.to_csv(output_path, index=False)

                except Exception as e:
                    logging.error(f"Error processing {csv_path}: {e}", exc_info=True)

def process_and_combine_iso_freq_files(visoqlr_no_cdna_dir, visoqlr_combined_bc_dir):
    """
    Process and combine all {sample}_{variant}_iso_freq.csv files in visoqlr_no_cdna_dir,
    group by 'Isoform', sum 'Number_of_Reads', recalculate 'Percentage', and save the result
    in visoqlr_combined_bc_dir.
    """
    for root, dirs, files in os.walk(visoqlr_no_cdna_dir):
        ### group CSVs by their {sample}/{variant} directories ###
        csv_files = [file for file in files if file.endswith('_iso_freq.csv')]

        if csv_files:
            variant = root.split('/')[-1]
            sample = root.split('/')[-2]

            combined_df = pd.DataFrame()

            ### concat all {sample}_{variant}_iso_freq.csv files in the same {sample}/{variant} directory ###
            for file in csv_files:
                csv_path = os.path.join(root, file)
                try:
                    df = pd.read_csv(csv_path)
                    df = df.drop(columns=['Isoform_id', 'Percentage'])
                    combined_df = pd.concat([combined_df, df], ignore_index=True)

                except Exception as e:
                    logging.error(f"Error processing {csv_path}: {e}", exc_info=True)

            ### group by 'Isoform', aggregate by summing 'Number_of_Reads' and taking the first 'Size' for each 'Isoform' ###
            combined_df = combined_df.groupby('Isoform').agg({
                'Size': 'first',
                'Number_of_Reads': 'sum'
            }).reset_index()

            ### recalculate 'Percentage' ###
            total_reads = combined_df['Number_of_Reads'].sum()
            combined_df['Percentage'] = (combined_df['Number_of_Reads'] / total_reads) * 100
            combined_df = combined_df.sort_values(by=['Percentage'], ascending=False)

            ### add 'Isoform_id' column (Iso1, Iso2, Iso3, etc.) ###
            combined_df['Isoform_id'] = ['Iso' + str(i + 1) for i in range(len(combined_df))]
            combined_df = combined_df[['Isoform_id', 'Isoform', 'Size', 'Number_of_Reads', 'Percentage']]

            ### define the output file path and save the combined df as {sample}_{variant}_iso_freq.csv ###
            output_file = f'{sample}_{variant}_iso_freq.csv'
            output_path = os.path.join(visoqlr_combined_bc_dir, output_file)
            combined_df.to_csv(output_path, index=False)


##################################
### (ADAPTED) VISOQLR PLOTTING ###                                                                
##################################

### Set constant bar height and fixed figure height per isoform for plots ###
BAR_HEIGHT = 0.4  
FIGURE_HEIGHT_PER_ISOFORM = 0.66 
BREAKPOINT_PLOT_HEIGHT = 3.5 

### Function to parse isoform ###
def parse_isoform(isoform):
    segments = isoform.split('_')
    return [(int(segment.split('-')[0]), int(segment.split('-')[1])) for segment in segments]

### Function to plot isoform data from CSV files ###
def plot_isoform_file(df, ax, filename):
    df['Segments'] = df['Isoform'].apply(parse_isoform)

    num_segments = len(set(segment for segments in df['Segments'] for segment in segments))
    color_palette = sns.color_palette("magma_r", num_segments)
    color_map = {segment: color_palette[i] for i, segment in enumerate(sorted(set(segment for segments in df['Segments'] for segment in segments)))}

    num_segments = len(set(segment for segments in df['Segments'] for segment in segments))
    color_palette = sns.color_palette("magma_r", num_segments)
    color_map = {}
    for i, segment in enumerate(sorted(set(segment for segments in df['Segments'] for segment in segments))):
        color_map[segment] = color_palette[i]

    df = df[::-1].reset_index(drop=True)

    for i, row in df.iterrows():
        y_pos = i + 1
        segments = sorted(row['Segments'])
        for segment in segments:
            start, end = segment
            ax.barh(y_pos, end - start, left=start, height=BAR_HEIGHT, color=color_map[segment])
        for j in range(len(segments) - 1):
            line_x = [segments[j][1], segments[j + 1][0]]
            line_y = [y_pos, y_pos]
            ax.plot(line_x, line_y, color='black', linewidth=0.7)

    ### Set consistent y-limits to ensure equal bar height across all plots ###
    ax.set_ylim(0.5, len(df) + 0.5)

    ### Set y-ticks and labels ###
    ax.set_yticks(range(1, len(df) + 1))
    ax.set_yticklabels([f'Isoform {len(df) - i + 1} | {row.Size}bp | {row.Number_of_Reads} | {row.Percentage:.2f}%' for i, row in enumerate(df.itertuples(), 1)], fontsize=8)
    ax.set_xticks([])
    ax.set_xticklabels([])

    ### Remove spines ###
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('lightgrey')
    ax.spines['bottom'].set_color('lightgrey')

    ### Add legend ###
    handles = [mpatches.Patch(color=color_map[segment], label=f'{segment[0]}-{segment[1]}') for segment in sorted(color_map)]
    ax.legend(handles=handles, title="Exon Coordinates", bbox_to_anchor=(1.02, 1.11), loc='upper left', fontsize=8, title_fontsize=9)

    ### Add labels ###
    ax.set_xlim(ref_start_pos - 21, ref_end_pos + 22)
    ax.set_ylabel('Isoform ID | Size | Count | Percentage', labelpad=19, fontsize=10)
    ax.set_title(f"Isoforms of {filename.split('_iso_freq.csv')[0].replace('_', ' ')}", fontsize=12, pad=20)

### Function to plot breakpoints from GFF3 ###
def plot_breakpoints_file(breakpoints_df, ax, sample, variant):
    colors = {'start': 'crimson', 'end': 'mediumblue'}
    ax.bar(breakpoints_df['Breakpoints'], breakpoints_df['Percentage'], color=breakpoints_df['Type'].map(colors), width=1.9)

    ### Set axes labels ###
    ax.set_xlabel('Position', fontsize=9)
    ax.set_ylabel('Percentage', fontsize=9)

    ### Set x-limits and tick labels ###
    ax.set_xlim(ref_start_pos - 21, ref_end_pos + 22)
    ax.set_ylim(0, 100)
    max_breakpoint = breakpoints_df['Breakpoints'].max()
    ax.set_xticks(range(0, max_breakpoint + 50, 50))
    ax.tick_params(axis='x', labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ### Plot dots on top of bars ### 
    for idx, row in breakpoints_df.iterrows():
        if row['Percentage'] > 10:
            dot_color = colors[row['Type']]
            ax.plot(row['Breakpoints'], row['Percentage'] + 1, 'o', color=dot_color, markersize=2)

    ax.tick_params(axis='y', labelsize=7)

    ### Add legend for the breakpoints ###
    handles = [mpatches.Patch(color='crimson', label='Start'),
               mpatches.Patch(color='mediumblue', label='End')]
    ax.legend(handles=handles, title="Breakpoint Type", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, title_fontsize=9)

### Main function to process files and create combined plots ###
def process_combined_plots(isoform_file, visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir):
    filename = os.path.basename(isoform_file)
    sample = filename.split('_')[0]
    variant = filename.split('_')[1]
    pdf_path = os.path.join(VIsoQLR_plots_dir, f'{sample}_{variant}.pdf')

    ### Define {sample}_{variant}'s breakpoints csv output path ###
    sample_variant_gff3_path = os.path.join(sample_variant_gff3_dir, f"{sample}_{variant}.csv")

    ### Check if the isoform file is empty before attempting to read it ###
    if os.stat(isoform_file).st_size == 0:
        print(f"Skipping {isoform_file}: File is empty.")
        return

    try:
        ### Plotting isoforms from isoform file ###
        df_isoform = pd.read_csv(isoform_file, sep=',')
        if df_isoform.empty:
            print(f"Skipping {isoform_file}: No data found.")
            return

        df_isoform = df_isoform[df_isoform['Percentage'] >= 1]

        num_isoforms = len(df_isoform)
        isoform_plot_height = max(FIGURE_HEIGHT_PER_ISOFORM, FIGURE_HEIGHT_PER_ISOFORM * num_isoforms)

        ### Total figure height = isoform plot height + constant breakpoint plot height ###
        total_figure_height = isoform_plot_height + BREAKPOINT_PLOT_HEIGHT

        ### Create subplots with the dynamic isoform plot height and fixed breakpoint plot height ###
        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(11, total_figure_height),
                                       gridspec_kw={'height_ratios': [isoform_plot_height, BREAKPOINT_PLOT_HEIGHT]})

        plot_isoform_file(df_isoform, ax1, filename)

    except pd.errors.EmptyDataError:
        print(f"Skipping {isoform_file}: No columns to parse.")
        return
    except Exception as e:
        print(f"Error processing {isoform_file}: {e}")
        return

    ### Process GFF3 files for the same {sample}/{variant} combination ###
    gff3_files = glob.glob(os.path.join(gmap_dir, f"{sample}/{variant}/gff3", "*.gff3"))
    if not gff3_files:
        print(f"No GFF3 files found for {sample}/{variant} in {gmap_dir}")
        return
    all_data = []
    for gff3_file in gff3_files:
        gff3 = pd.read_csv(gff3_file, sep='\t', header=None, comment='#', dtype=str)
        gff3.index = range(1, len(gff3) + 1)
        gff3 = gff3[gff3[2] == 'exon'][[3, 4, 8]]
        gff3.columns = ['start', 'end', 'id']
        gff3['id'] = gff3['id'].str.extract(r'Name=([^;]+)')
        gff3['start'] = gff3['start'].astype(int)
        gff3['end'] = gff3['end'].astype(int)
        all_data.append(gff3)

    combined_df = pd.concat(all_data, ignore_index=True)

    ### Calculate number of unique reads after combining data for each sample's variant ###
    number_of_unique_reads = combined_df['id'].nunique()
    combined_df = combined_df[["start", "end"]]

    ### Reshape data to long format ###
    starts = combined_df[['start']].rename(columns={'start': 'Breakpoints'})
    starts['Type'] = 'start'
    ends = combined_df[['end']].rename(columns={'end': 'Breakpoints'})
    ends['Type'] = 'end'
    breakpoints_df = pd.concat([starts, ends], ignore_index=True)

    ### Calculate frequency and percentage of each unique start and end breakpoint ###
    breakpoints_df = breakpoints_df.groupby(['Breakpoints', 'Type']).size().reset_index(name='Frequency')
    breakpoints_df['Percentage'] = (breakpoints_df['Frequency'] / number_of_unique_reads) * 100
    breakpoints_df.to_csv(sample_variant_gff3_path, index=False)

    ### Plot breakpoints plot on the bottom (ax2) ###
    plot_breakpoints_file(breakpoints_df, ax2, sample, variant)

    ### Save the combined plot to a PDF ###
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches='tight', pad_inches=1.31)
    plt.close(fig)

### Function to handle parallel processing of plots ###
def plot_iso_bp_parallel(visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir):
    isoform_files = [os.path.join(visoqlr_combined_bc_dir, file) for file in os.listdir(visoqlr_combined_bc_dir) if file.endswith('csv')]

    if isoform_files:
        print(f"Found {len(isoform_files)} CSV files to process.")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:  # Adjust the number of workers as needed
            futures = [executor.submit(process_combined_plots, isoform_file, visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir) for isoform_file in isoform_files]
            for future in futures:
                future.result()  # Ensure each task completes
    else:
        print(f"No CSV files found in directory: {visoqlr_combined_bc_dir}")

def plot_sample_breakpoints(gmap_dir, sample_gff3_dir, sample_breakpoint_freq_plot_dir):

    ### create a dictionary to group GFF3 files by sample ###
    sample_data = {}

    ### walk through each {sample}/{variant} directory ###
    for root, dirs, files in os.walk(gmap_dir):
        gff3_files = glob.glob(os.path.join(root, "*.gff3"))

        if gff3_files:
            ### extract sample name from the first GFF3 file ###
            sample = root.split(os.sep)[-3]

            ### process each GFF3 file in the current directory ###
            for gff3_file in gff3_files:
                ### check for and skip empty gff3 file ###
                if os.stat(gff3_file).st_size == 0:
                    print(f"Skipping empty file: {gff3_file}")
                    continue

                try:
                    gff3 = pd.read_csv(gff3_file, sep='\t', header=None, comment='#', dtype=str)
                    gff3.index = range(1, len(gff3) + 1)

                    ### filter and select columns ###
                    gff3 = gff3[gff3[2] == 'exon'][[3, 4, 8]]
                    gff3.columns = ['start', 'end', 'id']

                    ### extract the ID and convert start and end to integers ###
                    gff3['id'] = gff3['id'].str.extract(r'Name=([^;]+)')
                    gff3['start'] = gff3['start'].astype(int)
                    gff3['end'] = gff3['end'].astype(int)

                    ### append the data to the corresponding sample in the dictionary ###
                    if sample not in sample_data:
                        sample_data[sample] = []
                    sample_data[sample].append(gff3)

                except pd.errors.EmptyDataError:
                    print(f"Skipping {gff3_file}: No columns to parse.")
                    continue

    ### concat data and perform aggregation for each sample ###
    for sample, dataframes in sample_data.items():
        combined_df = pd.concat(dataframes, ignore_index=True)
        number_of_unique_reads = combined_df['id'].nunique()
        combined_df = combined_df[["start", "end"]]

        ### reshape data to long format ###
        starts = combined_df[['start']].rename(columns={'start': 'Breakpoints'})
        starts['Type'] = 'start'

        ends = combined_df[['end']].rename(columns={'end': 'Breakpoints'})
        ends['Type'] = 'end'

        ### combine start and end data into a single df ###
        breakpoints_df = pd.concat([starts, ends], ignore_index=True)

        ### calculate frequencies ###
        breakpoints_df = breakpoints_df.groupby(['Breakpoints', 'Type']).size().reset_index(name='Frequency')

        ### calculate percentage ###
        breakpoints_df['Percentage'] = (breakpoints_df['Frequency'] / number_of_unique_reads) * 100

        ### define the output file path and save aggregated df ###
        output_file = os.path.join(sample_gff3_dir, f"{sample}.csv")
        breakpoints_df.to_csv(output_file, index=False)

        ### use Plotly to create an interactive bar plot ###
        fig = make_subplots()

        ### separate data by type (start and end) for coloring ###
        for bp_type, color in zip(['start', 'end'], ['crimson', 'mediumblue']):
            type_data = breakpoints_df[breakpoints_df['Type'] == bp_type]
            fig.add_trace(go.Bar(
                x=type_data['Breakpoints'],
                y=type_data['Percentage'],
                name=f'{bp_type.capitalize()} Breakpoints',
                marker=dict(color=color),
                width=1.9,
                hovertemplate=
                    'Breakpoint: %{x}<br>' +
                    'Percentage: %{y:.2f}%<br>' +
                    f'Type: {bp_type.capitalize()}', 
            ))

            ### add dots on top of bars where "Percentage" > 10 ###
            high_percentage = type_data[type_data['Percentage'] > 10]
            fig.add_trace(go.Scatter(
                x=high_percentage['Breakpoints'],
                y=high_percentage['Percentage'], 
                mode='markers',
                marker=dict(color=color, size=6),
                hovertemplate=
                    'Breakpoint: %{x}<br>' +
                    'Percentage: %{y:.2f}%<br>' +
                    f'Type: {bp_type.capitalize()}', 
                name=f'{bp_type.capitalize()} Breakpoints', 
                showlegend=False  
            ))

        ### calculate the x-axis limits and ticks ###
        x_min = breakpoints_df["Breakpoints"].min() - 19
        x_max = breakpoints_df["Breakpoints"].max() + 19
        x_ticks = list(range(0, x_max + 50, 50))  

        ### customise layout ###
        fig.update_layout(
            title=f"Breakpoints Frequency for {sample}",
            xaxis=dict(
                title="Breakpoints",
                range=[x_min, x_max], 
                tickvals=x_ticks, 
                ticktext=[str(tick) for tick in x_ticks],  
            ),
            yaxis=dict(
                title="Percentage",
                range=[0, 100],  
            ),
            template="plotly_white",
            showlegend=True,
        )

        ### save plot as {sample}_breakpoint_frequency.html ###
        html_output_file = os.path.join(sample_breakpoint_freq_plot_dir, f"{sample}_breakpoint_frequency.html")
        fig.write_html(html_output_file)
        # print(f"Interactive plot saved to: {html_output_file}")

### codon_barcode file to remove ###
def remove_codon_file(root_output_dir, codon_filename="codon_barcodes_all_info.csv"):
    codon_file_path = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", codon_filename)
    if os.path.exists(codon_file_path):
        os.remove(codon_file_path)
    else:
        print(f"{codon_file_path} does not exist - nothing to remove.")

def main():
    create_directories()

    chopper_start_time = time.time()
    choppered_fastq_path = run_chopper(cDNA_fastq_file, chopper_cDNA_dir)
    chopper_end_time = time.time()
    print(f"Time taken to QC fastq with chopper: {chopper_end_time - chopper_start_time:.3f} seconds")

    demux_idx_start_time = time.time()
    demux_index(choppered_fastq_path, demux_index_dir, index_sequence_patterns)
    demux_idx_end_time = time.time()
    print(f"Time taken to demultiplex fastq by index: {demux_idx_end_time - demux_idx_start_time:.3f} seconds")

    demux_bc_start_time = time.time()
    process_barcode_fastq_directory(demux_index_dir, barcode_txt_file_path, variant_info_csv)
    demux_bc_end_time = time.time()
    print(f"Time taken to demultiplex sample.fastq by barcodes: {demux_bc_end_time - demux_bc_start_time:.3f} seconds")

    build_gmap_db_start_time = time.time()
    database_dir, database_name = setup_and_build_gmap_database(gmap_path, reference_fasta)
    build_gmap_db_end_time = time.time()
    print(f"Time taken to build gmap database: {build_gmap_db_end_time - build_gmap_db_start_time:.3f} seconds")

    gmap_start_time = time.time()
    process_alignments(sample_barcode_demux_dir, gmap_dir, database_dir, database_name)
    gmap_end_time = time.time()
    print(f"Time taken for gmap alignment: {gmap_end_time - gmap_start_time:.3f} seconds")

    id_iso_start_time = time.time()
    process_all_gff3_files(gmap_dir)
    id_iso_end_time = time.time()
    print(f"Time taken to ID and quantify isoforms per barcode.fastq: {id_iso_end_time - id_iso_start_time:.3f} seconds")

    filter_cDNA_start_time = time.time()
    process_iso_freq_csv_files(visoqlr_dir, visoqlr_no_cdna_dir)
    filter_cDNA_end_time = time.time()
    print(f"Time taken to filter cDNA from ID-ed isoforms and requantify %: {filter_cDNA_end_time - filter_cDNA_start_time:.3f} seconds")

    combine_var_iso_start_time = time.time()
    process_and_combine_iso_freq_files(visoqlr_no_cdna_dir, visoqlr_combined_bc_dir)
    combine_var_iso_end_time = time.time()
    print(f"Time taken to combine isoforms per variant and requantify %: {combine_var_iso_end_time - combine_var_iso_start_time:.3f} seconds")

    plot_iso_start_time = time.time()
    try:
        plot_iso_bp_parallel(visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir)
        plot_iso_end_time = time.time()
        print(f"Time taken to generate plot of isoforms per variant: {plot_iso_end_time - plot_iso_start_time:.3f} seconds")
    except Exception as e:
        print(f"Error in plot_iso_bp_parallel: {e}")

    plot_breakpoint_start_time = time.time()
    try:
        plot_sample_breakpoints(gmap_dir, sample_gff3_dir, sample_breakpoint_freq_plot_dir)
        plot_breakpoint_end_time = time.time()
        print(f"Time taken to generate plot of breakpoints per sample: {plot_breakpoint_end_time - plot_breakpoint_start_time:.3f} seconds")
    except:
        print(f"Error in plot_sample_breakpoints: {e}")

    delete_gmap_database(database_dir, database_name)
    remove_codon_file(root_output_dir)

if __name__ == "__main__":
    main()