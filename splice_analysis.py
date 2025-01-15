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
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import shlex
import logging
import multiprocessing
import threading


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
gmap_path = "/oak/stanford/groups/euan/projects/variant_effect_mapping/tools/GMAP/bin/"

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

def create_directories():
    """Create necessary directories if they don't exist"""
    global base_result_directory
    base_result_directory = os.path.join(root_output_dir, "splice_results")
    if not os.path.exists(base_result_directory):
        os.makedirs(base_result_directory, exist_ok=True)
    ### make these directories accessible globally ###
    global chopper_cDNA_dir, demux_index_dir, sample_barcode_demux_dir, gmap_dir, visoqlr_dir, visoqlr_no_cdna_dir, visoqlr_combined_bc_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir
    chopper_cDNA_dir = os.path.join(base_result_directory, "cDNA_QC_fastq")                          
    demux_index_dir = os.path.join(base_result_directory, "Sample_Fastq")
    sample_barcode_demux_dir = os.path.join(base_result_directory, "Sample_Barcode_Demux")
    gmap_dir = os.path.join(base_result_directory, "GMAP")
    visoqlr_dir = os.path.join(base_result_directory, "VIsoQLR")
    visoqlr_no_cdna_dir = os.path.join(base_result_directory, "VIsoQLR_no_cDNA")
    visoqlr_combined_bc_dir = os.path.join(base_result_directory, "VIsoQLR_sample_variant")
    VIsoQLR_plots_dir = os.path.join(base_result_directory, "VIsoQLR_plots")
    sample_variant_gff3_dir = os.path.join(base_result_directory, "sample_variant_breakpoints")
    for dir_path in [chopper_cDNA_dir, demux_index_dir, sample_barcode_demux_dir, gmap_dir, visoqlr_dir, visoqlr_no_cdna_dir, visoqlr_combined_bc_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir]:
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
    Searches and extracts a subsequence between specified 5' indices and 3' indices within each read sequence of input cDNA fastq, 
    returning the subsequence along with its quality scores.

    Parameters:
    * sequence (str): Each read seq of input cDNA fastq file
    * quality (str): The quality scores corresponding to each read seq
    * fwd_front (str): The 5' index sequence of the forward read (5' -> 3')
    * fwd_end (str): The 3' index sequence of the forward read (5' -> 3')
    * rev_front (str): The 5' index sequence of the reversed read (3' -> 5') ((i.e. the reverse complement of fwd_end))
    * rev_end (str): The 3' index sequence of the reversed read (3' -> 5') ((i.e. the reverse complement of fwd_front))

    Returns:
    * Tuple (extracted_seq (str), extracted_quality (str)):
        - The extracted cDNA subsequence and its quality scores if a match is found; otherwise, returns (None, None)

    Workflow:
    1. Contains inner is_close_match function that: 
        * checks if a subsequence matches a given pattern with up to a specified number of allowed mismatches
        * returns True if the number of mismatches is within the allowed limit, False otherwise

    2. Processing Forward Matches:
        * searches for `fwd_front` sequence starting from the beginning of `sequence`, storing index as `fwd_start_index` if 'fwd_front' match is found
        * searches for `fwd_end` sequence starting from the end of `sequence`, moving backwards and stores index as `fwd_end_index` if 'fwd_end' match is found
        * if both `fwd_front` and `fwd_end` are found and `fwd_start_index` < `fwd_end_index`, the subsequence and its quality scores between these indices are extracted and returned

    3. Processing Reverse Match (Executed on reads where forward matches are not found):
        * searches for `rev_front` sequence starting from the beginning of `sequence`, storing index as `rev_start_index` if 'rev_front' match is found
        * searches for `rev_end` sequence starting from the end of `sequence', moving backwards and stores index as `rev_end_index` if `rev_end` match is found
        * if both `rev_front` and `rev_end` are found and `rev_start_index` is < `rev_end_index`, the subsequence and its quality scores between these indices are extracted 
        * extracted sequences are reverse complemented and quality scores are reversed and returned

    4. Fallback:
        * if neither a forward nor a reverse match is found, the function returns `(None, None)` to indicate that no valid match was detected
    """ 
    def is_close_match(seq, pattern, allowed_mismatches=1):
        mismatches = sum(1 for a, b in zip(seq, pattern) if a != b)
        return mismatches <= allowed_mismatches

    ### Processing forward match ###
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

    ### Processing reverse match ###
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
    Processes a chunk of FASTQ data, trims sequences based on specified forward and reverse index sequences, and stores trimmed sequences and quality scores for each sample

    Parameters:
    * data_chunk (list of str): Chunk of FASTQ data containing multiple sequences and their quality scores
    * index_sequence_patterns (dict): Dictionary where each sample key contains a dict with the forward and reverse index sequences ('fwd_front', 'fwd_end', 'rev_front', 'rev_end').

    Returns:
    * results (dict): A dictionary with sample names as keys and lists of trimmed sequences and corresponding quality scores as values.

    Workflow:
    1. Initialize Output Dictionary:
        * create a dictionary `results` with keys corresponding to sample names from `index_sequence_patterns`, and initalise values as empty lists

    2. Process Each FASTQ Record:
        * convert the `data_chunk` list into an iterator `it` for sequential processing
        * enter a loop to process each FASTQ record within chunk

    3. Barcode Matching:
        * use `find_and_trim_sequence` function to search for specified index sequences in the `sequence`, extract the subsequence and corresponding quality scores
       
    4. Validation and Trimming:
        * ensure extracted sequence and quality scores are non-empty after finding a match
        * further trim 8 bases from the front and end of both `trimmed_seq` and `trimmed_quality` to obtain `final_trimmed_seq` and `final_trimmed_quality`
        * ensure that the final trimmed sequence and quality are still non-empty after trimming

    5. Store Results:
        * append FASTQ record (with `seqID`, `final_trimmed_seq`, `separator`, `final_trimmed_quality`) to list for sample in `results` dictionary if the final trimmed sequence and quality are valid

    6. Return Results:
        * return the `results` dictionary, which contains the demultiplexed and trimmed FASTQ entries for each sample 

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
    This function demultiplexes a FASTQ file containing cDNA sequences based on specified index sequence patterns, processing the file in chunks, 
    and uses parallel processing to improve efficiency, before finally writing the resulting sequences to separate FASTQ files for each sample i.e. demultiplexing

    Parameters:
    * cDNA_fastq_file (str): Path to the input FASTQ file containing cDNA sequences
    * demux_index_dir (str): Directory where the demultiplexed {sample}.fastq files will be saved
    * index_sequence_patterns (dict): A dictionary where keys are sample names and values are dictionaries containing the forward and reverse sequence patterns for demultiplexing; each value dictionary has the keys: fwd_front, fwd_end, rev_front, rev_end

    Workflow:
    1. Reads all lines in fastq file into 'data' list
    2. Splits the data list into smaller chunks of size chunk_size, where each chunk is processed independently
    3. Creates a temporary file for each sample to store the demultiplexed sequences before they are written to the final output files
    4. Uses ProcessPoolExecutor to process each chunk in parallel, where process_fastq_chunk function is submitted as a task for each chunk, and results are written to corresponding sample's temp file as each task completes
    5. Close and merge temporary files
        * closes each temp file
        * opens the final output file for each sample in append mode ('ab')
        * reads from the temp file and writes its contents to the final output file
        * deletes the temp file after merging its contents
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
    Reads barcode.txt file produced from process_variants.py containing barcode sequences to generate a dictionary of compiled regular expression patterns, 
    where each pattern is constructed by embedding each unique barcode between pre-defined prefix and suffix

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
    * data (str): A string containing FASTQ data, where each set of four lines represents a single read (seqID, sequence, separator, quality)
    * barcodes (dict): A dictionary where keys are barcodes (str) and values are compiled regex patterns (re.compile objects) to search within the sequences

    Returns:
    * outputs (defaultdict): A dictionary where keys are filenames ("{barcode}.fastq" or "wildtype.fastq") and values are lists of matched FASTQ entries

    Workflow:
    1. Creates a defaultdict of lists to store the classified FASTQ records where each key corresponds to a filename, and each value is a list of FASTQ records
    2. Splits input FASTQ data string into individual lines, stored in the lines list
    3. Iterates over the lines in chunks of four, extracting the seqID, sequence, '+', and quality scores for each read
    4. Match sequences to barcode and wildtype patterns:
        * initialises a matched flag to False
        * iterates over each barcode and its corresponding regex pattern in the barcodes dictionary
        * checks if the sequence matches the barcode regex using regex.search(sequence)
        * appends the FASTQ read to the corresponding output list (keyed by {barcode}.fastq), sets the matched flag to True, and breaks out of loop if barcode match is found
        * checks if the sequence matches the wildtype pattern using wt_regex.search(sequence) if no barcode match was found (matched is False)
        * appends the FASTQ read to the wildtype.fastq output list if wildtype match is found
    5. Returns the outputs dictionary containing demultipleded FASTQ reads
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
    Processes a FASTQ file in chunks to classify sequences based on pre-defined barcode patterns before writing the results to separate FASTQ files in specified output directory

    Parameters:
    * fastq_file (str): Path to the input FASTQ file containing sequence data
    * barcodes_file (str): Path to barcode.txt file containing barcodes that are each written on a separate line
    * output_dir (str): Directory where the demultiplexed FASTQ files will be saved

    Workflow:
    1. Calls read_barcodes function to read and compile regex patterns from barcode.txt file and store the resulting dictionary of barcode patterns in "barcodes" variable
    2. Creates a defaultdict of lists to store the classified FASTQ reads, where each key corresponds to a filename, and each value is a list of FASTQ reads
    3. Read and process FASTQ file in chunks:
        * opens the specified fastq_file in read mode  and reads the file in chunks of CHUNK_SIZE lines, storing each chunk in "data" variable
        * breaks out of the loop if data is empty (end of file)
        * calls the process_sample_fastq_chunk function to process the chunk of data using the barcodes dictionary
        * appends the resulting classified FASTQ reads to the corresponding lists in the final_outputs dictionary
    4. Write demultiplexed data to output files:
        * iterates over the final_outputs dictionary
        * constructs the output file path for each filename
        * opens the output file in append mode and writes the demuxed FASTQ reads
        * catches any exceptions that occur during the file reading and processing, printing an error message indicating the nature of the exception and the file being processed
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
    1. Reads variants_info.csv into a df
    2. Initialises an empty dictionary for mapping barcodes to variants
    3. Iterates over each row of the df and splits the 'BARCODE' column into individual barcodes
    4. Maps each barcode to the corresponding 'VARIANT' value in the dictionary
    5. Returns the barcode-to-variant mapping dictionary
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
    Organises FASTQ files in the output directory into subdirectories based on what variants they contain

    Parameters:
    * output_dir (str): Directory containing the FASTQ files to be organised
    * barcode_to_variant (dict): A dictionary mapping barcodes (str) to variants (str)

    Workflow:
    1. Defines specific directories for 'wildtype' and ensure they exist
    2. Iterates over each file in the output directory.
    3. Moves 'wildtype.fastq' files to the 'wildtype' subdirectory
    4. For other files, extracts the barcode from the filename
    5. Maps the barcode to the corresponding variant using the barcode_to_variant dictionary
    6. Move the file to a subdirectory named after the variant, creating the subdirectory if it doesn't exist
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
    This function processes all FASTQ files in an input directory using specified barcode patterns with parallel processing and organises the resulting files by variant

    Parameters:
    * input_dir (str): Path to the input directory containing FASTQ files
    * barcodes_file (str): Path to barcode.txt file containing barcodes that are each written on a separate line
    * variant_info_csv (str): Path to variants_info.csv containing variant information

    Workflow:
    1. Creates a barcode-to-variant mapping using create_barcode_variant_mapping
    2. Lists all FASTQ files in the input_dir and prepares tasks for processing each FASTQ file, ensuring that the output directories exist
    3. Uses a ProcessPoolExecutor to process each FASTQ file in parallel, submitting tasks to process_sample_fastq_file
    4. Organises resulting files by variant using organise_files_by_variant after processing each file
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
    1. Construct the output filenames for SAM, GFF3, and error files based on the input FASTQ filename
    2. Ensure the output directories exist by creating them if necessary
    3. Construct the GMAP command for SAM alignment
    4. Execute the SAM alignment command, redirecting output and error logs to the specified files
    5. Construct the GMAP command for GFF3 alignment
    6. Execute the GFF3 alignment command, appending error logs to the specified file
    7. Return a success message if both alignments complete successfully, or an error message if an exception occurs
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
    A wrapper function that unpacks arguments and calls the align_fastq_with_gmap function that serves parallel processing 
    where the executor.map function requires a function that accepts a single argument.

    Parameters:
    * args (tuple): A tuple containing the following elements:
        - input_fastq (str): Path to the input FASTQ file to be aligned
        -  gmap_sam_output_dir (str): Directory where the SAM output files will be saved
        - gmap_gff3_output_dir (str): Directory where the GFF3 output files will be saved
        - gmap_gff3_err_dir (str): Directory where the error files will be saved
        - database_dir (str): Path to the directory containing the GMAP database
        - database_name (str): Name of the GMAP database

    Returns:
    * The result of the align_fastq_with_gmap function call

    Workflow:
    1. Unpack the tuple args into individual variables
    2. Call the align_fastq_with_gmap function with the unpacked arguments
    3. Return the result of the align_fastq_with_gmap function
    """
    input_fastq, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name = args
    return align_fastq_with_gmap(input_fastq, gmap_sam_output_dir, gmap_gff3_output_dir, gmap_gff3_err_dir, database_dir, database_name)

def process_alignments(sample_barcode_demux_dir, gmap_dir, database_dir, database_name):
    """
    This function processes all FASTQ files for alignment using GMAP, except for 'unknown.fastq' files, 
    sets up the necessary output directories and uses parallel processing to align the FASTQ files

    Parameters:
    * sample_barcode_demux_dir (str): Directory containing the demultiplexed sample FASTQ files
    * gmap_dir (str): Base directory where the GMAP output files will be saved
    * database_dir (str): Path to the directory containing the GMAP database
    * database_name (str): Name of the GMAP database

    Workflow:
    1. Find all FASTQ files in the sample_barcode_demux_dir directory recursively
    2. Initialise a list to store tasks for parallel processing
    3. Iterate over each FASTQ file (excluding unknown.fastq):
        * extract the sample and variant information from the file path
        * construct the output directories for SAM files, GFF3 files, and error files
        * append a tuple of arguments to the tasks list for each FASTQ file
        * use ProcessPoolExecutor to process each FASTQ file in parallel by mapping the align_wrapper function to the tasks list
        * collect and store the results of the alignment tasks
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
    """
    Delete the GMAP database directory

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


###########################
### VIsoQLR PROCESSING #### 
###########################

def gff3_data_loading(input_path):
    """
    Load and process a GFF3 file, extracting and returning exon information along with the number of unique reads

    Parameters:
    * input_path (str): Path to GFF3 file to be loaded

    Returns:
    * gff3 (df): A df containing processed GFF3 data with the columns 'gene', 'start', 'end', 'id', and 'score'.
    * num_reads_post_trimming (int): The number of unique reads (based on 'id') after processing

    Workflow:
    1. Load GFF3 File: 
        * reads the GFF3 file, skipping comment lines and loading columns as strings
    2. Filter for Exons: 
        * retains only rows where the third column equals 'exon' and selects relevant columns (gene, start, end, id, score)
    3. Clean and Transform Data:
        * converts 'start' and 'end' columns to integers
        * extracts 'id' using regex
        * converts 'score' column to numeric
    4. Count Unique Reads by calculating number of unique 'id' values
    5. Returns the processed df and the unique read count
    """
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
    """
    Calculates the start and end breakpoints of exons, filters them based on frequency thresholds, and adjusts positions to avoid overlap between adjacent breakpoints

    Parameters:
    * raw_exons (df): df containing exon start and end positions (processed gff3)
    * num_reads_post_trimming (int): The number of unique reads after processing

    Returns:
    * start_position (df): df with start breakpoints and their left and right boundaries after adjustment
    * end_position (df): df with end breakpoints and their left and right boundaries after adjustment

    Workflow:
    1. Initial Count Filtering:
        * counts frequencies of each start and end position
        * discard positions < user-defined 'breakpoint_freq_threshold'
    2. Sorts the start and end breakpoints in ascending order and resets indices for clean referencing
    3. Processing Very Close Breakpoint:
        * if two breakpoints are within 'very_close_bp' and =< 'breakpoint_padding', the breakpoint with lower frequency is removed
    4. Adds left and right boundaries (padding) to each breakpoint
    5. Overlap Resolution:
        * adjusts positions to prevent overlapping of padding boundaries between consecutive breakpoints
        * uses midpoint to split overlapping ranges
    6. Returns cleaned and adjusted start and end dfs with 'breakp', 'left', and 'right' columns.
    """
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
    """
    Assigns breakpoints to exons based on their start and end positions relative to provided breakpoint ranges

    Parameters:
    * raw_exons (df): df containing exon start and end positions (processed gff3)
    * break_points_list (tuple containing 2 dfs):
        - start_position: df of start breakpoints and their corresponding left and right boundaries
        - end_position: df of end breakpoints and their corresponding left and right boundaries

    Returns:
    * bp_assigned_exons (df): df of exons with two new columns, 'start_tag' and 'end_tag', containing assigned breakpoint tags

    Workflow:
    1. Initialise Start and End Tags:
        * adds two new columns, 'start_tag' and 'end_tag', initialised as NaN
    2. Start Breakpoint Assignment:
        * checks if the exon start falls within the breakpoint's left and right boundaries for each start breakpoint in 'start_position' df
        * assigns corresponding 'breakp' value from `start_position` to the 'start_tag' column
    3. End Breakpoint Assignment:
        * checks if the exon end falls within the breakpoint's left and right boundaries for each end breakpoint in 'end_position'
        * assigns the corresponding 'breakp' value from `end_position` to the 'end_tag' column
    4. Fills any remaining NaN values in 'start_tag' and 'end_tag' with `pd.NA` and converts them to integer type (`Int64`)
    5. Returns the updated df `bp_assigned_exons` with 'start_tag' and 'end_tag' columns indicating assigned breakpoints
    """
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
    """
    Calculates detailed breakpoint information, including read frequencies and percentage increases, based on the assigned exons and trimmed reads

    Parameters:
    * break_points_list (tuple containing 2 dfs):
        - start_position_in (df): Start breakpoints with their left and right boundaries
        - end_position_in (df): End breakpoints with their left and right boundaries
    * num_reads_post_trimming (int): The total number of unique reads after trimming
    * bp_assigned_exons (df): df of exons with their start and end tags

    Returns:
    * start_output (df): Summary table for start breakpoints including breakpoints, lower and upper limits, number of reads (with percentages), and total number of reads (with percentages)
    * end_output (df): Summary table for end breakpoints including breakpoints, lower and upper limits, number of reads (with percentages), and total number of reads (with percentages)
    * start_end_position (df): Complete table including both start and end breakpoints with additional frequency and increase percentage columns

    Workflow:
    1. Initialise Start and End Breakpoints:
        * copies `break_points_list` and adds a 'type' column to indicate whether the breakpoint is a start or end
    2. Calculate Frequency of Start and End Breakpoints:
        * computes frequency ('Freq') for both start and end positions by counting occurrences of breakpoints
        * merges the counts with start and end breakpoints to append the frequencies to each breakpoint
    3. Assign Frequency to Start and End Tags:
        * calculates frequency that each breakpoint is tagged in exons for each start and end tag
        - merges the tag counts into respective start and end positions
    4. Concatenates the start and end breakpoint data, calculates the relative frequency (as %), and computes the % increase for intervals between breakpoints
    5. Assigns the gene name from the first exon and include it in result df
    6. Create Formatted Output:
        * formats read counts and % for both start and end positions
        * final output consists of formatted columns for breakpoint, lower/upper limits, and percentages for both the start and end tables 
    7. Returns three dfs:
        * 'start_output': summary of start breakpoints
        * 'end_output': summary of end breakpoints
        * 'start_end_position': detailed table containing both start and end breakpoints along with read frequencies, %, and % increases
    """
    start_position_in = break_points_list[0].copy()
    end_position_in = break_points_list[1].copy()

    start_position_in['type'] = 'start'
    end_position_in['type'] = 'end'

    ### Merge start positions ###
    start_counts = bp_assigned_exons['start'].value_counts().reset_index()
    start_counts.columns = ['breakp', 'Freq']
    start_position = start_position_in.merge(start_counts, on='breakp', how='left').fillna(0)

    ### Merge end positions ###
    end_counts = bp_assigned_exons['end'].value_counts().reset_index()
    end_counts.columns = ['breakp', 'Freq']
    end_position = end_position_in.merge(end_counts, on='breakp', how='left').fillna(0)

    ### Merge start tags ###
    start_tag_counts = bp_assigned_exons['start_tag'].value_counts().reset_index()
    start_tag_counts.columns = ['breakp', 'Freq']
    start_position = start_position.merge(start_tag_counts, on='breakp', how='left', suffixes=('_x', '_y')).fillna(0)

    ### Merge end tags ###
    end_tag_counts = bp_assigned_exons['end_tag'].value_counts().reset_index()
    end_tag_counts.columns = ['breakp', 'Freq']
    end_position = end_position.merge(end_tag_counts, on='breakp', how='left', suffixes=('_x', '_y')).fillna(0)

    ### Combine start and end positions ###
    start_end_position = pd.concat([start_position, end_position], ignore_index=True)
    start_end_position['relative_freq'] = round(start_end_position['Freq_x'] / num_reads_post_trimming * 100, 2)
    start_end_position['relative_freq_padding'] = round(start_end_position['Freq_y'] / num_reads_post_trimming * 100, 2)

    start_end_position['increase_percentage'] = round((start_end_position['Freq_y'] - start_end_position['Freq_x']) / start_end_position['Freq_x'] * 100, 2)
    start_end_position['gene'] = bp_assigned_exons['gene'].iloc[0]
    start_end_position = start_end_position[['type', 'gene', 'breakp', 'left', 'right', 'Freq_x', 'relative_freq', 'Freq_y', 'relative_freq_padding', 'increase_percentage']]
    start_end_position.columns = ['#type', 'gene', 'breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads_breakpoint', 'percent_of_reads_breakpoint', 'number_of_reads_interval', 'percent_of_reads_interval', 'num_of_reads_increase_percentage_interval']

    ### Create formatted read counts ###
    start_end_position['num_of_reads (%)'] = start_end_position['num_of_reads_breakpoint'].astype(str) + ' (' + start_end_position['percent_of_reads_breakpoint'].astype(str) + '%)'
    start_end_position['total_num_of_reads (%)'] = start_end_position['number_of_reads_interval'].astype(str) + ' (' + start_end_position['percent_of_reads_interval'].astype(str) + '%)'

    start_output = start_end_position[start_end_position['#type'] == 'start'][['breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads (%)', 'total_num_of_reads (%)']]
    end_output = start_end_position[start_end_position['#type'] == 'end'][['breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads (%)', 'total_num_of_reads (%)']]

    start_end_position = start_end_position[['#type', 'gene', 'breakpoint', 'lower_limit', 'upper_limit', 'num_of_reads_breakpoint', 'percent_of_reads_breakpoint', 'number_of_reads_interval', 'percent_of_reads_interval', 'num_of_reads_increase_percentage_interval']]

    return start_output, end_output, start_end_position

def exon_definition(bp_assigned_exons):
    """
    Filters and defines valid exons based on breakpoint tags and removes duplicates or erroneous entries

    Parameters:
    * bp_assigned_exons (df): df containing exons with their corresponding start and end breakpoint tags, including exon id, start_tag, and end_tag

    Returns:
    * defined_exons (df): df of exons with valid start and end tags, exons without duplicated tags, exons where the start tag <= end tag, and 'coordinates' column, which combines start and end tags into a single string ('start-end')
    * exon_id_filter (pd.Series): Series of exon IDs that failed any of the filters i.e. invalid exons

    Workflow:
    1. Filter Exons with Missing Values:
        * exon_id_filter1: IDs exons where either 'start_tag' or 'end_tag' = NaN
    2. Filter Exons with Start Greater Than End:
        * exon_id_filter2: IDs exons where the 'start_tag' is > 'end_tag'
    3. Filter Duplicated Exons:
        * dup1: finds exons with duplicated 'id' and 'start_tag' values
        * dup2: finds exons with duplicated 'id' and 'end_tag' values
        * exon_id_filter3: combines the duplicates into a single filter for exons
    4. Combines 'exon_id_filter1', 'exon_id_filter2' and 'exon_id_filter3' to create a full list of exon IDs to filter out invalid exons
    5. Select exons that are excluded by invalid exon filter list
    6. Adds a 'coordinates' column, concatenating the 'start_tag' and 'end_tag' as a string formatted as 'start-end'
    7. Returns two values:
        * defined_exons: df of valid exons including the new 'coordinates' column
        * exon_id_filter: series of exon IDs that were filtered out due to missing values, invalid tags, or duplication
    """
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
    """
    Calculates exon statistics based on defined exons

    Parameters:
    * defined_exons (df): df of exons with 'coordinates' and 'id'

    Returns:
    * exon_info_df (df) containing exon statistics, including:
        - 'exon': start and end positions of the exon ('start-end')
        - 'Size': length of the exon
        - 'number_of_reads': Number of reads of each exon
        - 'percent_of_reads': % of total reads of each exon

    Workflow:
    1. Count the occurrences of each exon coordinate and calculates the % of reads of each exon
    2. Computes the size of each exon by subtracting the start coordinate from the end coordinate.
    3. Sorts and re-indexes the result, keeping columns 'exon', 'Size', 'number_of_reads', and 'percent_of_reads'
    """
    exon_info_df = defined_exons['coordinates'].value_counts().reset_index()
    exon_info_df.columns = ['exon', 'number_of_reads']
    exon_info_df['percent_of_reads'] = round(exon_info_df['number_of_reads'] / defined_exons['id'].nunique() * 100, 1)
    exon_info_df = exon_info_df.sort_values('exon')

    exon_info_df['Size'] = exon_info_df['exon'].str.split('-', expand=True).apply(lambda x: int(float(x[1])) - int(float(x[0])) + 1, axis=1)
    exon_info_df.index = range(1, len(exon_info_df) + 1)
    exon_info_df = exon_info_df[["exon", "Size", "number_of_reads", "percent_of_reads"]]

    return exon_info_df

def isoform_definition(defined_exons):
    """
    Defines isoforms based on exon coordinates and calculates isoform frequencies

    Parameters:
    * defined_exons (df): df of exons with 'id' (read ID) and 'coordinates'

    Returns:
    * isoform_frequencies (df) containing isoform frequency statistics:
        - 'read_isoform': concatenated exon coordinates forming an isoform
        - 'Freq': number of reads of each isoform
        - 'perc': % of total reads of each isoform
    * read_isoform (np.ndarray): array of concatenated exon coordinates for each read (isoform definition)

    Workflow:
    1. Group exons by 'id' and concatenates sorted exon coordinates
    2. Calculate the frequency of each unique isoform and % of reads of each isoform
    """
    read_list = defined_exons.groupby('id')['coordinates'].apply(lambda x: '_'.join(sorted(x)))
    read_isoform = read_list.values

    isoform_frequencies = read_list.value_counts().reset_index()
    isoform_frequencies.columns = ['read_isoform', 'Freq']
    isoform_frequencies['perc'] = round(isoform_frequencies['Freq'] * 100 / isoform_frequencies['Freq'].sum(), 1)

    return isoform_frequencies, read_isoform


def isoform_full_length_filter(isoforms, break_points):
    """
    Filters isoforms that span the full-length region between the first and last breakpoints

    Parameters:
    * isoforms (tuple) containing:
        - isoform_info (df): df with columns 'read_isoform', 'Freq', and 'perc'
        - isoform_read (np.ndarray): array of isoform definitions for each read
    * break_points (tuple) containing start breakpoints df and end breakpoints df

    Returns:
    * keep_isoform_info (df): df of isoforms that span the full-length region (from start to end breakpoints)
    * keep_isoform_read (np.ndarray): array of isoform reads that span the full region
    * partial_reads (np.ndarray): array of isoform reads that do not span the full region

    Workflow:
    1. IDs the minimum start and maximum end breakpoints
    2. Filters isoforms that contain both the start and end breakpoints and 
       splits them into 'keep' (full-length) and 'discard' (partial-length) based on the presence of start and end coordinates
    3. Recalculates the % of reads for isoforms that span the full-length region.
    """
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
    """
    Sorts the segments of an isoform by their starting coordinate

    Parameters:
    * isoform (str): A string of exon segments in the format 'start-end' separated by underscores

    Returns:
    * sorted_isoform (str): The isoform string with segments sorted by their start positions

    Workflow:
    1. Splits the isoform string by underscores into individual segments
    2. Sorts the segments based on the integer value of the start position
    3. Rejoins the sorted segments back into a single string
    """
    segments = isoform.split('_')
    sorted_segments = sorted(segments, key=lambda x: int(x.split('-')[0]))
    return '_'.join(sorted_segments)

def calculate_size(isoform):
    """
    Calculates the total size (number of base pairs) of an isoform

    Parameters:
    * isoform (str): A string of exon segments in the format 'start-end' separated by underscores

    Returns:
    * size (int): The total size of the isoform by summing the lengths of all segments

    Workflow:
    1. Splits the isoform string into segments by underscores
    2. Calculates the size by subtracting the start from the end and adding 1 for each segment
    3. Sums the sizes of all segments to get the total size of the isoform.
    """

    segments = isoform.split('_')
    size = sum([(int(part.split('-')[1]) - int(part.split('-')[0]) + 1) for part in segments])
    return size

def isoform_info_fun(isoform_frequencies, n_inicial, n_post_trim, n_final, n_not_full_length=None):
    """
    Generates information for isoform groups and calculates their frequency and size

    Parameters:
    * isoform_frequencies (df): df of isoform frequencies with columns ['read_isoform', 'Freq', 'perc']
    * n_inicial (int): initial number of reads
    * n_post_trim (int): number of reads after trimming
    * n_final (int): number of reads after filtering
    * n_not_full_length (int, optional): number of partial length reads; default is None

    Returns:
    * all_groups (df) containing detailed isoform information, including:
        - 'Isoform_id': ID for each isoform
        - 'Isoform': sorted isoform string
        - 'Size': total size (bp) of isoform
        - 'Number_of_Reads': number of reads of each isoform
        - 'Percentage': % of reads of each isoform

    Workflow:
    1. Calculate the number of vector reads and reads without consensus breakpoints
    2. Assign unique 'Iso' IDs to each isoform
    3. Calculate the total size by summing the sizes of its segments for each isoform
    4. Concatenate the special read groups (vector, no-consensus, partial) with the isoform frequencies
    5. Format columns and sort segments before returning df
    """
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
    """
    Processes a GFF3 file to extract exon, breakpoint, and isoform data, saving the results

    Parameters:
    * gff3_path (str): path to the GFF3 file
    * sample (str): sample name
    * variant (str): variant name
    * barcode (str): 12bp barcode sequence

    Workflow:
    1. Loads exon data from the {barcode}.gff3 file and calculates the number of reads post-trimming
    2. Calculates breakpoints and assigns them to exons
    3. Defines exons, calculates isoform frequencies, and filters for full-length isoforms
    4. Classifies IDs based on full-length isoforms and generates a dictionary of isoform frequencies
    5. Saves the isoform frequency data to {sample}_{barcode}_iso_freq.csv for the given sample and barcode.

    Returns:
    * Tuple: (sample, variant, barcode, isoform dictionary) if successful, otherwise `None`
    """
    try:
        ### load the GFF3 data ###
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

        ### save iso_dict df as {sample}_{barcode}_iso_freq.csv in visoqlr_dir/{sample}/{variant} dir ###
        output_dir = os.path.join(visoqlr_dir, f"{sample}/{variant}")
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{sample}_{barcode}_iso_freq.csv")
        iso_dict.to_csv(output_path, index=False)

        return (sample, variant, barcode, iso_dict)

    except Exception as e:
        logging.error(f"Error processing {gff3_path}: {e}", exc_info=True)
        return None

def process_and_combine_csv_files(sample, variant, lock):
    """
    Combines {sample}_{barcode}_iso_freq.csv files for each sample variant, ensuring synchronization with a lock

    Parameters:
    * sample (str): sample name
    * variant (str): variant name
    * lock (multiprocessing.Lock): a lock object to ensure only one process performs file combination at a time

    Workflow:
    1. Acquire the lock to synchronize file combination between processes
    2. Find and read all {sample}_{barcode}_iso_freq.csv files into df for the given sample and variant in visoqlr_dir, skipping any empty file
    3. Drops columns 'Isoform_id' and 'Percentage' from each df
    4. Combines all non-empty dfs for each sample variant
    5. Ensures the 'Size' column exists and sums up "Number_of_Reads" for each unique isoform using groupby
    7. Calculates the total number of reads and 'Percentage' of each unique isoform
    8. Adds an 'Isoform_id' column to label each isoform sequentially
    9. Outputs combined df to new {sample}_{variant}_combined.csv file in same dir as the {sample}_{barcode}_iso_freq.csv files it was combined from

    Returns None; saves {sample}_{variant}_combined.csv in visoqlr_dir if valid data is found and processed
    """
    with lock:                                                                                                                  ### ensure only one process can combine files at a time ###
        variant_dir = os.path.join(visoqlr_dir, sample, variant)
        csv_files = glob.glob(os.path.join(variant_dir, f"{sample}_*_iso_freq.csv"))
        dataframes = []
        for file in csv_files:
            df = pd.read_csv(file)

            if df.empty:                                                                                                        ### skip if the df is empty and continue with other files ###
                print(f"Skipping empty file: {file}")
                continue

            for column in ["Isoform_id", "Percentage"]:                                                                         ### check and drop df["Isoform_id"] and df["Percentage"] if they exist ###
                if column in df.columns:
                    df = df.drop(columns=[column])

            dataframes.append(df)

        if dataframes:
            combined_df = pd.concat(dataframes, ignore_index=True)
            
            if 'Size' in combined_df.columns:                                                                                   ### ensure the 'Size' column exists before applying the groupby ###
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
    """
    Processes all GFF3 files in gmap_dir in parallel using multiprocessing to generate isoform data, and combines the CSV files by sample and variant.

    Parameters:
    * gmap_dir (str): Path to the directory containing the GFF3 files organized by sample and variant subdirectories

    Workflow:
    1. Create multiprocessing `Lock` to prevent race conditions during file combination
    2. Use `ProcessPoolExecutor` to process GFF3 files in parallel
    3. Extracts the sample, variant, and barcode information from dir path and submits process_gff3_file for each GFF3 file to executor
    4. Collect results in 'iso_results' after all GFF3 files are processed
    5. Calls process_and_combine_csv_files function to combine the {sample}_{barcode}_iso_freq.csv files for each sample-variant pair once all futures are completed 

    Returns None; saves {sample}_{variant}_combined.csv in visoqlr_dir if valid data is found and processed and logs errors if encountered
    """
    lock = multiprocessing.Lock()                                                                                               ### create a lock for synchronising the combination of CSV files ###
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

            for future in as_completed(futures):                                                                                ###  wait for all iso_dicts to be saved before proceeding ###
                result = future.result()
                if result:
                    iso_results.append(result)

        variants_processed = set((sample, variant) for sample, variant, _, _ in iso_results)                                    ### combine the {sample}_{barcode}_iso_freq.csv files for each variant once, ensuring no race conditions ###
        for sample, variant in variants_processed:
            process_and_combine_csv_files(sample, variant, lock)

    except Exception as e:
        logging.error(f"Error processing GFF3 files in {gmap_dir}: {e}", exc_info=True)

def process_iso_freq_csv_files(visoqlr_dir, visoqlr_no_cdna_dir):
    """
    Removes cDNA isoform from all {sample}_{barcode}_iso_freq.csv files in visoqlr_dir, recalculates % and stores output {sample}_{barcode}_iso_freq.csv in visoqlr_no_cdna_dir
    
    Parameters:
    * visoqlr_dir (str): Path to the directory containing the input {sample}_{barcode}_iso_freq.csv files
    * visoqlr_no_cdna_dir (str): Path to the directory where processed CSV files will be saved

    Workflow:
    1. Iterates through the input visoqlr_dir directory, processing files that match the '_iso_freq.csv' pattern
    2. Loads the {sample}_{barcode}_iso_freq.csv files into dfs, removes columns "Isoform_id" and "Percentage"
    3. Removes rows where 'Isoform' equals the reference start-end positions (i.e. cDNA), and recalculates 'Percentage'
    4. Saves the processed files into visoqlr_no_cdna_dir, maintaining the parent {sample}/{variant} directory structure

    Returns None; saves {sample}_{barcode}_iso_freq.csv in visoqlr_no_cdna_dir if valid data is found and processed and logs errors if encountered
    """
    for root, dirs, files in os.walk(visoqlr_dir):
        for file in files:
            if file.endswith('_iso_freq.csv'):
                csv_path = os.path.join(root, file)                                                                             ### extract file path and sample/barcode info ### 
                try:
                    df = pd.read_csv(csv_path)

                    df = df.drop(columns=['Isoform_id', 'Percentage'])

                    ### remove cDNA (Isoform = 1-1238) ###
                    df = df[df['Isoform'] != f'{ref_start_pos}-{ref_end_pos}']                                                  
                
                    ### recalculate the 'Percentage' column ###
                    total_reads = df['Number_of_Reads'].sum()
                    df['Percentage'] = (df['Number_of_Reads'] / total_reads) * 100
                    df = df.sort_values(by=['Percentage'], ascending=False)

                    ### add back an 'Isoform_id' column with values like 'Iso1', 'Iso2', etc ###
                    df['Isoform_id'] = ['Iso' + str(i + 1) for i in range(len(df))]
                    df = df[['Isoform_id', 'Isoform', 'Size', 'Number_of_Reads', 'Percentage']]

                    ### define the new output path in visoqlr_no_cdna_dir ###
                    relative_path = os.path.relpath(root, visoqlr_dir)
                    output_dir = os.path.join(visoqlr_no_cdna_dir, relative_path)
                    os.makedirs(output_dir, exist_ok=True)
                    output_path = os.path.join(output_dir, file)

                    df.to_csv(output_path, index=False)

                except Exception as e:
                    logging.error(f"Error processing {csv_path}: {e}", exc_info=True)

def process_and_combine_iso_freq_files(visoqlr_no_cdna_dir, visoqlr_combined_bc_dir):
    """
    Combines and process {sample}_{barcode}_iso_freq.csv files belonging to the same visoqlr_no_cdna_dir/{sample}/{variant} dir, 
    and saves the resulting {sample}_{variant}_iso_freq.csv files to visoqlr_combined_bc_dir

    Parameters:
    * visoqlr_no_cdna_dir (str): Path to the directory containing the input {sample}_{barcode}_iso_freq.csv files
    * visoqlr_combined_bc_dir (str): Path to the directory where the combined CSV files will be saved.

    Workflow:
    1. Iterates through the {sample}/{variant} directories in visoqlr_no_cdna_dir to look for {sample}_{variant}_iso_freq.csv files
    2. Loads each CSV file into a df and drops the 'Isoform_id' and 'Percentage' columns
    3. Concatenates dfs belonging to the same visoqlr_no_cdna_dir/{sample}/{variant} dir
    4. Groups the data by 'Isoform', sums 'Number_of_Reads', and takes the first 'Size' for each 'Isoform'
    5. Recalculates the 'Percentage' column based on the total 'Number_of_Reads' and add back 'Isoform_id' column
    6. Saves the combined df as {sample}_{variant}_iso_freq.csv file in visoqlr_combined_bc_dir

    Returns None; saves {sample}_{variant}_iso_freq.csv in visoqlr_combined_bc_dir if valid data is found and processed and logs errors if encountered
    """
    for root, dirs, files in os.walk(visoqlr_no_cdna_dir):
        ### group {sample}_{barcode}_iso_freq.csv by their {sample}/{variant} directories ###
        csv_files = [file for file in files if file.endswith('_iso_freq.csv')]

        if csv_files:
            variant = root.split('/')[-1]
            sample = root.split('/')[-2]

            combined_df = pd.DataFrame()

            ### concat all {sample}_{barcode}_iso_freq.csv in the same {sample}/{variant} directory ###
            for file in csv_files:
                csv_path = os.path.join(root, file)
                try:
                    df = pd.read_csv(csv_path)                                                                                  ### convert each csv to df ###
                    df = df.drop(columns=['Isoform_id', 'Percentage'])                                                          ### drop 'Isoform_id' and 'Percentage' columns ###
                    combined_df = pd.concat([combined_df, df], ignore_index=True)                                               ### concat the current df ###

                except Exception as e:
                    logging.error(f"Error processing {csv_path}: {e}", exc_info=True)

            combined_df = combined_df.groupby('Isoform').agg({                                                                  ### group by 'Isoform' ###
                'Size': 'first',                                                                                                ### take first 'Size' value for each unique 'Isoform' ###
                'Number_of_Reads': 'sum'                                                                                        ### sum 'Number_of_Reads' for each unique 'Isoform' ###
            }).reset_index()

            ### recalculate % ###
            total_reads = combined_df['Number_of_Reads'].sum()
            combined_df['Percentage'] = (combined_df['Number_of_Reads'] / total_reads) * 100
            combined_df = combined_df.sort_values(by=['Percentage'], ascending=False)

            ### add back 'Isoform_id' column (Iso1, Iso2, Iso3, etc.) ###
            combined_df['Isoform_id'] = ['Iso' + str(i + 1) for i in range(len(combined_df))]
            combined_df = combined_df[['Isoform_id', 'Isoform', 'Size', 'Number_of_Reads', 'Percentage']]

            # define the output file path and save the combined df as {sample}_{variant}_iso_freq.csv in visoqlr_combined_bc_dir ###
            output_file = f'{sample}_{variant}_iso_freq.csv'
            output_path = os.path.join(visoqlr_combined_bc_dir, output_file)
            combined_df.to_csv(output_path, index=False)


###################################################################################################################################################################################
                                                                                # VISOQLR PLOTTING #                                                                
###################################################################################################################################################################################

### Set constant variables for isoform plots ###
BAR_HEIGHT = 0.4                                                                                                                ### visual thickness of bars in barh plot ###
FIGURE_HEIGHT_PER_ISOFORM = 0.66                                                                                                ### height per isoform to determine overall figure height ###
BREAKPOINT_PLOT_HEIGHT = 3.5                                                                                                    ### constant figure height for breakpoint plot ###

### Function to parse isoform ###
def parse_isoform(isoform):
    segments = isoform.split('_')
    return [(int(segment.split('-')[0]), int(segment.split('-')[1])) for segment in segments]

### Function to plot isoform data from CSV files ###
def plot_isoform_file(df, ax, filename):
    """
    Visualises isoforms generated per {sample}_{variant} as barh plots

    Parameters:
    * df (df): df containing isoform data with columns 'Isoform', 'Size', 'Number_of_Reads', and 'Percentage'
    * ax (matplotlib.axes.Axes): matplotlib axis where the plot will be rendered
    * filename (str): name of file being plotted to generate plot title

    Workflow:
    1. Parses 'Isoform' column into segments and creates a color map based on the number of unique segments
    2. Iterates over each row, plotting horizontal bars representing the isoform segments and black lines between adjacent segments representing region between exons
    3. Sets y-limits to ensure consistent bar height across plots
    4. Configures axes ticks and tick labels, removes unnecessary spines, add exon legend, derive plot title from filename and other plot configurations
    """
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

    ### set consistent y-limits to ensure equal bar height across all plots ###
    ax.set_ylim(0.5, len(df) + 0.5)

    ### set y-ticks and labels ###
    ax.set_yticks(range(1, len(df) + 1))
    ax.set_yticklabels([f'Isoform {len(df) - i + 1} | {row.Size}bp | {row.Number_of_Reads} | {row.Percentage:.2f}%' for i, row in enumerate(df.itertuples(), 1)], fontsize=8)
    ax.set_xticks([])
    ax.set_xticklabels([])

    ### remove spines ###
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('lightgrey')
    ax.spines['bottom'].set_color('lightgrey')

    ### add legend ###
    handles = [mpatches.Patch(color=color_map[segment], label=f'{segment[0]}-{segment[1]}') for segment in sorted(color_map)]
    ax.legend(handles=handles, title="Exon Coordinates", bbox_to_anchor=(1.02, 1.11), loc='upper left', fontsize=8, title_fontsize=9)

    ### add labels ###
    ax.set_xlim(ref_start_pos - 21, ref_end_pos + 22)
    ax.set_ylabel('Isoform ID | Size | Count | Percentage', labelpad=19, fontsize=10)
    ax.set_title(f"Isoforms of {filename.split('_iso_freq.csv')[0].replace('_', ' ')}", fontsize=12, pad=20)

### Function to plot breakpoints from GFF3 ###
def plot_breakpoints_file(breakpoints_df, ax, sample, variant):
    """
    Plots {sampele}_{variant}'s isoform breakpoint data as vertical bar plot 

    Parameters:
    * breakpoints_df (df): df containing breakpoint data, with columns 'Breakpoints', 'Percentage', and 'Type'
    * ax (matplotlib.axes.Axes): matplotlib axis where the plot will be rendered
    * sample (str): sample name
    * variant (str): variant name

    Workflow:
    * Plots horizontal bars representing the percentage of breakpoints, colored based on the 'Type' ('start' (crimson) or 'end' (mediumblue))
    """
    colors = {'start': 'crimson', 'end': 'mediumblue'}
    ax.bar(breakpoints_df['Breakpoints'], breakpoints_df['Percentage'], color=breakpoints_df['Type'].map(colors), width=1.9)

    ### set axes labels ###
    ax.set_xlabel('Position', fontsize=9)
    ax.set_ylabel('Percentage', fontsize=9)

    ### set x-limits and tick labels ###
    ax.set_xlim(ref_start_pos - 21, ref_end_pos + 22)
    ax.set_ylim(0, 100)
    max_breakpoint = breakpoints_df['Breakpoints'].max()
    ax.set_xticks(range(0, max_breakpoint + 50, 50))
    ax.tick_params(axis='x', labelsize=6)

    ### configure spines ###
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ### plot dots on top of bars ### 
    for idx, row in breakpoints_df.iterrows():
        if row['Percentage'] > 10:
            dot_color = colors[row['Type']]
            ax.plot(row['Breakpoints'], row['Percentage'] + 1, 'o', color=dot_color, markersize=2)

    ax.tick_params(axis='y', labelsize=7)

    ### add legend for the breakpoints ###
    handles = [mpatches.Patch(color='crimson', label='Start'),
               mpatches.Patch(color='mediumblue', label='End')]
    ax.legend(handles=handles, title="Breakpoint Type", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, title_fontsize=9)

### Main function to process files and create combined plots ###
def process_combined_plots(isoform_file, visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir):
    """
    Processes {sample}_{variant}_iso_freq.csv file from visoqlr_combined_bc_dir and corresponding GFF3 files for a specific {sample}/{variant} combination, 
    generating a combined PDF with isoform and breakpoint plots

    Parameters:
    * isoform_file (str): path to {sample}_{variant}_iso_freq.csv file
    8 visoqlr_combined_bc_dir (str): dir of {sample}_{variant}_iso_freq.csv files
    * gmap_dir (str): dir of {barcode}.gff3 files organised by {sample}/{variant}
    * VIsoQLR_plots_dir (str): output directory for saving the pdf files containing isoform and breakpoint plots 
    * sample_variant_gff3_dir (str): output directory for saving processed breakpoint data as {sample}_{variant}.csv files in sample_variant_gff3_dir

    Workflow:
    1. Loads the isoform file, filtering for isoforms with a 'Percentage' >= 1 to be plotted
    2. Dynamically calculates the figure height based on the number of isoforms remmaining
    2. Calls plot_isoform_file function to create subplots, with the first subplot on the top of the pdf to visualise isoforms generated per {sample}_{variant}
    3. Loads all GFF3 files for the given {sample}/{variant}, extracting 'start', 'end', and 'id' fields from exons
    4. Combines the start and end breakpoints and calculate each unique isoform's frequencies and percentages 
    5. Saves the processed breakpoint data as {sample}_{variant}.csv in sample_variant_gff3_dir and calls plot_breakpoints_file function to plot the 2nd subplot, breakpoint plot, on the bottom of pdf file
    6. Outputs the combined isoform and breakpoint plot as {sample}_{variant}.pdf in VIsoQLR_plots_dir
    """
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

        # Pass filename as an argument to plot_isoform_file
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
    """
    Handles parallel processing to generate {sample}_{variant}.pdf containing isoform and breakpoint plots to improve speed and efficiency

    Parameters:
    * visoqlr_combined_bc_dir (str): dir {sample}_{variant}_iso_freq.csv files
    * gmap_dir (str): dir containing {barcode}.gff3} files organised by {sample}/{variant}
    * VIsoQLR_plots_dir (str): output directory for saving {sample}_{variant}.pdf
    * sample_variant_gff3_dir (str): output directory for saving processed breakpoint data as {sample}_{variant}.csv in sample_variant_gff3_dir

    Workflow:
    1. Searches for {sample}_{variant}.csv files in visoqlr_combined_bc_dir and prints the total number of files to process
    2. Utilises `ProcessPoolExecutor` to submit tasks in parallel, where each task processes a single csv file using the 'process_combined_plots' function
    3. Ensures each plot is processed by calling `.result()` on all futures to avoid skipping any tasks
    4. Prints a message if no {sample}_{variant}.csv files are found
    5. Outputs {sample}_{variant}.pdf containing isoform and breakpoint plots for each {sample}_{variant} in VIsoQLR_plots_dir
    """
    isoform_files = [os.path.join(visoqlr_combined_bc_dir, file) for file in os.listdir(visoqlr_combined_bc_dir) if file.endswith('csv')]

    if isoform_files:
        print(f"Found {len(isoform_files)} CSV files to process.")

        with ProcessPoolExecutor(max_workers=max_workers) as executor: 
            futures = [executor.submit(process_combined_plots, isoform_file, visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir) for isoform_file in isoform_files]
            for future in futures:
                future.result()  
    else:
        print(f"No CSV files found in directory: {visoqlr_combined_bc_dir}")

### Function to plot {sample} breakpoint plots ###
def plot_sample_breakpoints(gmap_dir, sample_gff3_dir, sample_plot_dir):
    """
    Processes {barcode}.gff3 files in {sample}/{variant} directories to generate {sample} breakpoint frequency plots

    Parameters:
    * gmap_dir (str): dir containing {barcode}.gff3} files organised by {sample}/{variant}
    * sample_gff3_dir (str): dir to save the aggregated breakpoint data as {sample}.csv files in sample_gff3_dir
    * sample_plot_dir (str): dir to save the breakpoint frequency plots as {sample}.html in sample_plot_dir

    Workflow:
    1. Walks through gmap_dir to find {barcode}.gff3 files in {sample}/{variant} subdirectories
    2. Extracts exon coordinates (start and end) from each {barcode}.gff3 file
    3. Aggregates data by 'start' and 'end' coordinates to calculate the frequency and % of each unique breakpoint
    4. Saves the aggregated data as {sample}.csv files in sample_gff3_dir
    5. Plots the breakpoint frequencies using Plotly, where start breakpoints are crimson bars, end breakpoints are mediumblue bars and dots are added atop bars if % of breakpoint > 10
    6. Configure hover labels and legend
    7. Saves the interactive plot as {sample}.html sample_plot_dir
    """
    ### create dictionary to group GFF3 files by sample ###
    sample_data = {}

    ### walk through each {sample}/{variant} directory in gmap_dir ###
    for root, dirs, files in os.walk(gmap_dir):
        gff3_files = glob.glob(os.path.join(root, "*.gff3"))

        if gff3_files:
            ### extract sample name ###
            sample = root.split(os.sep)[-3]

            ### process each {barcode}.gff3 file in current dir ###
            for gff3_file in gff3_files:
                ### check if file is empty before processing and continue if yes ###
                if os.stat(gff3_file).st_size == 0:
                    print(f"Skipping empty file: {gff3_file}")
                    continue

                ### process each {barcode}.gff3 file in current dir ###
                try:
                    gff3 = pd.read_csv(gff3_file, sep='\t', header=None, comment='#', dtype=str)
                    gff3.index = range(1, len(gff3) + 1)
                    
                    ### filter for exon info and select relevant columns ###
                    gff3 = gff3[gff3[2] == 'exon'][[3, 4, 8]]
                    gff3.columns = ['start', 'end', 'id']
                    
                    ### extract exon ID and convert start and end to integers ###
                    gff3['id'] = gff3['id'].str.extract(r'Name=([^;]+)')
                    gff3['start'] = gff3['start'].astype(int)
                    gff3['end'] = gff3['end'].astype(int)

                    ### append data to the corresponding sample in the dictionary ###
                    if sample not in sample_data:
                        sample_data[sample] = []
                    sample_data[sample].append(gff3)

                except pd.errors.EmptyDataError:
                    print(f"Skipping {gff3_file}: No columns to parse.")
                    continue

    ### concatenate the data and perform aggregation for each sample ###
    for sample, dataframes in sample_data.items():
        combined_df = pd.concat(dataframes, ignore_index=True)
        number_of_unique_reads = combined_df['id'].nunique()
        combined_df = combined_df[["start", "end"]]

        ### reshape data to long format i.e. stack up 'start' and 'end' values in new column 'Breakpoints' ###
        starts = combined_df[['start']].rename(columns={'start': 'Breakpoints'})
        starts['Type'] = 'start'

        ends = combined_df[['end']].rename(columns={'end': 'Breakpoints'})
        ends['Type'] = 'end'

        ### combine start and end data into a single df ###
        breakpoints_df = pd.concat([starts, ends], ignore_index=True)

        ### calculate frequency of each unique breakpoint and type combo and its percentage ###
        breakpoints_df = breakpoints_df.groupby(['Breakpoints', 'Type']).size().reset_index(name='Frequency')
        breakpoints_df['Percentage'] = (breakpoints_df['Frequency'] / number_of_unique_reads) * 100
        
        ### define the output file path and save as sample}.csv to sample_gff3_dir ###
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

            ### add dots on top of bars where "Percentage" of breakpoints > 10 ###
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
        
        ### configure the x-axis limits and ticks and customise layout ###
        x_min = breakpoints_df["Breakpoints"].min() - 19
        x_max = breakpoints_df["Breakpoints"].max() + 19
        x_ticks = list(range(0, x_max + 50, 50))  

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

        ### save the plot as {sample}.html file in sample_plot_dir
        html_output_file = os.path.join(sample_plot_dir, f"{sample}.html")
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

    plot_start_time = time.time()
    plot_iso_bp_parallel(visoqlr_combined_bc_dir, gmap_dir, VIsoQLR_plots_dir, sample_variant_gff3_dir)
    plot_end_time = time.time()
    print(f"Time taken to generate plot of isoforms per variant: {plot_end_time - plot_start_time:.3f} seconds")

    delete_gmap_database(database_dir, database_name)
    
    remove_codon_file(root_output_dir)

if __name__ == "__main__":
    main()
