import os
import gzip
import re
import subprocess
import pandas as pd
import glob
import cstag
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed 
import pysam
import csv
from itertools import combinations
from collections import Counter
import warnings 
import matplotlib.pyplot as plt 
import time     
import numpy as np
import shutil
import yaml
from __future__ import annotations
from pathlib import Path

# silence future warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

### Determine the full path of which pipeline is in ###
vesper_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(vesper_dir, "vesper_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from vesper_config.yaml ###
tile_number = config["tile_number"]

baseq_threshold = config["baseq_threshold"]
mapq_threshold = config["mapq_threshold"]
read_length_filter = config["read_length_filter"]

### T1 start and end ###
start_bp = config["start_bp"]
end_bp = config["end_bp"]

cDNA_file = config["cDNA_file"]
raw_fastqgz_dir = config["raw_fastqgz_dir"]
ref_fasta = config["ref_fasta"]
root_output_dir = config["root_output_dir"]
os.makedirs(root_output_dir, exist_ok=True)
trim_len = config["trim_length"]

hiseq_fastq_dir = os.path.dirname(os.path.dirname(raw_fastqgz_dir))
trimmed_fastq_dir = os.path.join(hiseq_fastq_dir, "trimmed_unzipped", f"T{tile_number}")

minimap2_sam_dir = os.path.join(root_output_dir, "minimap2_sam")
filtered_sam_dir = os.path.join(root_output_dir, "filtered_sam", "unzipped")
filtered_samgz_dir = os.path.join(root_output_dir, "filtered_sam", "gzipped")

intermediate_results_output_dir = os.path.join(root_output_dir, 'intermediate_files')
main_results_output_dir = os.path.join(root_output_dir, 'main_results')
raw_output_dir = os.path.join(intermediate_results_output_dir, 'raw_variants_per_sample')

variants_per_sample_dir = os.path.join(main_results_output_dir, 'variants_per_sample')
snv_per_sample_dir = os.path.join(main_results_output_dir, 'SNV_per_sample')
nm_counts_dir = os.path.join(main_results_output_dir, 'nm_counts')
nm_histograms_dir = os.path.join(main_results_output_dir, 'nm_histograms')

for dir_path in [minimap2_sam_dir, filtered_sam_dir, filtered_samgz_dir, intermediate_results_output_dir, main_results_output_dir, 
                 raw_output_dir, variants_per_sample_dir, snv_per_sample_dir, nm_counts_dir, nm_histograms_dir]:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)

max_workers = os.cpu_count()
print(f"executing pipeline with system-detected number of threads: {max_workers}")


###############################################################################################################################
                ### UNZIP *.FASTQ.GZ FILES AND TRIM OFF PRIMER REGION (20 BP FROM 5' & 3' END OF EACH READ) ###
###############################################################################################################################

def trim_and_unzip_file(args):
    raw_fastqgz_path, trimmed_fastq, trim_len = args
    try:
        with gzip.open(raw_fastqgz_path, 'rt') as f_in, open(trimmed_fastq, 'w') as f_out:
            while True:
                header = f_in.readline()
                if not header:
                    break  # EOF
                seq = f_in.readline().strip()
                plus = f_in.readline()
                qual = f_in.readline().strip()

                # Trim both 5' and 3' ends
                trimmed_seq = seq[trim_len:-trim_len] if len(seq) > 2 * trim_len else ''
                trimmed_qual = qual[trim_len:-trim_len] if len(qual) > 2 * trim_len else ''

                if trimmed_seq and trimmed_qual:
                    f_out.write(f"{header}{trimmed_seq}\n{plus}{trimmed_qual}\n")

        return f"✔ Trimmed and saved: {os.path.basename(trimmed_fastq)}"
    except Exception as e:
        return f"✘ Error with {os.path.basename(raw_fastqgz_path)}: {e}"

def parallel_trim_and_unzip(raw_fastqgz_dir, trimmed_fastq_dir):
    os.makedirs(trimmed_fastq_dir, exist_ok=True)

    tasks = []
    for filename in os.listdir(raw_fastqgz_dir):
        if filename.endswith('.fastq.gz') or filename.endswith('.fq.gz'):
            raw_fastqgz_path = os.path.join(raw_fastqgz_dir, filename)
            output_filename = filename.replace('.fastq.gz', '.fastq').replace('.fq.gz', '.fastq')
            trimmed_fastq = os.path.join(trimmed_fastq_dir, output_filename)
            tasks.append((raw_fastqgz_path, trimmed_fastq, trim_len))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(trim_and_unzip_file, task) for task in tasks]
        for future in as_completed(futures):
            print(future.result())


###############################################################################################################################
                                ### PERFORM ALIGNMENT OF TRIMMED *.FASTQ FILES WITH MINIMAP2 ###                                
###############################################################################################################################

def run_minimap2(trimmed_fastq_dir):
    """
    Runs minimap2 on paired FASTQ files in the given directory.

    Parameters:
        trimmed_fastq_dir (str): Directory containing *_1.fastq and *_2.fastq files.

    Output:
        Creates a SAM file for each paired FASTQ input.
    """
    # Regex pattern to parse sample, rep, tile, read
    pattern = re.compile(r"(?P<sample>.+?)(?P<rep>\d+)_(?P<tile>T\d)_(?P<read>[12])\.fastq$")

    # Get all R1 files (ending in _1.fastq)
    fastq_files = [f for f in os.listdir(trimmed_fastq_dir) if f.endswith("_1.fastq")]

    for r1_filename in fastq_files:
        match = pattern.match(r1_filename)
        if not match:
            print(f"Skipping unrecognized file: {r1_filename}")
            continue

        sample = match.group("sample")
        rep = match.group("rep")
        tile = match.group("tile")

        r1_path = os.path.join(trimmed_fastq_dir, r1_filename)
        r2_filename = r1_filename.replace("_1.fastq", "_2.fastq")
        r2_path = os.path.join(trimmed_fastq_dir, r2_filename)

        if not os.path.exists(r2_path):
            print(f"Missing R2 for {r1_filename}, skipping.")
            continue

        # Construct output SAM filename
        minimap2_sam_filename = f"{rep}_{sample}_{tile}.sam"
        minimap2_sam_path = os.path.join(minimap2_sam_dir, minimap2_sam_filename)

        # Build and run minimap2 command
        cmd = [
            "minimap2", "-ax", "sr", "-t", str(max_workers),
            "--cs", "--sam-hit-only", "--secondary=no",
            ref_fasta, r1_path, r2_path
        ]

        with open(minimap2_sam_path, "w") as sam_out:
            subprocess.run(cmd, stdout=sam_out, check=True)


###############################################################################################################################
                ### FILTER SAM FILES TO RETAIN ONLY PRIMARY ALIGNMENTS > MAPQ_THRESHOLD & READ_LENGTH_FILTER ###                                
###############################################################################################################################

def filter_minimap2_sam(minimap2_sam_path, filtered_sam_path):
    """Filters SAM file based on read length and mapping quality."""
    print(f"Filtering: {minimap2_sam_path}")

    try:
        with pysam.AlignmentFile(minimap2_sam_path, "r", check_sq=False) as infile, open(filtered_sam_path, "w") as outfile:
            # Write headers
            for header_line in infile.text.splitlines():
                outfile.write(header_line + "\n")

            # Apply filtering
            for read in infile.fetch(until_eof=True):
                if (
                    not read.is_secondary and
                    not read.is_supplementary and
                    not read.is_unmapped and
                    (len(read.query_sequence) == read_length_filter if read_length_filter else True) and
                    read.mapping_quality >= mapq_threshold
                ):
                    outfile.write(read.to_string() + "\n")
    except Exception as e:
        print(f"Error processing {minimap2_sam_path}: {e}")

def compress_filtered_sam(filtered_sam_path, filtered_samgz_path):
    """Compress a SAM file into .gz format."""
    print(f"Compressing: {filtered_sam_path}")
    with open(filtered_sam_path, 'rb') as f_in, gzip.open(filtered_samgz_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def process_single_sam(file_name):
    """Processes a single SAM or gzipped SAM file."""
    minimap2_sam_path = os.path.join(minimap2_sam_dir, file_name)

    filtered_sam_path = os.path.join(filtered_sam_dir, os.path.basename(minimap2_sam_path))
    filtered_samgz_path = os.path.join(filtered_samgz_dir, os.path.basename(minimap2_sam_path) + ".gz")

    # Step 1: Filter SAM
    filter_minimap2_sam(minimap2_sam_path, filtered_sam_path)
    # Step 2: Compress filtered file
    compress_filtered_sam(filtered_sam_path, filtered_samgz_path)
    print(f"Finished processing {file_name}")

def process_sam_files():
    """Processes SAM files in parallel using ProcessPoolExecutor."""
    files = [f for f in os.listdir(minimap2_sam_dir) if f.endswith('.sam')]

    #max_workers = min(len(files), os.cpu_count())  # Limit workers to available CPU cores

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_sam, file_name): file_name for file_name in files}

        for future in futures:
            try:
                future.result()  # Raise exceptions if any occur during processing
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")

###############################################################################################################################
                            ### GET READ COUNTS TO CALCULATE VARIANT FREQUENCIES ###                                
###############################################################################################################################

def get_read_counts(filtered_samgz_path):
    """Process a single SAM file and calculate read quality scores without using large DataFrames."""
    print(f"Processing Read Counts for: {filtered_samgz_path}")

    position_counts = Counter()
    total_scores = 0
    num_positions = 0

    try:
        with pysam.AlignmentFile(filtered_samgz_path, "r", check_sq=False) as samgz_file:
            for read in samgz_file:
                start_position = read.reference_start
                if read.query_qualities is not None:
                    for i, quality in enumerate(read.query_qualities):
                        position_value = start_position + i
                        if quality >= baseq_threshold:
                            position_counts[position_value] += 1

        # Only keep values within the range
        position_counts = {pos: count for pos, count in position_counts.items() if start_bp <= pos <= end_bp}

        # Compute average
        total_scores = sum(position_counts.values())
        num_positions = len(position_counts)

        avg_score = total_scores / (end_bp - start_bp + 1) if num_positions > 0 else 0
        avg_score = round(avg_score)

        return os.path.basename(filtered_samgz_path).replace('.sam.gz', ''), avg_score

    except Exception as e:
        print(f"Error processing {filtered_samgz_path}: {e}")
        return os.path.basename(filtered_samgz_path).replace('.sam.gz', ''), None

def get_read_counts_parallel(filtered_samgz_dir):
    """Process SAM files in parallel using ProcessPoolExecutor."""
    sam_gz_files = [os.path.join(filtered_samgz_dir, f) for f in os.listdir(filtered_samgz_dir) if f.endswith('.sam.gz')]

    read_counts = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(get_read_counts, filtered_samgz_path): filtered_samgz_path for filtered_samgz_path in sam_gz_files}

        for future in futures:
            try:
                filename, avg_score = future.result()
                read_counts[filename] = avg_score
                print(f"Processed {filename}: {avg_score}")
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")

    return read_counts


###############################################################################################################################
                                            ### MAP VARIANTS TO ITS CDNA POSITION ###                                
###############################################################################################################################

# Create dictionary to map bp to cDNA
def map_pos_to_cdna(cDNA_file):
    cdna_df = pd.read_csv(cDNA_file)
    cnda_df = cdna_df[["bp", "cDNA"]]
    cdna_df['cDNA'] = cnda_df['cDNA'].astype(str)
    cdna_df = cdna_df.drop_duplicates(subset=['bp']).reset_index(drop=True)
    bp_to_cdna = dict(zip(cdna_df['bp'], cdna_df['cDNA']))
    return bp_to_cdna
bp_to_cdna = map_pos_to_cdna(cDNA_file)

# Function to convert each variant in 'cs_long_sub_baseq30' to 'cdna_variants'
def convert_to_cdna(variants):
    cdna_variants = []
    for variant in variants.split(', '):
        # Split position and mutation part (e.g., "168065*ag" -> "168065" and "ag")
        try:
            pos, mutation = variant.split('*')
            pos = int(pos)  # Convert position to an integer to look up in bp_to_cdna
            cdna_pos = bp_to_cdna.get(pos, pos)  # Use cDNA position if available; otherwise keep original
            cdna_variants.append(f"{cdna_pos}*{mutation}")
        except ValueError:
            # Handle any unexpected format errors
            cdna_variants.append(variant)  # Keep variant as is if it cannot be processed
    return ', '.join(cdna_variants)

# Function to reformat each variant in 'cdna_variants'
def reformat_variant(variants):
    formatted_variants = []
    for variant in variants.split(', '):
        try:
            cdna_position, bases = variant.split('*')
            ref_base, obs_base = bases[0].upper(), bases[1].upper()  # Capitalize the bases
            formatted_variants.append(f"{cdna_position}{ref_base}>{obs_base}")
        except ValueError:
            # If there's a format error, skip this variant
            continue
    # Join the formatted variants back into a string
    return ', '.join(formatted_variants)


###############################################################################################################################
                                    ### CALCULATE SNVS FROM CS SHORT AND LONG TAGS ###                                
###############################################################################################################################
def open_sam_file(sam_file_path):
    """Open a SAM file, supporting both .sam and .sam.gz formats."""
    if sam_file_path.endswith('.gz'):
        return gzip.open(sam_file_path, 'rt')
    else:
        return open(sam_file_path, 'r')

def calculate_cs_short_sub(row):
    cs_tag = row['minimap2_cs']
    current_position = row['read_start']
    substitutions = []

    # Regular expression to match each component in cs_tag
    pattern = re.compile(r'(\d+)|(\*[a-z]{2})|(\+[a-z]+)|(-[a-z]+)')

    for match in pattern.finditer(cs_tag):
        if match.group(1):                                                                                          # Exact matches (e.g., ":50")
            current_position += int(match.group(1))
        elif match.group(2):                                                                                        # Substitutions (e.g., "*ct")
            sub_base = match.group(2)[1:]                                                                           # Extract 'ct' after '*'
            substitutions.append(f"{current_position}*{sub_base}")
            current_position += 1                                                                                   # Move to the next position after substitution
        elif match.group(3):                                                                                        # Insertions (e.g., "+at")
            pass                                                                                                    # Ignore insertion bases for position calculation
        elif match.group(4):                                                                                        # Deletions (e.g., "-gtac")
            del_length = len(match.group(4)) - 1                                                                    # Subtract 1 to ignore the '-' sign
            current_position += del_length                                                                          # Account for deleted bases

    # Join all substitutions for the row, separated by commas
    return ', '.join(substitutions)

def calculate_cs_long_sub(row):
    cs_masked = row['cs_masked']
    current_position = row['read_start']
    substitutions = []

    # Regular expression to match each component in cs_masked
    pattern = re.compile(r'(=[A-Za-z]+)|(\*[a-z]{2})|(\+[a-z]+)|(-[a-z]+)')

    for match in pattern.finditer(cs_masked):
        if match.group(1):                                                                                          # Matches (e.g., "=CGATCG")
            # Extract matched segment and increment position for each base except '='
            match_segment = match.group(1)[1:]                                                                      # Ignore '='
            current_position += len(match_segment)
        elif match.group(2):                                                                                        # Substitutions (e.g., "*ct")
            sub_base = match.group(2)[1:]                                                                           # Extract 'ct' after '*'
            substitutions.append(f"{current_position}*{sub_base}")
            current_position += 1                                                                                   # Move to the next position after substitution
        elif match.group(3):                                                                                        # Insertions (e.g., "+gtc")
            pass                                                                                                    # Ignore insertion bases for position calculation
        elif match.group(4):                                                                                        # Deletions (e.g., "-ata")
            del_length = len(match.group(4)) - 1                                                                    # Subtract 1 to ignore the '-' sign
            current_position += del_length                                                                          # Account for deleted bases

    # Join all substitutions for the row, separated by commas
    return ', '.join(substitutions)


###############################################################################################################################
                                    ### REMOVE VARIANTS THAT ARE LESS THAN BASEQ_THRESHOLD ###                                
###############################################################################################################################

# Function to remove variants with baseq < baseq_threshold
def remove_n_variants(variants):
    # Split by comma and filter out any variant containing 'n'
    filtered_variants = [variant for variant in variants.split(', ') if 'n' not in variant]
    # Join the filtered list back into a string
    return ', '.join(filtered_variants)

# Function to filter variants based on position range
def filter_variants(variants):
    filtered_variants = []
    for variant in variants.split(', '):
        # Split the position and mutation part
        try:
            pos, mutation = variant.split('*')
            pos = int(pos)  # Convert position to an integer
            # Keep only if position is within the specified range
            if start_bp <= pos <= end_bp:
                filtered_variants.append(f"{pos}*{mutation}")
        except ValueError:
            # If there's a format error, skip this variant
            continue
    # Join the filtered list back into a string
    return ', '.join(filtered_variants)


###############################################################################################################################
                            ### EXTRACT FEATURES FROM FILTERED *.SAM.GZ FILE AND PARSE THEM INTO DF ###                                
###############################################################################################################################

def sam_to_dataframe(sam_file_path):
    """
    Convert a SAM file into a DataFrame with columns:
    'read_id', 'flag', 'read_start', 'seq', 'base_qualities', 'minimap2_cs', 'cigar', "minimap2_cs".
    """
    data = []

    with pysam.AlignmentFile(sam_file_path, "rb") as samfile:
        for read in samfile:
            # Extract the required fields
            read_id = read.query_name
            flag = read.flag
            read_start = read.reference_start + 1  # SAM is 1-based
            seq = read.query_sequence
            base_qualities = read.qual
            cigar = read.cigarstring
            length = read.query_alignment_length

            # Extract 'CS' tags
            minimap2_cs = None

            for tag, value in read.get_tags():
                if tag == "cs":
                    minimap2_cs = value
                    break

            data.append([read_id, flag, read_start, seq, base_qualities, minimap2_cs, cigar, length])

    # Create a DataFrame from the list of data
    df = pd.DataFrame(data, columns=['read_id', 'flag', 'read_start', 'seq', 'base_qualities', 'minimap2_cs', 'cigar', 'length'])

    # Sort by read_id
    df.sort_values(by='read_id', inplace=True)

    # Apply cstag.lengthen to create a new column 'cs_long'
    df['cs_long'] = df.apply(
        lambda row: cstag.lengthen(row['minimap2_cs'], row['cigar'], row['seq'])
        if row['minimap2_cs'] and row['cigar'] and row['seq'] else None,
        axis=1
    )

    # Create a new DataFrame (nm0_df) with rows where 'cs_long' does not contain '*'
    nm0_df = df[~df['cs_long'].str.contains(r'\*', na=False)].copy()

    # Remove rows where 'cs_tag' contains no '*'
    df = df[df['cs_long'].str.contains(r'\*')]

    # Count the number of unique read_id values in nm0_df
    removed_row_count = nm0_df['read_id'].nunique()

    # Reset the DataFrame index
    df.reset_index(drop=True, inplace=True)

    # Apply cstag.mask to create a new column 'cs_masked'
    df['cs_masked'] = df.apply(
        lambda row: cstag.mask(row['cs_long'], row['cigar'], row['base_qualities'], baseq_threshold)
        if row['cs_long'] and row['cigar'] and row['base_qualities'] else None,
        axis=1
    )

    df = df[["read_id", "flag", "read_start","minimap2_cs", "cs_masked"]]
    df['cs_short_sub'] = df.apply(calculate_cs_short_sub, axis=1)                                                   ### get variants from cs_short tags ###
    df['cs_long_sub'] = df.apply(calculate_cs_long_sub, axis=1)                                                     ### get variants from cs_long tags ###
    df['cs_long_sub_baseq30'] = df['cs_long_sub'].apply(remove_n_variants)                                          ### remove variants with baseq < 30 ###                   
    df['cs_long_sub_baseq30'] = df['cs_long_sub_baseq30'].apply(filter_variants)                                    ### remove variants outside of position range ###

    # Create a copy of the DataFrame for processing 'cs_long_sub'
    processed_df = df[["read_id", "flag", "cs_long_sub_baseq30"]].copy()

    # Remove rows where 'cs_long_sub' is NaN or an empty string after stripping whitespace
    processed_df['cs_long_sub_baseq30'] = processed_df['cs_long_sub_baseq30'].fillna('').str.strip()
    processed_df = processed_df[processed_df['cs_long_sub_baseq30'] != ""]

    processed_df['cdna_variants'] = processed_df['cs_long_sub_baseq30'].apply(convert_to_cdna)          ### convert to cDNA variants ###
    processed_df['cdna_variants'] = processed_df['cdna_variants'].apply(reformat_variant)               ### reformat cdna variants ###  

    ################ FIND COMMON VARIANTS SHARING THE SAME READ_ID ################
    # Drop 'flag' and 'cs_long_sub_baseq30' columns
    processed_df = processed_df[['read_id', 'cdna_variants']]

    # Define a function to find common variants 
    def find_common_variants(variants_list):
        # Split each string of variants by ", ", convert to a set, and find the intersection
        sets = [set(variants.split(', ')) for variants in variants_list]
        common_variants = set.intersection(*sets) if sets else set()
        # Join the common variants back into a single comma-separated string
        return ', '.join(sorted(common_variants))

    processed_df_grouped = processed_df.groupby('read_id').agg({
        'cdna_variants': find_common_variants
    }).reset_index()

    # Save the processed DataFrame to a separate CSV file
    raw_output_file = os.path.basename(sam_file_path).replace('.sam.gz', '.csv')
    raw_output_path = os.path.join(raw_output_dir, raw_output_file)
    processed_df_grouped.to_csv(raw_output_path, index=False)

    return processed_df_grouped, removed_row_count

def process_single_sam_file(sam_file):
    """
    Processes a single SAM file and returns the filename, processed DataFrame, and removed row count.
    """
    filename = os.path.basename(sam_file).split('.')[0]
    processed_df_grouped, removed_row_count = sam_to_dataframe(sam_file)
    return filename, processed_df_grouped, removed_row_count

def process_sam_cs_variants(filtered_samgz_dir):
    """
    Process each .sam.gz file in a directory, convert to a grouped DataFrame, and collect removed row counts.
    """
    sam_files = [os.path.join(filtered_samgz_dir, f) for f in os.listdir(filtered_samgz_dir) if f.endswith('.sam.gz')]
    dfs = []
    removed_counts = {}

    # Use ProcessPoolExecutor to run file processing in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_sam_file, sam_file): sam_file for sam_file in sam_files}

        # Collect results as they complete
        for future in as_completed(futures):
            sam_file = futures[future]
            try:
                filename, processed_df_grouped, removed_row_count = future.result()
                dfs.append(processed_df_grouped)
                removed_counts[filename] = removed_row_count
            except Exception as exc:
                print(f"{sam_file} generated an exception: {exc}")

    return dfs, removed_counts


###############################################################################################################################
                                    ### GET COUNTS OF EACH UNIQUE CDNA VARIANT COMBINATION ###                                
###############################################################################################################################
def get_variants(file_path, variants_per_sample_dir, snv_per_sample_dir, nm_histograms_dir, removed_counts, read_counts_dict):
    """
    Process each CSV file in raw_output_dir by dropping all columns except 'cdna_variants',
    grouping by unique 'cdna_variants' values, and counting occurrences.

    Parameters:
    - file_path: File path to each CSV file in raw_output_dir
    - variants_per_sample_dir: Path to the directory where processed files with all cDNA_variants per read will be saved.
    - snv_per_sample_dir: Path to directory where total counts for each SNV per sample will be saved
    - removed_counts: Dictionary containing removed row counts with sample IDs as keys.
    - nm_histograms_dir: Path to the directory where NM histograms will be saved as PDF files.
    """
    filename = os.path.basename(file_path)
    sample_ID = filename.split('.csv')[0]

    # Read the CSV into a DataFrame and keep only 'cdna_variants'
    df = pd.read_csv(file_path, usecols=["cdna_variants"])

    # Group by unique 'cdna_variants' and count occurrences
    processed_df = df.groupby("cdna_variants").size().reset_index(name='count')

    # Sort the DataFrame by count in descending order
    processed_df = processed_df.sort_values(by='count', ascending=False)

    # Rename the count column to include the sample_ex_replicate number
    processed_df.rename(columns={'count': f'count_{sample_ID}'}, inplace=True)

    # Add depth to df
    processed_df[f'depth_{sample_ID}'] = read_counts_dict[sample_ID] if sample_ID in read_counts_dict else print(f"{sample_ID} not found in read_counts_dict")

    # Save the processed DataFrame to the output directory
    output_file_path = os.path.join(variants_per_sample_dir, filename)
    processed_df.to_csv(output_file_path, index=False)

    # Calculate SNV counts
    explode_df = processed_df.assign(cdna_SNV=processed_df["cdna_variants"].str.split(", ")).explode("cdna_SNV")
    snv_df = explode_df.groupby("cdna_SNV", as_index=False).agg({f"count_{sample_ID}": "sum", f"depth_{sample_ID}": "first"})
    snv_df["pos"] = snv_df["cdna_SNV"].str.extract(r"([-]?\d+)").astype(int)
    snv_df = snv_df.sort_values(by="pos").drop(columns=["pos"]).reset_index(drop=True)
    snv_output_file = os.path.join(snv_per_sample_dir, f"{sample_ID}_SNV.csv")
    snv_df.to_csv(snv_output_file, index=False)

    # Calculate number of mutations (NM) for each cdna_variant
    processed_df['NM'] = processed_df['cdna_variants'].apply(lambda x: 1 + x.count(','))
    processed_df = processed_df.drop(columns=['cdna_variants', f'depth_{sample_ID}'])

    # Get counts of unique NMs per sample_rep
    grouped_df = processed_df.groupby('NM', as_index=False).agg({f'count_{sample_ID}': 'sum'})

    # Use removed_counts to add a row to group_df where NM=0
    nm_zero_count = removed_counts.get(sample_ID, 0)
    nm_zero_row = pd.DataFrame({'NM': [0], f'count_{sample_ID}': [nm_zero_count]})
    grouped_df = pd.concat([grouped_df, nm_zero_row], ignore_index=True)
    nm_counts_file_path = os.path.join(nm_counts_dir, f"{sample_ID}_NM_counts.csv")
    grouped_df.to_csv(nm_counts_file_path, index=False)

    # Plot histogram to visualise NM distribution per sample_rep
    plt.figure(figsize=(10, 6))
    plt.bar(grouped_df['NM'], grouped_df[f'count_{sample_ID}'], width=0.8, align='center', color='#DDD5F3')

    ### Dynamically set tick intervals ###
    max_count = grouped_df[f'count_{sample_ID}'].max()
    tick_interval = 10 ** (len(str(max_count)) - 2)

    # Set y limits
    plt.ylim(0, max_count*1.1)
    plt.xlabel("Number of Mutations per Variant", labelpad=10, fontsize=10)
    plt.ylabel("Count", labelpad=10, fontsize=10)
    plt.title(f"Distribution of Number of Mutations per Variant for {sample_ID}", fontsize=13, pad=13)
    plt.xticks(grouped_df['NM'], fontsize=8)
    plt.xticks(range(0, grouped_df['NM'].max() + 1, 1), fontsize=8)
    plt.yticks(range(0, max_count + tick_interval, tick_interval), fontsize=8)

    # plot grid behind bars
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().set_axisbelow(True)

    # Save the histogram as a PDF
    histogram_path = os.path.join(nm_histograms_dir, f"{sample_ID}_NM_histogram.pdf")
    plt.savefig(histogram_path, format='pdf', bbox_inches='tight')
    plt.close()

def get_variants_parallel(raw_output_dir, variants_per_sample_dir, snv_per_sample_dir, nm_histograms_dir, removed_counts, read_counts_dict):
    """
    Process CSV files in parallel using ProcessPoolExecutor.
    """
    csv_files = [os.path.join(raw_output_dir, f) for f in os.listdir(raw_output_dir) if f.endswith('.csv')]

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(get_variants, file_path, variants_per_sample_dir, snv_per_sample_dir, nm_histograms_dir, removed_counts, read_counts_dict): file_path for file_path in csv_files}

        for future in as_completed(futures):
            try:
                filename = future.result()
                # print(f"Processed {filename}")
            except Exception as e:
                print(f"Exception processing {futures[future]}: {e}")


def main():
    start_time = time.time()

    trim_and_unzip_fq_start_time = time.time()
    parallel_trim_and_unzip(raw_fastqgz_dir, trimmed_fastq_dir)
    trim_and_unzip_fq_end_time = time.time()
    print(f"Time elapsed to unzip and trim input fastq.gz files: {trim_and_unzip_fq_end_time - trim_and_unzip_fq_start_time:.2f} seconds")

    minimap2_start_time = time.time()
    run_minimap2(trimmed_fastq_dir)
    minimap2_end_time = time.time()
    print(f"Time elapsed for minimap2 alignment: {minimap2_end_time - minimap2_start_time:.2f} seconds")

    filter_sam_start_time = time.time()
    process_sam_files()
    filter_sam_end_time = time.time()
    print(f"Time elapsed to QC and filter SAM files: {filter_sam_end_time - filter_sam_start_time:.2f} seconds")

    read_counts_dict_start_time = time.time()
    read_counts_dict = get_read_counts_parallel(filtered_samgz_dir)
    read_counts_dict_end_time = time.time()
    print(f"Time elapsed to calculate depths: {read_counts_dict_end_time - read_counts_dict_start_time:.2f} seconds")

    variant_calling_start_time = time.time()
    processed_dfs, removed_counts = process_sam_cs_variants(filtered_samgz_dir)
    variant_calling_end_time = time.time()
    print(f"Time elapsed to call variants: {variant_calling_end_time - variant_calling_start_time:.2f} seconds")

    count_variants_start_time = time.time()
    get_variants_parallel(raw_output_dir, variants_per_sample_dir, snv_per_sample_dir, nm_histograms_dir, removed_counts, read_counts_dict)
    count_variants_end_time = time.time()
    print(f"Time elapsed to calculate variant statistics: {count_variants_end_time - count_variants_start_time:.2f} seconds")


if __name__ == "__main__":
    main()