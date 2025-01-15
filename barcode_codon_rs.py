import os
import yaml
import pandas as pd
import numpy as np
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

### Determine the full path of which minigene_splicing_assay dir was copied to ###
minigene_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(minigene_dir, "minigene_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from config ###
root_output_dir = config["root_output_dir"]
qc_report = config["qc_report"]
mut_region_start = config["mut_region_start"]

### Define variants_info.csv path ###
variant_info_csv =  os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "codon_barcodes_all_info.csv")

### Define output file paths ###
barcode_codon_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "codon_barcodes_info.csv")
non_unique_bc_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "non_unique_codon_barcdoes.txt")
unique_bc_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "unique_codon_barcodes.txt")
codon_depths_pdf_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "variant_codon_depth_plots.pdf")

### Define SNV-related file paths to remove ###
non_unique_SNVs_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "non_unique_SNV_barcodes.txt")
unique_SNVs_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "unique_SNV_barcodes.txt")
SNV_barcode_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "SNV_barcode_info.csv")
SNV_plots_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "SNV_depths_plots.pdf")

### Read variants_info.csv ###
df = pd.read_csv(variant_info_csv)

### Initialise a dictionary to store variant depthsper {barcode}.bam ###
variant_barcode_depth_dict = {}

### Iterate through each row in the variant_info_csv ###
for _, row in df.iterrows():
    variant = row['VARIANT']
    barcodes = row['BARCODE'].split(',')
    depths = [int(depth.strip()) for depth in row['VARIANT_ALLELIC_DEPTH'].split(', ')]

    ### Create dictionary entries with (variant, barcode) as keys and depth as values ###
    for barcode, depth in zip(barcodes, depths):
        variant_barcode_depth_dict[(variant, barcode)] = depth

### Create the dictionary mapping position to reference base ###
pos_to_ref_base = {}

for variant in df['VARIANT']:
    pos = ''.join(filter(str.isdigit, variant))                                                                                                     ### extract position (numeric part) ###
    ref_base = ''.join(filter(str.isalpha, variant.split('>')[0]))                                                                                  ### extract reference base ###
    pos_to_ref_base[int(pos)] = ref_base

### Split and explode 'BARCODE' column ###
df = df.assign(BARCODE=df['BARCODE'].str.split(',')).explode('BARCODE')

### Group by BARCODE and join the VARIANTS ###
df = df.groupby('BARCODE')['VARIANT'].apply(', '.join).reset_index()
df.columns = ['BARCODE', 'VARIANTS']

### Function to determine the codon number ###
def get_codon_number(position):
    return (position - mut_region_start) // 3 + 1

### Function to determine the positions corresponding to a codon number ###
def codon_positions(codon_number):
    start_pos = mut_region_start + (codon_number - 1) * 3
    return [start_pos, start_pos + 1, start_pos + 2]

### Function to create the VARIANT_CODON string ###
def create_variant_codon(variants, codon_number):
    pos_to_obs_base = {int(variant.split('>')[0][:-1]): variant.split('>')[1] for variant in variants.split(', ')}
    positions = codon_positions(codon_number)
    variant_codon = ''
    for pos in positions:
        if pos in pos_to_obs_base:
            variant_codon += pos_to_obs_base[pos]
        else:
            variant_codon += pos_to_ref_base[pos]
    return variant_codon

### Process each row to determine codon numbers and variant codons ###
result_rows = []

for index, row in df.iterrows():
    barcode = row['BARCODE']
    variants = row['VARIANTS'].split(', ')
    codon_dict = {}

    ### Organise variants by codon ###
    for variant in variants:
        pos = int(variant[:-3])                                                                                                                     ### extract the position from the variant string ###
        codon_num = get_codon_number(pos)
        if codon_num not in codon_dict:
            codon_dict[codon_num] = []
        codon_dict[codon_num].append(variant)

    ### Add each group of variants belonging to the same codon as a new row in result_rows ###
    for codon_num, variant_list in codon_dict.items():
        variant_codon = create_variant_codon(', '.join(variant_list), codon_num)
        result_rows.append({
            'BARCODE': barcode,
            'VARIANTS': ', '.join(variant_list),
            'CODON_NUMBER': codon_num,
            'VARIANT_CODON': variant_codon
        })

### Create the unfiltered df from the list of result rows ###
result_df = pd.DataFrame(result_rows)

### Load the QC file ###
qc_xls = pd.ExcelFile(qc_report)
qc_df = pd.read_excel(qc_xls, "VariantProportion")

### Make a copy of qc_df and retain only "AA Position" and "variant_codon" columns ###
qc_df1 = qc_df.copy()
qc_df1 = qc_df1[['AA Position', 'variant_codon']]

### Retain only "wt_codon" and "wt_aa" columns in qc_df###
qc_df = qc_df[["AA Position", "wt_codon", "wt_aa"]]

### Rename AA Position to CODON_NUMBER in qc_df ###
qc_df = qc_df.rename(columns={'AA Position': 'CODON_NUMBER', 'wt_codon': 'WT_CODON', 'wt_aa': 'WT_AA'})

### Create a dictionary that maps AA Position to all its affiliated variant_codon values ###
aa_pos_to_variant_codon = qc_df1.groupby('AA Position')['variant_codon'].apply(set).to_dict()

### Filter result_df based on the dictionary ###
filtered_rows = []

for index, row in result_df.iterrows():
    codon_number = row['CODON_NUMBER']
    variant_codon = row['VARIANT_CODON']

    ### Check if codon_number corresponds to an AA Position in the dictionary ###
    if codon_number in aa_pos_to_variant_codon:
        ### Check if the variant_codon exists in the dictionary for the codon_number ###
        if variant_codon in aa_pos_to_variant_codon[codon_number]:
            filtered_rows.append(row)

### Create the filtered df that omits rows where variant codon is not found in qc file ###
filtered_df = pd.DataFrame(filtered_rows)

grouped_df = filtered_df.groupby(['VARIANTS', 'CODON_NUMBER', 'VARIANT_CODON'])['BARCODE'].apply(', '.join).reset_index()

### Rename columns ###
grouped_df.columns = ['VARIANTS', 'CODON_NUMBER', 'VARIANT_CODON', 'BARCODES']
grouped_df['NUM_BARCODES'] = grouped_df['BARCODES'].apply(lambda x: x.count(',') + 1)

### Split "BARCODES" column into individual barcodes and count frequencies ###
barcode_series = grouped_df['BARCODES'].str.split(', ').explode()
barcode_counts = barcode_series.value_counts()

### Identify non-unique barcodes that exist more than once and save to non_unique_bc.txt file ###
non_unique_barcodes = barcode_counts[barcode_counts > 1]
with open(non_unique_bc_filepath, 'w') as file:
    file.write("Barcodes that exist more than once:\n")
    file.write("BARCODES\n")
    non_unique_barcodes.to_string(file, header=False)

### Identify unique barcodes that exist only once and save to unique_bc.txt ###
unique_barcodes = barcode_counts[barcode_counts == 1]
with open(unique_bc_filepath, 'w') as file:
    file.write("Barcodes that exist only once:\n")
    file.write("BARCODES\n")
    unique_barcodes.to_string(file, header=False)

### Codon to amino acid map ###
codon_to_aa = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

### Function to map codon to amino acid ###
def map_codon_to_aa(codon):
    return codon_to_aa.get(codon.upper(), 'X')                                                                                                      ### return 'X' for unknown codons ###

### Apply the function to "VARIANT_CODON" column to create "VARIANT_AA" column ###
grouped_df['VARIANT_AA'] = grouped_df['VARIANT_CODON'].apply(map_codon_to_aa)

### Create a dictionary from qc_df with CODON_NUMBER as keys and WT_CODON and WT_AA as values to fill in WT info in grouped_df ###
codon_to_wt_info = {row['CODON_NUMBER']: {'WT_CODON': row['WT_CODON'], 'WT_AA': row['WT_AA']} 
                    for _, row in qc_df.iterrows()}

grouped_df['WT_CODON'] = np.nan                                                                                                                     ### add empty columns "WT_AA" and "WT_CODON" to grouped_df ###
grouped_df['WT_AA'] = np.nan

### Fill in the "WT_AA" and "WT_CODON" columns using the dictionary and grouped_df["CODON_NUMBER"] ###
for index, row in grouped_df.iterrows():
    codon_number = row['CODON_NUMBER']
    if codon_number in codon_to_wt_info:
        grouped_df.at[index, 'WT_CODON'] = codon_to_wt_info[codon_number]['WT_CODON']
        grouped_df.at[index, 'WT_AA'] = codon_to_wt_info[codon_number]['WT_AA']

grouped_df = grouped_df[["CODON_NUMBER", "VARIANTS", "BARCODES", "NUM_BARCODES", "WT_CODON", "VARIANT_CODON", "WT_AA", "VARIANT_AA"]]

def calculate_variant_depth(variants, barcodes, depth_dict):
    """
    Calculates the total depth for each variant based on its associated barcodes and variant_barcode_depth_dict dictionary

    Parameters:
    * variants (str): string containing variants separated by ', '
    * barcodes (str): string containing barcodes separated by ', ' corresponding to the variants it's associated with
    * depth_dict (dict): variant_barcode_depth_dict where keys are tuples of (variant, barcode) and values are the variant allelic depths seen per {barcode}.bam

    Workflow:
    1. Split the 'variants' and 'barcodes' strings into a lists ('variants_list' and 'barcodes_list' respectively) using ', ' as delimiter 
    2. Initialise an empty list 'variant_depths' to store depth values for each variant
    3. Iterate over each variant in 'variants_list':
        * initialise a variable 'total_depth' to 0 for the current variant
        * iterate over each barcode in 'barcodes_list'
            - use variant_barcode_depth_dict to get the depth value for the specific variant-barcode pair
            - adds depth value to 'total_depth' if the pair exists in variant_barcode_depth_dict; 0 if pair does not exist
        * append formatted string "variant: depth" to 'variant_depths' list.
    4. Join all the elements in 'variant_depths' using ', ' and return resulting string
    """
    ### Split the variants and barcodes by ', ' ###
    variants_list = variants.split(', ')
    barcodes_list = barcodes.split(', ')
    
    variant_depths = []
    
    ### Iterate over each variant ###
    for variant in variants_list:
        total_depth = 0
        ### Sum the depths for the current variant using the barcodes ###
        for barcode in barcodes_list:
            total_depth += depth_dict.get((variant, barcode), 0)
        ### Append the result as "variant: depth" ###
        variant_depths.append(f"{variant}: {total_depth}")
        
    ### Return the joined string of depths ###
    return ', '.join(variant_depths)

### Apply calculate_variant_depth function to create new "VARIANTS_DEPTH" column with 'variant_depths' string values ###
grouped_df['VARIANTS_DEPTH'] = grouped_df.apply(lambda row: calculate_variant_depth(row['VARIANTS'], row['BARCODES'], variant_barcode_depth_dict), axis=1)

def find_smallest_depth(variants_depth):
    ### Split the variants and depths by ', ' and get the depth values ###
    depths = [int(variant_depth.split(': ')[1]) for variant_depth in variants_depth.split(', ')]
    ### Return the smallest value ###
    return min(depths)

### Apply find_smallest_depth function to create new "VARIANT_CODON_DEPTH" column with minimum depth found for that variant combination ###
grouped_df['VARIANT_CODON_DEPTH'] = grouped_df['VARIANTS_DEPTH'].apply(find_smallest_depth)

### Rearrange columns and save as barcode_codon_info.csv ###
grouped_df_filt = grouped_df[["CODON_NUMBER", "VARIANTS", "BARCODES", "NUM_BARCODES", "WT_CODON", "VARIANT_CODON", "WT_AA", "VARIANT_AA", "VARIANT_CODON_DEPTH"]]
grouped_df_filt.to_csv(barcode_codon_filepath, index=False)

### Add new column 'CODON_#_VARIANT_CODON' with values that represent codon_number | variant codon ###
grouped_df_filt['CODON_#_VARIANT_CODON'] = grouped_df_filt['CODON_NUMBER'].astype(str) + " | " + grouped_df_filt['VARIANT_CODON']

### Plot variant codon depth information ###
with PdfPages(codon_depths_pdf_filepath) as pdf:
    ### histogram of depths vs codon_number | variant codon ###

    ### colour palette for histogram ###
    color_cycle = ['pink', 'mediumvioletred', 'darkmagenta']                                                                                        ### main 3 colors to repeat for the codons ###                                                                                           
    bar_colors = [color_cycle[(codon_num - 1) % len(color_cycle)] for codon_num in grouped_df_filt['CODON_NUMBER']]                                 ### assign colors based on the codon number modulo the length of the color cycle ###

    ### plot histogram with the bars coloured by the codon number ###
    fig, ax = plt.subplots(figsize=(50, 20))
    bars = ax.bar(grouped_df_filt['CODON_#_VARIANT_CODON'], grouped_df_filt['VARIANT_CODON_DEPTH'], color=bar_colors)

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('Codon Number | Variant Codon', fontsize=40, labelpad=30)
    ax.set_xticks(np.arange(len(grouped_df_filt)))                                                                                                
    ax.set_xticklabels(grouped_df_filt['CODON_#_VARIANT_CODON'], rotation=90, fontsize=2)
    ax.set_xlim(-1, len(grouped_df_filt))

    ### configure y-axis label and tick labels ###
    ax.set_ylabel('Variant Codon Depth', fontsize=40, labelpad=30)                                                                         
    yticks = np.arange(0, grouped_df_filt['VARIANT_CODON_DEPTH'].max() + 1000, 500)                                                                 
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=30)

    ### set title ###
    ax.set_title('Variant Codon Depth Distribution', fontsize=50, pad=50)

    ### plot gridlines behind bars ###
    ax.grid(axis='y', linestyle='--', alpha=0.3)                                                                                                   
    ax.set_axisbelow(True)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.15)                                                                                

    ### generate the legend labels dynamically based on codon colour pattern ###
    legend_labels = []
    for i, color in enumerate(color_cycle):
        ### calculate which codons correspond to this color (1-based indexing) ###
        codons = [str(codon) for codon in range(i + 1, max(grouped_df_filt['CODON_NUMBER']) + 1, len(color_cycle))]
        legend_labels.append(f"Codons {', '.join(codons)}")

    ### create legend patches based on the color cycle and corresponding codons ###
    legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(color_cycle, legend_labels)]

    ### place the legend at the bottom center of the plot ###
    legend = ax.legend(
        handles=legend_patches, 
        loc='upper center', 
        bbox_to_anchor=(0.5, -0.1),  
        fontsize=19,
        ncol=len(color_cycle),
        frameon=False
    )

    ### save histogram to PDF and close ###
    pdf.savefig(fig)
    plt.close(fig)

    ### Plot heatmap visualising read depths per variant codon ###
    depth_heatmap_data = grouped_df_filt.pivot(index='VARIANT_CODON', columns='CODON_NUMBER', values='VARIANT_CODON_DEPTH')
    fig, ax = plt.subplots(figsize=(50, 35))
    depth_heatmap = sns.heatmap(depth_heatmap_data, cmap='magma_r', annot=False, cbar_kws={'label': 'Variant Codon Depth'}, ax=ax)

    ### access colorbar from heatmap and configure tick and label sizes ###
    depth_cbar = depth_heatmap.collections[0].colorbar
    depth_cbar.ax.tick_params(labelsize=30) 
    depth_cbar.set_label('Variant Codon Depth', size=35, labelpad=40) 

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('Codon Number', fontsize=40, labelpad=30)
    ax.set_xticks(np.arange(0.5, len(depth_heatmap_data.columns), 1))
    ax.set_xticklabels(depth_heatmap_data.columns, rotation=90, fontsize=20)

    ### configure y-axis label and tick labels ###
    ax.set_ylabel('Variant Codon', fontsize=40, labelpad=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=30)

    ### set title ###
    ax.set_title('Heatmap of Variant Codon Depth Across Codon Numbers', fontsize=50, pad=50)
    
    plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.1)

    ### save heatmap to same PDF and close ###
    pdf.savefig(fig)
    plt.close(fig)

    ### Plot heatmap visualising number of barcodes per variant codon ###
    barcode_heatmap_data = grouped_df_filt.pivot(index='VARIANT_CODON', columns='CODON_NUMBER', values='NUM_BARCODES')
    fig, ax = plt.subplots(figsize=(50, 35))
    barcode_heatmap = sns.heatmap(barcode_heatmap_data, cmap='magma_r', annot=False, cbar_kws={'label': 'Number of Barcodes'}, ax=ax)

    ### access colorbar from heatmap and configure tick and label sizes ###
    barcode_cbar = barcode_heatmap.collections[0].colorbar
    barcode_cbar.ax.tick_params(labelsize=30)
    barcode_cbar.set_label('Number of Barcodes', size=35, labelpad=40)

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('Codon Number', fontsize=40, labelpad=30)
    ax.set_xticks(np.arange(0.5, len(barcode_heatmap_data.columns), 1))
    ax.set_xticklabels(barcode_heatmap_data.columns, rotation=90, fontsize=20)

    ### configure y-axis label and tick labels ###
    ax.set_ylabel('Variant Codon', fontsize=40, labelpad=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=30)

    ax.set_title('Heatmap of Number of Barcodes Across Codon Numbers', fontsize=50, pad=50)
    
    plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.1)

    ### save heatmap to same PDF and close ###
    pdf.savefig(fig)
    plt.close(fig)

### Remove SNV-related files if they exist ###
file_paths = [non_unique_SNVs_filepath, unique_SNVs_filepath, SNV_barcode_filepath, SNV_plots_filepath]

for file_path in file_paths:
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print(f"{file_path} does not exist - nothing to remove")

