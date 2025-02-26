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
mut_regions = config["mut_regions"]

### Define variants_info.csv path ###
variant_info_csv =  os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "codon_barcodes_all_info.csv")

### Define output file paths ###
barcode_codon_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "codon_barcodes_results.csv")
bc_more_than_1_codon_file = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "barcodes_in_more_than_1_codon.csv")
codon_depths_pdf_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "variant_codon_plots.pdf")
unique_barcodes_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "unique_barcodes.txt")

### Define SNV-related file paths to remove ###
non_unique_SNVs_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "non_unique_SNV_barcodes.txt")
unique_SNVs_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "unique_SNV_barcodes.txt")
SNV_barcode_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "SNV_barcode_info.csv")
SNV_plots_filepath = os.path.join(root_output_dir, "variant_barcode_results", "postprocess_variants", "SNV_depths_plots.pdf")


#####################################################################################################################################################################
                                                ### CREATE DICTIONARY THAT STORES THE DEPTH OF EACH SNV_BARCODE ###
#####################################################################################################################################################################

df = pd.read_csv(variant_info_csv)

### Initialise a dictionary to store variant depths per {barcode}.bam {('SNV', 'barcode'): depth, ...} ###
variant_barcode_depth_dict = {}

### Iterate through each row in the variant_info_csv ###
for _, row in df.iterrows():
    variant = row['VARIANT']
    barcodes = row['BARCODE'].split(',')
    depths = [int(depth.strip()) for depth in row['VARIANT_ALLELIC_DEPTH'].split(', ')]

    ### Create dictionary entries with (variant, barcode) as keys and depth as values ###
    for barcode, depth in zip(barcodes, depths):
        variant_barcode_depth_dict[(variant, barcode)] = depth


###################################################################################################################################################################
                                        ### IDENTIFY UNIQUE SNVS OF EACH BARCODE (SHOULD BE SNVS WITHIN THE SAME CODON) ###
###################################################################################################################################################################

### Split and explode 'BARCODE' column ###
df = df.assign(BARCODE=df['BARCODE'].str.split(',')).explode('BARCODE')

### Group by BARCODE and join the VARIANTS ###
df = df.groupby('BARCODE')['VARIANT'].apply(', '.join).reset_index()
df.columns = ['BARCODE', 'VARIANTS']

###################################################################################################################################################################
                                                ### GET CDNA POSITIONS OF EACH CODON IN EACH MUTAGENESIS REGION ###
###################################################################################################################################################################

### Read QC file into df ###
qc_xls = pd.ExcelFile(qc_report)
qc_df = pd.read_excel(qc_xls, "VariantProportion")

qc_df = qc_df[['AA Position', 'wt_codon']]
qc_df.rename(columns={'AA Position': 'AA_Position'}, inplace=True)
qc_df.drop_duplicates(subset='AA_Position', keep='first', inplace=True)
qc_df = qc_df.reset_index(drop=True)

### Assign number to each mutagenesis region ###
qc_df['Region'] = (qc_df['AA_Position'].diff() > 1).cumsum() + 1

### Precompute ranges for each region ###
mut_regions_ranges = {i + 1: range(start, end + 1) for i, (start, end) in enumerate(mut_regions)}

### Initialize a tracker for start index in each region ###
region_trackers = {region: 0 for region in mut_regions_ranges}

### Function to assign unique Variant_Positions ###
def assign_variant_positions(row):
    region = row["Region"]
    start_index = region_trackers[region]
    region_range = list(mut_regions_ranges[region])
    positions = region_range[start_index:start_index + 3]
    region_trackers[region] += 3                                                                                                ### move the tracker forward by 3 (each codon has 3 bases) ###
    return ", ".join(map(str, positions))

### Add Variant_Positions column ###
qc_df["Variant_Positions"] = qc_df.apply(assign_variant_positions, axis=1)

### Verify Variant_Positions corresponds to range in mut_regions ###
for region, range_values in mut_regions_ranges.items():
    last_row = qc_df[qc_df["Region"] == region].iloc[-1]
    last_variant_positions = list(map(int, last_row["Variant_Positions"].split(", ")))
    assert range_values[-1] == last_variant_positions[-1], (
        f"Region {region}: Expected last value {range_values[-1]}, but got {last_variant_positions[-1]}"
    )

print("\nVerification completed: All last Variant_Positions values match their respective range end values.")

###################################################################################################################################################################
                                                        ### OBTAIN AA_POSITION AND WT_CODON OF THE SNVS ###
###################################################################################################################################################################

position_mapping = {}

for _, row in qc_df.iterrows():
    positions = map(int, row['Variant_Positions'].split(', '))
    for pos in positions:
        position_mapping[pos] = (row['AA_Position'], row['wt_codon'])

def map_variants_to_positions(variants):
    positions = [int(var[:-3]) for var in variants.split(', ')]                                                                 ### extract {pos} from each variant ({pos}{ref}>{obs}) ###
    aa_positions = set()
    wt_codons = set()

    for pos in positions:
        if pos in position_mapping:
            aa_position, wt_codon = position_mapping[pos]
            aa_positions.add(aa_position)
            wt_codons.add(wt_codon)

    return ', '.join(map(str, sorted(aa_positions))), ', '.join(sorted(wt_codons))

df[['AA_Position', 'wt_codon']] = df['VARIANTS'].apply(
    lambda variants: pd.Series(map_variants_to_positions(variants))
)

#####################################################################################################################################################################
                                ### CREATE DICTIONARY THAT MAPS EACH CDNA POSITION IN MUTAGENESIS REGION(S) TO ITS REFERENCE BASE ###
#####################################################################################################################################################################

### Create the dictionary mapping position to reference base ###
pos_to_ref_base = {}

for _, row in qc_df.iterrows():
    positions = list(map(int, row['Variant_Positions'].split(', ')))                                                            ### split Variant_Positions into a list of integers ###
    for pos, letter in zip(positions, row['wt_codon']):                                                                         ### iterate through each position and its corresponding letter in wt_codon ###
        pos_to_ref_base[pos] = letter

#####################################################################################################################################################################
                            ### IDENTIFY ANY MISSING SNVS PER CODON AND ADD THAT INTO VARIANTS_1 AS {MISSING_SNV_POS}{REF}>{REF} ### 
#####################################################################################################################################################################

### Function to parse through df and add missing cdna_pos values ###
def add_missing_cdna_pos_to_variants(df, qc_df, pos_to_ref_base):
    ### Create a dictionary mapping AA_Position to Variant_Positions ###
    aa_to_variant_positions = {}
    for _, row in qc_df.iterrows():
        aa_position = row["AA_Position"]
        variant_positions = list(map(int, row["Variant_Positions"].split(", ")))
        aa_to_variant_positions[aa_position] = variant_positions

    ### Iterate through each row in df ###
    for index, row in df.iterrows():
        variants = row["VARIANTS"].split(", ")
        aa_positions = list(map(int, row["AA_Position"].split(", ")))

        ### Parse existing cdna_pos from variants ###
        existing_cdna_pos = set(
            int(variant.split(">")[0][:-1]) for variant in variants
        )

        ### Iterate through each AA_Position ###
        for aa_pos in aa_positions:
            if aa_pos in aa_to_variant_positions:
                ### Get the associated cdna_pos values for current AA_Position ### 
                associated_cdna_pos = set(aa_to_variant_positions[aa_pos])

                ### Find missing cdna_pos values ###
                missing_cdna_pos = associated_cdna_pos - existing_cdna_pos

                ### Construct and add missing variants ###
                for missing_pos in missing_cdna_pos:
                    ref_base = pos_to_ref_base.get(missing_pos, "N")
                    missing_variant = f"{missing_pos}{ref_base}>{ref_base}"
                    variants.append(missing_variant)

        ### Update the row in df ###
        df.at[index, "VARIANTS_1"] = ", ".join(variants)

    return df

updated_df = add_missing_cdna_pos_to_variants(df, qc_df, pos_to_ref_base)

#####################################################################################################################################################################
                        ### SORT THE SNV1, SNV2, SNV3(, SNV4, SNV5, SNV6, ...) IN VARIANTS_1 IN ASCENDING ORDER OF THEIR CDNA POSITIONS ###
#####################################################################################################################################################################

### Sort variants by cdna_pos ###
def sort_variants(variant_string):
    variants = variant_string.split(", ")
    sorted_variants = sorted(variants, key=lambda x: int(x.split(">")[0][:-1]))
    return ", ".join(sorted_variants)

updated_df["VARIANTS_1"] = updated_df["VARIANTS_1"].apply(sort_variants)

#####################################################################################################################################################################
                                                    ### DETERMINE VARIANT CODON BY JOINING {OBS} OF EACH SNV ###
#####################################################################################################################################################################

### Derive VARIANT_CODON ###
def generate_variant_codon(variants):
    variant_list = variants.split(", ")
    variant_codons = []

    ### Iterate through chunks of 3 consecutive variant values ###
    for i in range(0, len(variant_list), 3):
        triplet = variant_list[i:i+3]
        #### Extract the {obs} values and join them into a codon ###
        codon = "".join([variant.split(">")[1] for variant in triplet])
        variant_codons.append(codon)

    return ", ".join(variant_codons)

### Call generate_variant_codon to create VARIANT_CODON column ###
updated_df["VARIANT_CODON"] = updated_df["VARIANTS_1"].apply(generate_variant_codon)

updated_df.rename(columns={'wt_codon': 'WT_CODON', 'AA_Position': 'AA_POSITION'}, inplace=True)

#####################################################################################################################################################################
                                ### REGENERATE WILDTYPE CODON TO MATCH THE ORDER OF SORTED VARIANTS (VARIANTS_1) AND VARIANT_CODON ###
#####################################################################################################################################################################

### Derive WILDTYPE_CODON ###
def generate_wt_codon(variants):
    variant_list = variants.split(', ')
    wt_codons = []
    for i in range(0, len(variant_list), 3):                                                                                    ### Process 3 variants with consecutive positions belonging to same codom at a time ###
        codon = "".join([v.split('>')[0][-1] for v in variant_list[i:i+3]])                                                     ### Extract and join {ref} ###
        wt_codons.append(codon)
    return ", ".join(wt_codons)

### Call generate_wt_codon function create WILDTYPE_CODON column ###
updated_df['WILDTYPE_CODON'] = updated_df['VARIANTS_1'].apply(generate_wt_codon)

updated_df = updated_df[["BARCODE", "VARIANTS", "VARIANTS_1", "WILDTYPE_CODON", "VARIANT_CODON"]]


#####################################################################################################################################################################
                            ### REGENERATE AA_POSITION (CODON_NUMBER) TO MATCH THE ORDER OF SORTED VARIANTS (VARIANTS_1) AND VARIANT_CODON ###
#####################################################################################################################################################################

AApos_codonpos_dict = {}

for _, row in qc_df.iterrows():
    aa_pos = row['AA_Position']
    variant_positions = list(map(int, row['Variant_Positions'].split(', ')))                                                    ### split Variant_Positions into a list of integers ###                            
    ### Append or create the entry for the wt_codon key ###
    if aa_pos not in AApos_codonpos_dict:
        AApos_codonpos_dict[aa_pos] = variant_positions
    else:
        AApos_codonpos_dict[aa_pos].extend(variant_positions)

### Reverse the dictionary to create a mapping of cdna_pos to AA_Position ###
cdna_to_AApos = {pos: aa_pos for aa_pos, positions in AApos_codonpos_dict.items() for pos in positions}

def find_codon_numbers(variants, cdna_to_AApos):
    """Find codon numbers for sets of 3 cdna_pos values in variants."""
    ### Extract cdna_pos values from the variants ###
    cdna_positions = [int(variant.split('>')[0][:-1]) for variant in variants.split(', ')]

    ### Group cdna_positions into sets of 3 ###
    codon_numbers = []
    for i in range(0, len(cdna_positions), 3):
        ### Get a set of 3 positions ###
        codon_set = cdna_positions[i:i+3]

        ### Check if all 3 positions map to the same AA_Position ###
        aa_positions = {cdna_to_AApos.get(pos) for pos in codon_set if pos in cdna_to_AApos}
        if len(aa_positions) == 1:  
            codon_numbers.append(next(iter(aa_positions)))
        else:
            codon_numbers.append(', '.join(map(str, sorted(aa_positions))))

    ### Join codon numbers as a string separated by ', ' ###
    return ', '.join(map(str, codon_numbers))

updated_df['CODON_NUMBER'] = updated_df['VARIANTS_1'].apply(lambda x: find_codon_numbers(x, cdna_to_AApos))

#####################################################################################################################################################################
                    ### CREATE DICTIONARY OF CODON TO AMINO ACIDS + DERIVE WILDTYPE_AA FROM WILDTYPE_CODON AND VARIANT_AA FROM VARIANT_CODON ###
#####################################################################################################################################################################

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

### Function to map codons to amino acids ###
def codon_to_aa_mapping(codon_str):
    codons = codon_str.split(', ')
    amino_acids = [codon_to_aa.get(codon, '?') for codon in codons]
    return ', '.join(amino_acids)

### Call codon_to_aa_mapping to create WILDTYPE_AA and VARIANT_AA columns ###
updated_df["WILDTYPE_AA"] = updated_df["WILDTYPE_CODON"].apply(codon_to_aa_mapping)
updated_df["VARIANT_AA"] = updated_df["VARIANT_CODON"].apply(codon_to_aa_mapping)
filtered_df = updated_df[["CODON_NUMBER", "BARCODE", "VARIANTS", "WILDTYPE_CODON", "VARIANT_CODON", "WILDTYPE_AA", "VARIANT_AA"]]

### Sort the dataframe by CODON_NUMBER ###
filtered_df = filtered_df.sort_values(by='CODON_NUMBER').reset_index(drop=True)

#####################################################################################################################################################################
                                                    ### REMOVE BARCODES THAT IS TAGGED TO DIFFERENT CODONS ###
#####################################################################################################################################################################

### ID barcodes tagged to more than 1 codon and save ###
non_unique_bc = filtered_df[filtered_df['CODON_NUMBER'].str.contains(',', na=False)]
non_unique_bc = non_unique_bc.reset_index(drop=True)
non_unique_bc.to_csv(bc_more_than_1_codon_file, index=False)

### filter out barcodes tagged to more than 1 codon before further analyses ###
filtered_df = filtered_df[~filtered_df['CODON_NUMBER'].str.contains(',')].reset_index(drop=True)

#####################################################################################################################################################################
                                                ### FIND NUMBER OF UNIQUE BARCODES PER VARIANT_CODON ###
#####################################################################################################################################################################

grouped_df = filtered_df.groupby(['VARIANTS', 'CODON_NUMBER', 'WILDTYPE_CODON', 'VARIANT_CODON', 'WILDTYPE_AA', 'VARIANT_AA'])['BARCODE'].apply(', '.join).reset_index()
grouped_df['NUM_BARCODES'] = grouped_df['BARCODE'].apply(lambda x: x.count(',') + 1)
grouped_df['CODON_NUMBER'] = grouped_df['CODON_NUMBER'].astype(int)

#####################################################################################################################################################################
                                                                ### GET DEPTH OF VARIANT CODON ###
#####################################################################################################################################################################

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
grouped_df['VARIANTS_DEPTH'] = grouped_df.apply(lambda row: calculate_variant_depth(row['VARIANTS'], row['BARCODE'], variant_barcode_depth_dict), axis=1)

def find_smallest_depth(variants_depth):
    ### Split the variants and depths by ', ' and get the depth values ###
    depths = [int(variant_depth.split(': ')[1]) for variant_depth in variants_depth.split(', ')]
    ### Return the smallest value ###
    return min(depths)

### Apply find_smallest_depth function to create new "VARIANT_CODON_DEPTH" column with minimum depth found for that variant combination ###
grouped_df['VARIANT_CODON_DEPTH'] = grouped_df['VARIANTS_DEPTH'].apply(find_smallest_depth)

### Re-read QC file into df to make dictionary to map AA_position global to codon number ###
qc_xls_report = pd.ExcelFile(qc_report)
qc_df_report = pd.read_excel(qc_xls, "VariantProportion")

qc_df_report = qc_df_report[["AA Position", "aa_position_global"]]
qc_df_report = qc_df_report.drop_duplicates(subset=["AA Position"], keep='first')                                               ### drop duplicates based on 'AA Position' and keep the first occurrence ###
qc_df_report = qc_df_report.reset_index(drop=True)

### create dictionary with AA Position as key and aa_position_global as value ###
aa_global_dict = qc_df_report.set_index('AA Position')['aa_position_global'].to_dict()

### create AA_CONSEQUENCE column ###
grouped_df_filt['AA_POSITION'] = grouped_df_filt['CODON_NUMBER'].map(aa_global_dict)
grouped_df_filt['AA_CONSEQUENCE'] = 'p.' + grouped_df_filt['WILDTYPE_AA'] + grouped_df_filt['AA_POSITION'].astype(str) + grouped_df_filt['VARIANT_AA']

### Rearrange columns and save as barcode_codon_info.csv ###
grouped_df_filt = grouped_df[["CODON_NUMBER", "VARIANTS", "AA_CONSEQUENCE", "BARCODE", "NUM_BARCODES", "WILDTYPE_CODON", "VARIANT_CODON", "WILDTYPE_AA", "VARIANT_AA", "VARIANT_CODON_DEPTH"]]

grouped_df_filt.to_csv(barcode_codon_filepath, index=False)

### Save unique barcodes to a text file ###
all_barcodes = [barcode.strip() for barcodes in grouped_df_filt["BARCODE"] for barcode in barcodes.split(", ")]                              ### flatten the barcodes into a list ###
unique_barcodes = set()
duplicates = set()

for barcode in all_barcodes:
    if barcode in unique_barcodes:
        duplicates.add(barcode)
    else:
        unique_barcodes.add(barcode)

### write unique barcodes to a text file ###
with open(unique_barcodes_filepath, "w") as f:
    for barcode in sorted(unique_barcodes):  
        f.write(barcode + "\n")

#####################################################################################################################################################################
                                                                ### PLOT HISTOGRAM AND HEATMAPS ###
#####################################################################################################################################################################

### Add new column 'CODON_#_VARIANT_CODON' with values that represent aa_position | variant codon ###
grouped_df_filt['CODON_#_VARIANT_CODON'] = grouped_df_filt['AA_POSITION'].astype(str) + " | " + grouped_df_filt['VARIANT_CODON']

### Plot variant codon depth information ###
with PdfPages(codon_depths_pdf_filepath) as pdf:
    ### histogram of depths vs aa_position | variant codon ###

    ### colour palette for histogram ###
    color_cycle = ['pink', 'mediumvioletred', 'darkmagenta']                                                                                        ### main 3 colors to repeat for the codons ###                                                                                           
    bar_colors = [color_cycle[(codon_num - 1) % len(color_cycle)] for codon_num in grouped_df_filt['AA_POSITION']]                                 ### assign colors based on the codon number modulo the length of the color cycle ###

    ### plot histogram with the bars coloured by the codon number ###
    fig, ax = plt.subplots(figsize=(50, 20))
    bars = ax.bar(grouped_df_filt['CODON_#_VARIANT_CODON'], grouped_df_filt['VARIANT_CODON_DEPTH'], color=bar_colors)

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('AA Position | Variant Codon', fontsize=40, labelpad=30)
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
        codons = [str(codon) for codon in range(i + 1, max(grouped_df_filt['AA_POSITION']) + 1, len(color_cycle))]
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
    depth_heatmap_data = grouped_df_filt.pivot(index='VARIANT_CODON', columns='AA_POSITION', values='VARIANT_CODON_DEPTH')
    fig, ax = plt.subplots(figsize=(50, 35))
    depth_heatmap = sns.heatmap(depth_heatmap_data, cmap='magma_r', annot=False, cbar_kws={'label': 'Variant Codon Depth'}, ax=ax)

    ### access colorbar from heatmap and configure tick and label sizes ###
    depth_cbar = depth_heatmap.collections[0].colorbar
    depth_cbar.ax.tick_params(labelsize=30)
    depth_cbar.set_label('Variant Codon Depth', size=35, labelpad=40)

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('AA Position', fontsize=40, labelpad=30)
    ax.set_xticks(np.arange(0.5, len(depth_heatmap_data.columns), 1))
    ax.set_xticklabels(depth_heatmap_data.columns, rotation=90, fontsize=20)

    ### configure y-axis label and tick labels ###
    ax.set_ylabel('Variant Codon', fontsize=40, labelpad=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=30)

    ### set title ###
    ax.set_title('Heatmap of Variant Codon Depth Across Amino Acid Positions', fontsize=50, pad=50)

    plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.1)

    ### save heatmap to same PDF and close ###
    pdf.savefig(fig)
    plt.close(fig)

    ### Plot heatmap visualising number of barcodes per variant codon ###
    barcode_heatmap_data = grouped_df_filt.pivot(index='VARIANT_CODON', columns='AA_POSITION', values='NUM_BARCODES')
    fig, ax = plt.subplots(figsize=(50, 35))
    barcode_heatmap = sns.heatmap(barcode_heatmap_data, cmap='magma_r', annot=False, cbar_kws={'label': 'Number of Barcodes'}, ax=ax)

    ### access colorbar from heatmap and configure tick and label sizes ###
    barcode_cbar = barcode_heatmap.collections[0].colorbar
    barcode_cbar.ax.tick_params(labelsize=30)
    barcode_cbar.set_label('Number of Barcodes', size=35, labelpad=40)

    ### configure x-axis label and tick labels ###
    ax.set_xlabel('AA Position', fontsize=40, labelpad=30)
    ax.set_xticks(np.arange(0.5, len(barcode_heatmap_data.columns), 1))
    ax.set_xticklabels(barcode_heatmap_data.columns, rotation=90, fontsize=20)

    ### configure y-axis label and tick labels ###
    ax.set_ylabel('Variant Codon', fontsize=40, labelpad=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=30)

    ax.set_title('Heatmap of Number of Barcodes Across Amino Acid Positions', fontsize=50, pad=50)

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