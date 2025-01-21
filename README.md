### Minigene Pipeline
Single nucleotide variants (SNVs) were introduced into the mutagenesis region, spanning both introns and exons of a minigene plasmid containing the target exon. Additionally, random 12 bp barcodes, flanked by a constant 8 bp prefix and suffix, were also incorporated into the plasmid's exon. The modified minigene plasmid was then transfected into iPSCs and HEK cells and sequenced post-transcription using ONT. The pipeline has 3 features. The first feature, at the most foundational level, identifies barcode-SNV relationship, while the second feature identifies barcode-amino acid relationship. The third feature permits splicing analysis, which includes processing the ONT cDNA data, produces statistical reports in CSVs and visualisations in PDFs of the different isoforms found as a result of different SNVs. The workflow is illustrated below.

<img width="1386" alt="Screenshot 2024-06-26 at 1 54 09 PM" src="https://github.com/AshleyLab/VxE_Map/assets/96602087/956a3e16-a8a3-4c5c-8070-8db5cdf986ec">

#### Installation of Minigene Pipeline on Sherlock
1. Log in to Sherlock and navigate to directory where results are to be stored i.e. target_dir
2. Copy minigene_pipeline directory into your target directory using the below command  
   ```  
   cp -r /oak/stanford/groups/euan/projects/variant_effect_mapping/minigene_pipeline /path/to/target_dir
#### Setting Up Minigene Pipeline's Environment
3. Navigate into copied minigene_pipeline directory
   ```
   cd /path/to/target_dir/minigene_pipeline
4. Ensure that Python 3.12 is installed on your system. Install it otherwise before proceeding.
5. Install singularity image of DeepVariant (v1.6.1 (CPU version))
   ```
   BIN_VERSION="1.6.1"
   singularity pull docker://google/deepvariant:"${BIN_VERSION}"
6. Create and activate minigene_env 
   ```
   python3.12 -m venv minigene_env      # create new python virtual env
   source minigene_env/bin/activate     # activate minigene_env
7. Install packages
   ```
   python3 setup_minigene_env.py 
#### Execution
8. Edit minigene_config.yaml file to suit your data
   ```
   vim minigene_config.yaml
   
   # press 'i' to edit and make changes; when done, press 'esc', ':wq' and hit 'enter'
9. Navigate into log directory
   ```
   cd log

   # full path of log dir = "/path/to/target_dir/minigene_pipeline/log
10. Execute minigene pipeline
    ```
    # For analysis of only barcode-SNV relationship, do:
    sbatch /path/to/target_dir/minigene_pipeline/minigene_pipeline.sbatch /path/to/target_dir/minigene_pipeline "" false
    
    # For analysis of barcode-codon relationship, do:
    sbatch /path/to/target_dir/minigene_pipeline/minigene_pipeline.sbatch /path/to/target_dir/minigene_pipeline barcode_codon_rs.py
    
    # For full analysis including isoform processing, do:
    sbatch /path/to/target_dir/minigene_pipeline/minigene_pipeline.sbatch /path/to/target_dir/minigene_pipeline
#### Tools
samtools: [samtools](https://github.com/samtools/samtools)  
chopper: [chopper](https://github.com/wdecoster/chopper)  
minimap2: [minimap2](https://github.com/lh3/minimap)  
deepvariant: [deepvariant](https://github.com/google/deepvariant)   
GMAP: [GMAP](http://research-pub.gene.com/gmap)