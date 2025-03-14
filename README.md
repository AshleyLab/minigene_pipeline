### Minigene Pipeline
Single nucleotide variants (SNVs) were introduced into the mutagenesis region, spanning both introns and exons of a minigene plasmid containing the target exon. Additionally, random 12 bp barcodes, flanked by a constant 8 bp prefix and suffix, were also incorporated into the plasmid's exon. The modified minigene plasmid was then transfected into iPSCs and HEK cells and sequenced post-transcription using ONT. The pipeline has 3 features. The first feature, at the most foundational level, identifies barcode-SNV relationship, while the second feature identifies barcode-amino acid relationship. The third feature permits splicing analysis, which includes processing the ONT cDNA data, produces statistical reports in CSVs and visualisations in PDFs of the different isoforms found as a result of different SNVs. The workflow is illustrated below.

<img width="1711" alt="Image" src="https://github.com/user-attachments/assets/fef6b66e-def9-451a-9d4c-22338d8a1f76" />

#### Pre-Requisites
* Python 3.12
    * Ensure that Python 3.12 is installed on your system. Install it otherwise before proceeding.
* Singularity image of DeepVariant (v1.6.1 (CPU version))
    ```
    BIN_VERSION="1.6.1"
    singularity pull docker://google/deepvariant:"${BIN_VERSION}"
#### Installation of Minigene Pipeline 
1. Navigate into directory where minigene pipeline is to be  installed in
   ```
   wget https://github.com/AshleyLab/minigene_pipeline/releases/download/v1.0/minigene_pipeline.tar.gz
   tar -xvzf minigene_pipeline.tar.gz
#### Set Up Minigene Pipeline's Environment
2. Navigate into minigene_pipeline directory
   ```
   cd minigene_pipeline
3. Create and activate minigene_env 
   ```
   python3.12 -m venv minigene_env      # create new python virtual env
   source minigene_env/bin/activate     # activate minigene_env
4. Install packages
   ```
   python3 setup_minigene_env.py 
#### Configure minigene_config.yaml file
5. Edit minigene_config.yaml file to suit your data
   ```
   vim minigene_config.yaml
   
   # press 'i' to edit and make changes; when done, press 'esc', ':wq' and hit 'enter'
#### Execution
6. For splice analysis:
   ```
   python3 run_minigene_splice.py
#### Tools
samtools: [samtools](https://github.com/samtools/samtools)  
chopper: [chopper](https://github.com/wdecoster/chopper)  
minimap2: [minimap2](https://github.com/lh3/minimap)  
deepvariant: [deepvariant](https://github.com/google/deepvariant)   
GMAP: [GMAP](http://research-pub.gene.com/gmap)