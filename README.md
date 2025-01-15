### Minigene Pipeline
Single nucleotide variants (SNVs) were introduced into the mutagenesis region, spanning both introns and exons of a minigene plasmid containing the target exon. Additionally, random 12 bp barcodes, flanked by a constant 8 bp prefix and suffix, were also incorporated into the plasmid's exon. The modified minigene plasmid was then transfected into iPSCs and HEK cells and sequenced post-transcription using ONT. The pipeline has 3 features. The first feature, at the most foundational level, identifies barcode-SNV relationship, while the second feature identifies barcode-amino acid relationship. The third feature permits splicing analysis, which includes processing the ONT cDNA data, produces statistical reports in CSVs and visualisations in PDFs of the different isoforms found as a result of different SNVs. The workflow is illustrated below.

<img width="1386" alt="Screenshot 2024-06-26 at 1 54 09 PM" src="https://github.com/AshleyLab/VxE_Map/assets/96602087/956a3e16-a8a3-4c5c-8070-8db5cdf986ec">

#### Installation of Minigene Pipeline
1. Log in to Sherlock and navigate to directory where results are to be stored i.e. target_dir
2. Copy minigene_pipeline directory into your target directory using the below command  
   ```  
   cp -r /oak/stanford/groups/euan/projects/variant_effect_mapping/minigene_splicing_assay /path/to/target_dir 

#### Setting Up Minigene Pipeline's Environment
3. Navigate into copied minigene_splicing_assay directory
   ```
   cd /path/to/target_dir/minigene_splicing_assay
4. Ensure that Python 3.12 is installed on your system. Install it otherwise before proceeding.
5. Create and activate minigene_env 
   ```
   python3.12 -m venv minigene_env       # create new python virtual env
   source minigene_env/bin/activate      # activate minigene_env
6. Install packages
   ```
   pip install --upgrade pip
   pip install -r requirements.txt

   # Install Samtools (v1.16.1)
   wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
   tar -xvjf samtools-1.16.1.tar.bz2 
   cd samtools-1.16.1/
   ./configure --prefix=/path/to/target_dir/minigene_splicing_assay/minigene_splicing_assay/minigene_env --disable-lzma
   make
   make install
   cd ..                                                                  
   samtools --version                    # affirm that samtools was installed properly      

   # Install Chopper (v0.7.0)
   wget https://github.com/wdecoster/chopper/releases/download/v0.7.0/chopper-musl.zip
   unzip chopper-musl.zip -d chopper_bin
   mv chopper_bin/chopper minigene_env/bin
   chmod +x minigene_env/bin/chopper
   chopper --version                     # affirm that chopper was installed properly

   # Install minimap2 (2.28 (r1209))
   wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2
   tar -xvf minimap2-2.28_x64-linux.tar.bz2
   mv minimap2-2.28_x64-linux/minimap2 minigene_env/bin
   minimap2 --version                    # affirm that minimap2 was installed properly

   # Clean up
   rm samtools-1.16.1.tar.bz2 chopper-musl.zip minimap2-2.28_x64-linux.tar.bz2
   rm -r samtools-1.16.1 chopper_bin minimap2-2.28_x64-linux

#### Execution
7. Edit minigene_config.yaml file to suit your data
   ```
   vim minigene_config.yaml
   
   # press 'i' to edit and make changes; when done, press 'esc', ':wq' and hit 'enter'
8. Navigate into log directory
   ```
   cd log

   # full path of log dir = "/path/to/target_dir/minigene_splicing_assay/log
9. Execute minigene pipeline
   ```
   # For analysis of only barcode-SNV relationship, do:
   sbatch /path/to/target_dir/minigene_splicing_assay/minigene_pipeline.sbatch /path/to/target_dir/minigene_splicing_assay "" false

   # For analysis of barcode-codon relationship, do:
   sbatch /path/to/target_dir/minigene_splicing_assay/minigene_pipeline.sbatch /path/to/target_dir/minigene_splicing_assay barcode_codon_rs.py

   # For full analysis including isoform processing, do:
   sbatch /path/to/target_dir/minigene_splicing_assay/minigene_pipeline.sbatch /path/to/target_dir/minigene_splicing_assay
   
#### Tools:
samtools: [samtools](https://github.com/samtools/samtools)  
chopper: [chopper](https://github.com/wdecoster/chopper)
minimap2: [minimap2](https://github.com/lh3/minimap)  
deepvariant: [deepvariant](https://github.com/google/deepvariant)
