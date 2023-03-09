## Snakemake pipeline to process Metagenomic sequencing data and do straintracking
The pipelline generates the following things:
- Get a relative abundance table of species genome bins (SGBs) present in each sample (MetaPhlAn 4)
- Get stratified and unstratified pathway abundances and gene families (HUMAnN 3.6)
- Get information on the strain-sharing between samples (StrainPhlAn 4)

## In order to run the pipeline the following files/maps should be present:
### Conda environments
- mg4.yml

### Scripts:
- humann_strainphlan.snmk.py
- concat_strainsharing_thresholds.py 
	- Forked from merged_metaphlan_tables.py: https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/merge_metaphlan_tables.py)
- strainphlan_calc_threshold.R 
	- Forked from the tutorial of Valles-Colomer: https://github.com/biobakery/MetaPhlAn/wiki/Strain-Sharing-Inference / https://doi.org/10.1038/s41586-022-05620-1
- pyphlan scripts (see step 3 in the procedure of setting up the pipeline for the first time)
 
 
### Databases:
- hg19 bowtie2 index (see step 2 of the procedure; ~3.4 gb) 
	- Used to filter out human reads
- metaphlan database (able to download via metaphlan --force_download and metaphlan --install) 
	- For more info see https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4
- Chocophlan database (full_chocophlan.v201901_v31; see step 5 in the procedure)
- Protein database (uniref90; see step 5 in the procedure) 
	- For more info: https://github.com/biobakery/humann

### Other file:
- A dataframe (in my case called "df_for_thresholds.csv")
	- The first column of this file needs to be your samples (for me this was named "SAMPLE_ID")
	- For me, the dataframe looks as following:
  - ![image](https://user-images.githubusercontent.com/57362809/223940063-489ba2ef-992f-43ef-b5f5-af044a9ed65b.png)
	- Within the Rscript, this file is used to determine the relatedness to define the best threshold for normalised phylogenetic distance matrix
		- I used this pipeline for a Fecal Microbiota Transplantation (FMT) project
		- Based on your type of data (two timepoints without an FMT for instance), the R-script needs to be adjusted a bit

## Procedure when setting up for the first time: 
1. Install the conda environment:
	- "conda env create -f mg4.yml" Note that if there are conflicts found when installing the conda environments, you will have to set the channel_priority to false:
"conda config --set channel_priority false"

2. Download the "hg19" bowtie2 database:
	- https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	- Make sure the path is set correctly in the snakefile 
		- "humann_strainphlan.snmk.py" line 65; I used path "~/bowtie2_indexes/hg19/hg19"

3. Download pyphlan in the scripts folder:
	- Locate to the map where this repository is downloaded to
	- "mkdir scripts/pyphlan | git clone https://github.com/SegataLab/pyphlan scripts/pyphlan"

4. If you have not done so yet, download the metaphlan database:
	- "conda activate mg4"
	- "metaphlan --install" Note that I used metaphlan v4.0.5. at the time of creating the pipeline (database: mpa_vJan21_CHOCOPhlAnSGB_202103)

5. Download the nucleotide and protein databases for HUMAnN:
	- "humann_databases --download chocophlan full $INSTALL_LOCATION"
	- "humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION"

6. Set the configuration file to the locations of the downloaded databases:
	- "humann_config --update database_folders nucleotide $DIR_chocohplanfull"
	- "humann_config --update database_folders protein $DIR_uniref90"

7. Once everything is set-up, you can run the snakemake file:
	- Activate the conda environment (if you haven't done so already): "conda activate mg4"
	- "snakemake -s humann_strainphlan.snmk.py --cores 64"
		- Make sure that the "df_for_thresholds.csv" or its equivalent is present in the map
		- Make sure to have edited the path to the raw fastq files
		- lines 31 and 32 of the snakefile ("humann_strainphlan.snmk.py")
