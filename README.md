# equine_mirna_project_pennvet
A collction of scripts and workflow for studying miRNA in equine


![Image](https://github.com/user-attachments/assets/7f5013cd-c327-4985-b216-142a72bfe90a)


## Description
This repo houses the code used for analyzing miRNA from Dr.Kyla Ortved's Equine extracellular vesicles(EV) project

## Abstract
(PLACEHOLDER)

## Dependencies (Tested working for this pipeline)
- Java (v18~v24)
- Docker (v27.5.1)
- Nextflow (v24.10.5)
- R (v4.5.0)
- ggplot2 (v3.5.2)
- tidyverse (v2.0.0)
- isomiRs (v1.32.1)
- pheatmap (v1.0.13)
- Python (v3.10.11)
- Bedtools (v2.31.0)

## Directory Structure

```
equine_mirna_project_pennvet/
├── project			
│   ├── samplesheet.csv
│   ├── contrasts.csv
│  └── smallrnaseq_samplesheet_4da.csv
├── project_2/
│   ├── samplesheet_2.csv
│   ├── contrasts_2.csv
│  └── smallrnaseq_samplesheet_4da_2.csv
├── READMd	# This file
├─LICENSE	# License
├── nfcore_smrnaseq.sh	# Main script for running smrnaseq pipeline
├── nfcore_smrnaseq.sh	# Main script for running differentialabundance pipeline
├── process_mirdeep2.py	# Script for processing per sample miRDeep2 results
├─volcano_plotting.R	# Script for plotting differentally expressed genes
└── post_anals.sh	# Script for processing isoMir output

```

## Credits

- Mana Okudaira,Lauren K Olenick,Alexandra IJ Usimaki,Hoda Elkhenany,Jillian Bastidas,Renata L Linardi,Angela M Gaesser,Shannon S Connard,Luca Musante,Rui Xiao,Daniel Beiting,Kyla F Ortved
- The pipeline is developed and implemented by Rui Xiao

## Data availability
- The small RNAseq data generated for this study is available through (Bioproject place holder).


## Citations
(PLACE HOLDER FOR MANUSCRIPT DOI)

