#!/usr/bin/env bash

# Usage:
#   chmod u+x nfcore_smrnaseq.sh
#   ./nfcore_smrnaseq.sh
#   Make sure to have proper fastq files in place


nextflow run nf-core/smrnaseq \
    -r dev \
    --input sampleshee.csv \
    --save_reference \
    --outdir project \
    --filter_contamination \
    --mirgenedb \
    --fasta ec3_genome.fa \
    --cdna ec3_cDNA.fas \
    --ncrna ec3_lncRNA.fas \
    --mirgenedb_species "Eca" \
    --mirgenedb_gff eca.gff \
    --mirgenedb_mature eca_mature.fas \
    --mirgenedb_hairpin eca_precursor.fas \
    -profile docker,illumina \
    -bg