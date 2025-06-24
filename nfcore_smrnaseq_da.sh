#!/usr/bin/env bash

# Usage:
#   chmod u+x nfcore_smrnaseq_da.sh
#   ./nfcore_smrnaseq_da.sh
#   Make sure to have proper mirTOP output files in place
#	--deseq2_vst_nsub 100 somce only ~400 total miRNA annotated in equine genome 3.0


nextflow run nf-core/differentialabundance \
    --input smallrnaseq_samplesheet_4da.csv \
    --study_name project_mirTOP_DA \
    --contrasts conctasts.csv \
    --matrix mirtop_counts.tsv \
    --outdir project/DA_mirTOP  \
    --gtf eca.gtf \
    -profile rnaseq,docker \
    -c nfcore.config \
    --deseq2_vst_nsub 100