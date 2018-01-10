#!/bin/bash

GENOME_PATH="/users/lserrano/mweber/Translation_model/Ribosome_profiling_data/bowtie2_indexed_genome/Ecoli"
GENOME_NAME="ecoli_K12_MG1655_NC_000913.3"

nohup /users/lserrano/mweber/Research_Dropbox/RNA-seq_data_analysis/src/pipeline_roesti_3.py \
    --fastq-files *read[12].fastq* \
    --indexed-ref-genome "${GENOME_PATH}/${GENOME_NAME}" \
    --rRNA-bedfile "${GENOME_PATH}/${GENOME_NAME}_rRNA.bed" \
    --genome-bedfile "${GENOME_PATH}/${GENOME_NAME}.genome" \
    --genome-CDS-bedfile "${GENOME_PATH}/${GENOME_NAME}_CDS.bed" \
    --seq-end paired-end \
    --library-type rna-seq \
    --pipeline-name roesti_3_PE &> run_PE.out&
