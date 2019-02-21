#!/bin/bash

if [ $# != 7 ]; then
    echo "ERROR: command line contains $# arguments should be 7"
    exit 1
fi

SAMPLE=$1
INPUT_PATH=$2
OUTPUT_PATH=$3
USE_LOCAL_DISK=$4 # true or false
SCRIPT_PATH=$5
GENOME_BED_FILE=$6
NTHREADS=$7
# SCRIPT_PATH="/users/lserrano/mweber/Research_Dropbox/Translational_regulation_model/RNA-seq_data_analysis"
# SAMPLE="MPN_RP4_P2_FP_13860_ACTGAT"
# INPUT_PATH_LOCAL="."
# OUTPUT_PATH_LOCAL="."
# USE_LOCAL_DISK=false

if [ "$USE_LOCAL_DISK" = true ] ; then
    echo "Using local disk for computations"
    # Copy input file to local disk on node
    #cp ${INPUT_PATH}/${SAMPLE}_sorted.bam ${TMPDIR}/${SAMPLE}_sorted.bam
    echo "Copy done."
    INPUT_PATH_LOCAL=${TMPDIR}
    OUTPUT_PATH_LOCAL=${TMPDIR}
else
    echo "Using normal filesystem for computations"
    INPUT_PATH_LOCAL=${INPUT_PATH}
    OUTPUT_PATH_LOCAL=${OUTPUT_PATH}
fi

# Select ribosome footprints based on insert size.
# Write the output to new BED files.
python3 ${SCRIPT_PATH}/pipeline_roesti_extract_footprints.py ${INPUT_PATH_LOCAL} ${OUTPUT_PATH_LOCAL} ${SAMPLE}.filtered.bed ${SAMPLE}

# Count nb of footprints reads
wc -l ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.bed > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.nreads
wc -l ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.bed > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.nreads

# Compute per base coverage of footprints
# we use 0-based index
bedtools genomecov -i ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.bed -g ${GENOME_BED_FILE} -d -strand + \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.strandp_coverage.bed
bedtools genomecov -i ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.bed -g ${GENOME_BED_FILE} -d -strand - \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_narrow.strandm_coverage.bed
bedtools genomecov -i ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.bed -g ${GENOME_BED_FILE} -d -strand + \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.strandp_coverage.bed
bedtools genomecov -i ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.bed -g ${GENOME_BED_FILE} -d -strand - \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.footprints_wide.strandm_coverage.bed

if [ "$USE_LOCAL_DISK" = true ] ; then
    echo "Copying files from node local disk back to filesystem"
    #cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_F1R2 ${OUTPUT_PATH}/
fi


