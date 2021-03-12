#!/bin/bash
if [ $# -lt 7 ]
  then
    echo "ERROR: command line contains $# arguments should be 7 or 8"
    exit 1
fi

BED_FILE=$1
SAMPLE=$2
INPUT_PATH=$3
OUTPUT_PATH=$4
USE_LOCAL_DISK=$5 # true or false
SCRIPT_PATH=$6
GENOME_BED_FILE=$7
GENOME_CDS_BED_FILE=$8
NREADS_BED=$9
RRNA_BED_FILE=$10
NTHREADS=1

# Test variables
# $BED_FILE $SAMPLE $INPUT_PATH $OUTPUT_PATH $USE_LOCAL_DISK $SCRIPT_PATH bedfile.bed
# SAMPLE="test"
# INPUT_PATH="Task05_filter_alignments"
# OUTPUT_PATH="Task07_genome_coverage_fragment_count"
# BED_FILE=${INPUT_PATH}/"test.rRNA_removed.bed"
# GENOME_CDS_BED_FILE="NC_000912.1/NC_000912.1_CDS.bed"
# RRNA_BED_FILE="NC_000912.1/NC_000912.1_rRNA.bed"
# USE_LOCAL_DISK=false
# SCRIPT_PATH="/users/lserrano/mweber/Research_cloud/RNA-seq_data_analysis/src"
# OUTPUT_PATH_LOCAL=${OUTPUT_PATH}

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

# Run only if argument GENOME_CDS_BED_FILE is given
if [[ ! -z "$8" ]]; then

    # Compute fragments count for each CDS
    # -c    For each entry in A, report the number of hits in B while restricting to -f.
    #       Reports 0 for A entries that have no overlap with B. Restricted by -f and -r.
    # -s    Force “strandedness”. That is, only report hits in B that overlap A on the same strand.
    #       By default, overlaps are reported without respect to strand.
    # -F    Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
    echo "Computing count per CDS."
    bedtools intersect -c -s -F 0.5 -sorted -g ${GENOME_BED_FILE} -a ${GENOME_CDS_BED_FILE} -b ${BED_FILE} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.CDS_fragment_count.bed
fi

if [[ ! -z "$9" ]]; then
    echo "Computing count per CDS for rRNA."
    bedtools intersect -c -s -F 0.5 -sorted -g ${GENOME_BED_FILE} -a ${RRNA_BED_FILE} -b ${BED_FILE} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.rRNA_fragment_count.bed
fi

# Compute per base coverage
echo "Computing per base coverage for plus strand."
bedtools genomecov -i ${BED_FILE} -g ${GENOME_BED_FILE} -d -strand + \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.strandp_coverage.bed
echo "Computing per base coverage for minus strand."
bedtools genomecov -i ${BED_FILE} -g ${GENOME_BED_FILE} -d -strand - \
> ${OUTPUT_PATH_LOCAL}/${SAMPLE}.strandm_coverage.bed

# Compute total integral of per-base coverage for each CDS
# ... python script
echo "Computing CDS read count measures."
python3 "${SCRIPT_PATH}/pipeline_roesti_mean_coverage_per_CDS.py" \
    "${OUTPUT_PATH_LOCAL}/${SAMPLE}.strandp_coverage.bed" \
    "${OUTPUT_PATH_LOCAL}/${SAMPLE}.strandm_coverage.bed" \
    "${GENOME_CDS_BED_FILE}" \
    "${OUTPUT_PATH_LOCAL}" \
    "${SAMPLE}" \
    "${NREADS_BED}" \
    "${RRNA_BED_FILE}"

if [ "$USE_LOCAL_DISK" = true ] ; then
    echo "Copying files from node local disk back to filesystem"
    #cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_F1R2 ${OUTPUT_PATH}/
fi


