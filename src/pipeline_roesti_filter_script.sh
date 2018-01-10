#!/bin/bash

if [ $# != 10 ]; then
    echo "ERROR: command line contains $# arguments should be 10"
    exit 1
fi



SAMPLE=$1
INPUT_PATH=$2
OUTPUT_PATH=$3
QUALITY_THRESHOLD=$4
USE_LOCAL_DISK=$5 # true or false
RRNA_BEDFILE=$6
RRNA_TRNA_BEDFILE=$7
SCRIPT_PATH=$8
REMOVE_RRNA=$9
SEQ_END=${10}
#RRNA_BEDFILE="/users/lserrano/mweber/Research_Dropbox/Mycoplasma_pneumoniae_experimental_data/mpn_annotation_features/mpn_rRNA.bed"
#RRNA_TRNA_BEDFILE="/users/lserrano/mweber/Research_Dropbox/Mycoplasma_pneumoniae_experimental_data/mpn_annotation_features/mpn_rRNA_tRNA.bed"
#SCRIPT_PATH="/users/lserrano/mweber/Research_Dropbox/Translational_regulation_model/RNA-seq_data_analysis"
NTHREADS=8
# SAMPLE="MPN_RP4_P2_FP_13860_ACTGAT"
# INPUT_PATH_LOCAL="."
# OUTPUT_PATH_LOCAL="."
# USE_LOCAL_DISK=false

if [ "$USE_LOCAL_DISK" = true ] ; then
    echo "Using local disk for computations"
    # Copy input file to local disk on node
    echo "Copying input BAM file to local disk..."
    cp ${INPUT_PATH}/${SAMPLE}_sorted.bam ${TMPDIR}/${SAMPLE}_sorted.bam
    echo "Copy done."
    INPUT_PATH_LOCAL=${TMPDIR}
    OUTPUT_PATH_LOCAL=${TMPDIR}
else
    echo "Using normal filesystem for computations"
    INPUT_PATH_LOCAL=${INPUT_PATH}
    OUTPUT_PATH_LOCAL=${OUTPUT_PATH}
fi

# Filter samples by quality, keep only primary alignments
# -q INT        Drop alignments with MAPQ smaller than INT [0].
# -F 4          Drop unmapped reads
# -F 256        Keep only primary alignments
FILTER_FLAGS="-bu -F 4 -F 256 -q ${QUALITY_THRESHOLD}"
if [[ "$SEQ_END" = "single-end" ]]; then
    FILTER_FLAGS="${FILTER_FLAGS}"
elif [[ "$SEQ_END" = "paired-end" ]]; then
    # -f 2          Keep only properly paired-end reads
    FILTER_FLAGS="${FILTER_FLAGS} -f 2"
fi

samtools view ${FILTER_FLAGS} ${INPUT_PATH_LOCAL}/${SAMPLE}_sorted.bam \
| samtools sort -@ ${NTHREADS} -o ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam

samtools index ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam

# Check strand orientation for read 1 and read 2.
# Forward1-Reverse2 should correspond to strand -
# Reverse1-Forward2 should correspond to strand + (rRNA)
# Use high count of reads for rRNA region to determine strand orientation.

# Extract first rrna region from BED file
read REF_GENOME RRNA1_START RRNA1_END RRNA1_STRAND <<< "$(awk '(NR==1) {print $1, $2, $3, $6}' ${RRNA_BEDFILE})"

# counting reads in the first rRNA region, reads 1 forward
if [[ "$SEQ_END" = "single-end" ]]; then
    # -c             only count reads
    # -F 16          keep only reads that are NOT reverse complemented
    FILTER_FLAGS="-c -F 16"
elif [[ "$SEQ_END" = "paired-end" ]]; then
    # -c             only count reads
    # -F 16          keep only reads that are NOT reverse complemented
    # -f 64          keep only read 1 in read pair
    FILTER_FLAGS="-c -F 16 -f 64"
fi
samtools view ${FILTER_FLAGS} ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam ${REF_GENOME}:${RRNA1_START}-${RRNA1_END} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_F

# counting reads in the first rRNA region, reads 1 reverse
if [[ "$SEQ_END" = "single-end" ]]; then
    # -c             only count reads
    # -f 16          keep only reads that are reverse complemented
    FILTER_FLAGS="-c -f 16"
elif [[ "$SEQ_END" = "paired-end" ]]; then
    # -c             only count reads
    # -f 16          keep only reads that are reverse complemented
    # -f 64          keep only read 1 in read pair
    FILTER_FLAGS="-c -f 16 -f 64"
fi
samtools view ${FILTER_FLAGS} ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam ${REF_GENOME}:${RRNA1_START}-${RRNA1_END} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_R

# Import read counts and compare them
typeset -i nreads_rRNA1_F=$(cat ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_F)
typeset -i nreads_rRNA1_R=$(cat ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_R)
STRAND=0
if [[ "$RRNA1_STRAND" = "+" ]]; then
    if [ $nreads_rRNA1_F -ge $nreads_rRNA1_R ]
    then
        # Note: we also include here the case where the nb of reads are equal
        echo "sample: ${SAMPLE}  F1R2 is strand +, R1F2 is strand -"
        STRAND=+1
        STRAND_COLUMN=9
    elif [ $nreads_rRNA1_F -lt $nreads_rRNA1_R ]
    then
        echo "sample: ${SAMPLE}  R1F2 is strand +, F1R2 is strand -"
        STRAND=-1
        STRAND_COLUMN=10
    fi
elif [[ "$RRNA1_STRAND" = "-" ]]; then
    if [ $nreads_rRNA1_F -ge $nreads_rRNA1_R ]
    then
        # Note: we also include here the case where the nb of reads are equal
        echo "sample: ${SAMPLE}  F1R2 is strand -, R1F2 is strand +"
        STRAND=-1
        STRAND_COLUMN=10
    elif [ $nreads_rRNA1_F -lt $nreads_rRNA1_R ]
    then
        echo "sample: ${SAMPLE}  R1F2 is strand -, F1R2 is strand +"
        STRAND=+1
        STRAND_COLUMN=9
    fi

fi
echo ${STRAND} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.F.strand


if [[ "$SEQ_END" = "single-end" ]]; then
    # Convert BAM file reads to BED format (single-end reads).
    samtools view -h ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam \
    | samtools sort -n -l 0 \
    | bedtools bamtobed -i stdin > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed

elif [[ "$SEQ_END" = "paired-end" ]]; then
    # Convert BAM file reads to BEDPE format (BED format for paired end)
    # BAM file has to be sorted by read *name* in order for bedtools to find read pairs (normal sorting is by chromosomal position)
    # Then we cut the BEDPE format to form a simple BED format, keeping the start pos of read1 and end pos of read2.
    # Note that bedtools converts forward reads to + strand and reverse reads to - strand.
    # In paired-end reads, what determines the strand is the order of F1R2 or R1F2.
    # We use the option for bamtobed: -mate1 When writing BEDPE (-bedpe) format, always report mate one as the first BEDPE “block”.
    # This way, we will have in BEDPE file:
    # column 9 = strand read 1, column 10 = strand read 2
    # Note that in the BEDPE file, read 1 may have the start position after the read 2, we will have to parse
    # the file to convert the start and end positions from BEDPE to BED.
    samtools view -h ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bam \
    | samtools sort -n -l 0 \
    | bedtools bamtobed -i stdin -bedpe -mate1 \
    | cut -f 1,2,3,4,5,6,7,8,${STRAND_COLUMN} \
    > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bedpe

    # Compute start and end position for inserts and convert BEDPE to BED file
    python3 ${SCRIPT_PATH}/pipeline_roesti_bedpe2bed.py ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bedpe
fi


# Remove aligned reads that overlap with the rRNA or tRNA regions.
if [ "$REMOVE_RRNA" = true ] ; then
    bedtools intersect -v -s -a ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed -b ${RRNA_TRNA_BEDFILE} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed
else
    mv ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed
fi

# Count nb of reads
if [ "$REMOVE_RRNA" = true ] ; then
    wc -l < ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed.nreads
else
    wc -l < ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed.nreads
fi
wc -l < ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed.nreads

# Sort BED file
# Sort is VERY SLOW FOR LARGE FILES. When we remove rRNA, the file /${SAMPLE}.filtered.bed is much reduced to ~40M,
# so in this case sorting take litlle time.
# We need to sort the BED file, otherwise bedtools intersect will use huge memory (70G) when computing the number of
# reads per gene and is much slower.
# We can use the optimized and parallel implementation of sort provided by GNU coreutils.
# I installed it at /users/lserrano/mweber/local/bin/sort
# See http://davetang.org/muse/2013/11/20/sorting-a-huge-bed-file/

# Testing on file /users/lserrano/mweber/Translation_model/Ribosome_profiling_data/Mpn/RNA-seq/roesti_3_PE/Task05_filter_alignments/MPN_RP10_P2_UD_15116_AGTTCC.filtered.bed
# /usr/bin/time -v head -100000000 MPN_RP10_P2_UD_15116_AGTTCC.filtered.bed | sort -k 2,2n -k 3,3n > MPN_RP10_P2_UD_15116_AGTTCC.filtered.bed.sorted
# /usr/bin/time -v head -100000000 MPN_RP10_P2_UD_15116_AGTTCC.filtered.bed | /users/lserrano/mweber/local/bin/sort -k 2,2n -k 3,3n --parallel=2 > MPN_RP10_P2_UD_15116_AGTTCC.filtered.bed.sorted
# For 1e8 lines, system sort, 1 thread, 11mn28s, 2.7M peak_mem
# For 1e8 lines, sort (GNU coreutils) 8.4, 8 threads, 1mn5s, 2.7M peak_mem
# For 1e8 lines, sort (GNU coreutils) 8.4, 16 threads, 0mn54s, 2.7M peak_mem
/users/lserrano/mweber/local/bin/sort ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed -k 1,1 -k 2,2n -k 3,3n --parallel=${NTHREADS} > ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed.sorted

while [ ! -f ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed.sorted ]
do
  sleep 2
done

mv ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed.sorted ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed

if [ "$USE_LOCAL_DISK" = true ] ; then
    echo "Copying files from node local disk back to filesystem"
    cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_F ${OUTPUT_PATH}/
    cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.nreads_rRNA1_R ${OUTPUT_PATH}/
    cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.F.strand ${OUTPUT_PATH}/
    cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.filtered.bed ${OUTPUT_PATH}/
    if [ "$REMOVE_RRNA" = true ] ; then
        cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bed ${OUTPUT_PATH}/
    fi
    cp ${OUTPUT_PATH_LOCAL}/${SAMPLE}.bedpe ${OUTPUT_PATH}/
fi
