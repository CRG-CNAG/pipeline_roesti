#!/bin/bash

if [ $# -lt 1 ]
  then
    echo "ERROR: command line contains $# arguments should be 1"
    exit 1
fi

GENOME_CDS_BED_FILE=$1
echo ${GENOME_CDS_BED_FILE}
NTHREADS=1

mv ${GENOME_CDS_BED_FILE} ${GENOME_CDS_BED_FILE}.unsorted
while [ ! -f ${GENOME_CDS_BED_FILE}.unsorted ]
do
  sleep 2
done

/users/lserrano/mweber/local/bin/sort ${GENOME_CDS_BED_FILE}.unsorted -k 1,1 -k 2,2n -k 3,3n > ${GENOME_CDS_BED_FILE}

while [ ! -f ${GENOME_CDS_BED_FILE} ]
do
  sleep 2
done