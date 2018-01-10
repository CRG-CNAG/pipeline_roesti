#!/bin/bash

# Filter bed file by feature size MINSIZE <= L <= MAXSIZE

BEDFILENAME=$1
MINSIZE=$2
MAXSIZE=$3

awk -v a="${MINSIZE}" -v b="${MAXSIZE}" '(a <= ($3-$2)) && (($3-$2) <= b)' "${BEDFILENAME}"
