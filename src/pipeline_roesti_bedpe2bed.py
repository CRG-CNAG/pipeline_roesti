#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("bedpe_filename")
args = parser.parse_args()
print("BEDPE filename:",args.bedpe_filename)
bed_filename = re.sub(r'\.bedpe', r'.bed', args.bedpe_filename)
print("BED filename:",bed_filename)

with open(args.bedpe_filename,'r') as bedpe_file:
    with open(bed_filename,'w') as bed_file:
        for line in bedpe_file:
            chr1, start1, end1, chr2, start2, end2, name, mapq, strand = line.rstrip().split('\t')
            start1 = int(start1)
            start2 = int(start2)
            end1 = int(end1)
            end2 = int(end2)
            insert_start = start1 if start1 <= start2 else start2
            insert_end = end1 if end1 >= end2 else end2
            line2 = "\t".join((chr1, str(insert_start), str(insert_end), name, mapq, strand)) + '\n'
            bed_file.write(line2)
