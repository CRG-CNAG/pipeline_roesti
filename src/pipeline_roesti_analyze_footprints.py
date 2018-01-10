# -*- coding: utf-8 -*-
import argparse
import re
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("darkgrid")
seaborn.set_style("whitegrid")

parser = argparse.ArgumentParser()
parser.add_argument("bed_footprints_filename")
args = parser.parse_args()
print(args)
print("bed_footprints_filename:",args.bed_footprints_filename)
MPN_RP4_P2_FP_13860_ACTGAT.footprints_wide.bed
footprints_name = re.search(r'[/^](.+)\.bed', args.bed_footprints_filename).group(1)
print("footprints_name",footprints_name)

#bed_footprints_narrow_filename = re.sub(r'\.rRNA_removed\.bed', r'.footprints_narrow.bed', args.bed_filename)
#print("BED footprints narrow filename:",bed_footprints_narrow_filename)
#bed_footprints_wide_filename = re.sub(r'\.rRNA_removed\.bed', r'.footprints_wide.bed', args.bed_filename)
#print("BED footprints wide filename:",bed_footprints_wide_filename)

# Import all genes annotations


# Select all footprints that overlap with the CDS
bedtools intersect -wa -s -a MPN_RP4_P2_FP_13860_ACTGAT.footprints_wide.bed -b mpn_CDS_1.bed


