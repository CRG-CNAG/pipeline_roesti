#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("darkgrid")
seaborn.set_style("whitegrid")

parser = argparse.ArgumentParser()
parser.add_argument("input_path")
parser.add_argument("output_path")
parser.add_argument("bed_filename")
parser.add_argument("sample_name")

parser.add_argument("--footprint-min-size", dest="footprint_min_size", type=int, default=28,
                    help="footprint minimum size (inclusive)")
parser.add_argument("--footprint-max-size", dest="footprint_max_size", type=int, default=30,
                    help="footprint maximum size (inclusive)")
parser.add_argument("--footprint-min-size-wide", dest="footprint_min_size_wide", type=int, default=20,
                    help="footprint minimum size (inclusive)")
parser.add_argument("--footprint-max-size-wide", dest="footprint_max_size_wide", type=int, default=40,
                    help="footprint maximum size (inclusive)")
args = parser.parse_args()
print(args)
bed_filename = args.bed_filename
input_path = args.input_path
bed_filepath = str(Path(input_path) / bed_filename)
output_path = args.output_path
sample_name = args.sample_name
print("bed_filepath:", bed_filepath)
print("sample_name:", sample_name)

bed_footprints_narrow_filename = str(Path(output_path) / (sample_name + '.footprints_narrow.bed'))
print("BED footprints narrow filename:",bed_footprints_narrow_filename)
bed_footprints_wide_filename = str(Path(output_path) / (sample_name + '.footprints_wide.bed'))
print("BED footprints wide filename:",bed_footprints_wide_filename)

i = 0
size_data = []
with open(bed_filepath,'r') as bed_file:
    with open(bed_footprints_narrow_filename,'w') as bed_footprints_narrow_file:
        with open(bed_footprints_wide_filename,'w') as bed_footprints_wide_file:
            for line in bed_file:
                ref, start, end, name, mapq, strand = line.rstrip().split('\t')
                size = int(end) - int(start)
                size_data.append((size,strand))
                if args.footprint_min_size <= size <= args.footprint_max_size:
                    i += 1
                    bed_footprints_narrow_file.write(line)
                if args.footprint_min_size_wide <= size <= args.footprint_max_size_wide:
                    bed_footprints_wide_file.write(line)
                #line2 = "\t".join((chr1, insert_start, insert_end, name, mapq, strand)) + '\n'
                #if i == 100:
                    #break

# Plot distribution of insert sizes
sizeDf = pd.DataFrame(size_data, columns=['insert_size','strand'])
print("Total nb of inserts: ",len(sizeDf))


#hist, bins = np.histogram(size_data, list(range(0,500)))
#hist_filename = re.sub(r'\.rRNA_removed\.bed', r'.insert_size_dist.csv', args.bed_filename)
#pd.Series(hist, bins[:-1]).to_csv(hist_filename)
#x = [(a+bins[i+1])/2.0 for i,a in enumerate(bins[0:-1])]
#print(x)
#hist = pd.Series(hist, x)
#print(hist)

# Write the insert size distribution in the input path

plot_filename = str( Path(output_path) / (sample_name + '.insert_size_dist.png'))
print(plot_filename)
fig = plt.figure(figsize=(14,12))
ax = fig.add_subplot(111)
#hist.plot.bar(ax=ax, label='insert_size')
sizeDf.hist(by='strand', ax=ax, alpha=1, bins=list(range(10,60)), layout=(2,1))
#sizeDf.groupby('strand').hist(by='strand', alpha=0.5, bins=list(range(0,50)))
ax.set_title(sample_name + ' insert size')
plt.savefig(plot_filename)
