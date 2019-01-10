#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("whitegrid")



parser = argparse.ArgumentParser()
parser.add_argument("bed_cov_strand_p_file")
parser.add_argument("bed_cov_strand_m_file")
parser.add_argument("cds_bed_file")
parser.add_argument("output_path")
parser.add_argument("sample_name")
parser.add_argument("rRNA_bedfile", default=None)
args = parser.parse_args()

verbose = 2
bedCols = ['ref', 'start', 'end', 'name', 'mapq', 'strand', 'value']

test = False
if test:
    args.bed_cov_strand_p_file = '/users/lserrano/mweber/Ribosome_profiling_data/Ecoli/RNA-seq/roesti_3_PE/Task07_genome_coverage_fragment_count/Eco_RP4_S1_13104_GCCAAT.strandp_coverage.bed'
    args.bed_cov_strand_m_file = '/users/lserrano/mweber/Ribosome_profiling_data/Ecoli/RNA-seq/roesti_3_PE/Task07_genome_coverage_fragment_count/Eco_RP4_S1_13104_GCCAAT.strandm_coverage.bed'
    args.cds_bed_file = '/users/lserrano/mweber/Ribosome_profiling_data/bowtie2_indexed_genome/Ecoli/ecoli_K12_MG1655_NC_000913.3_CDS.bed'
    args.output_path = '/users/lserrano/mweber/Ribosome_profiling_data/Ecoli/RNA-seq/roesti_3_PE/Task07_genome_coverage_fragment_count'
    args.sample_name = 'Eco_RP4_S1_13104_GCCAAT'


method = 'memory'

print("Computing average coverage per base for all CDS...")
outputFilename = str(Path(args.output_path) / "{}.CDS_average_coverage.bed".format(args.sample_name))

if method == 'file':
    # WARNING: RRNA NOT IMPLEMENTED IN FILE MODE !!!!!!!!!!!!!!!!!
    with open(outputFilename, 'w') as avgCovFile:
        with open(args.cds_bed_file, 'r') as cdsFile:

            for lineCDS in cdsFile:
                sumCovInCDS = np.float64(0)
                refCDS, start, end, name, mapq, strand = lineCDS.rstrip().split('\t')
                if verbose >= 2: print(name)
                start = int(start)
                end = int(end)
                CDSlength = end - start

                if strand == '+':
                    covFileName = args.bed_cov_strand_p_file
                elif strand == '-':
                    covFileName = args.bed_cov_strand_m_file

                with open(covFileName, 'r') as covFile:
                    for lineCov in covFile:
                        refCov, pos, cov = lineCov.rstrip().split('\t')
                        pos = int(pos)
                        if refCov == refCDS and start <= pos < end:
                            sumCovInCDS += np.float64(cov)
                        elif refCov == refCDS and pos >= end:
                            # Assume that both coverage files are sorted by position
                            break
                    avgCovInCDS = sumCovInCDS / CDSlength
                    # print("name", name, "start", start, "end", end, "strand", strand, "sumCovInCDS", sumCovInCDS,
                    #       "CDSlength", CDSlength, "avgCovInCDS", avgCovInCDS)

                    avgCovFile.write("\t".join((refCDS, str(start), str(end), name, mapq, strand, str(avgCovInCDS))) + '\n')

elif method == 'memory':

    covColList = ['ref', 'pos', 'cov']
    covPlusDf = (pd.read_table(args.bed_cov_strand_p_file, header=None, names=covColList)
                 .sort_values('pos'))
    covMinusDf = (pd.read_table(args.bed_cov_strand_m_file, header=None, names=covColList)
                  .sort_values('pos'))
    if Path(args.cds_bed_file).is_file():
        CDSDf = (pd.read_table(args.cds_bed_file, header=None, names=bedCols)
                 .sort_values(['start', 'end']))
        CDSDf['length'] = CDSDf['end'] - CDSDf['start']
    else:
        CDSDf = pd.DataFrame(columns=bedCols + ['length'])

    if Path(args.rRNA_bedfile).is_file():
        rRNADf = (pd.read_table(args.rRNA_bedfile, header=None, names=bedCols)
                  .sort_values(['start', 'end']))
        rRNADf['length'] = rRNADf['end'] - rRNADf['start']
        rRNAList = rRNADf['name'].tolist()
    else:
        rRNADf = pd.DataFrame(columns=bedCols + ['length'])
        rRNAList = []
    CDSDf = pd.concat([CDSDf, rRNADf])


    with open(outputFilename, 'w') as avgCovFile:
        for index, cds in CDSDf.iterrows():
            if verbose >= 2: print(cds['name'])
            # ['ref', 'start', 'end', 'name', 'mapq', 'strand', 'value']
            start = int(cds['start']) if not pd.isnull(cds['start']) else None
            end = int(cds['end']) if not pd.isnull(cds['end']) else None
            CDSlength = int(cds['length']) if not pd.isnull(cds['length']) else None

            if start is not None and end is not None and CDSlength is not None:
                if cds['strand'] == '+':
                    covDf = covPlusDf
                elif cds['strand'] == '-':
                    covDf = covMinusDf

                avgCovInCDS = covDf[(start <= covDf['pos']) & (covDf['pos'] < end)]['cov'].sum() / CDSlength

                avgCovFile.write("\t".join((cds['ref'], str(start), str(end), str(cds['name']),
                                            str(cds['mapq']), cds['strand'], str(avgCovInCDS))) + '\n')


print("Computing average coverage per base for all CDS... finished")



print("Plotting and saving average coverage per base and fragment count...")

avgCovDf = pd.read_table(outputFilename, header=None, names=bedCols)
print("avgCovDf\n", avgCovDf)
if len(avgCovDf) > 0:
    avgCovDf.rename(columns={'value':'avg_coverage'}, inplace=True)

    fragmentCountDf = pd.read_table(str(Path(args.output_path) / "{}.CDS_fragment_count.bed".format(args.sample_name)),
                                    header=None, names=bedCols)
    fragmentCountDf.rename(columns={'value':'fragment_count'}, inplace=True)
    if len(fragmentCountDf) > 0:
        fragmentCountDf['fragment_count_per_base'] = fragmentCountDf.apply(lambda x:
                                                                           (x['fragment_count'] /
                                                                            (x['end'] - x['start'])),
                                                                           axis=1)
    # print(avgCovDf.head())
    # print(fragmentCountDf.head())

    filename = str(Path(args.output_path) / "{}.rRNA_fragment_count.bed".format(args.sample_name))
    if Path(filename).exists():
        rRNAfragmentCountDf = pd.read_table(filename, header=None, names=bedCols)
        rRNAfragmentCountDf.rename(columns={'value':'fragment_count'}, inplace=True)
        if len(rRNAfragmentCountDf) > 0:
            rRNAfragmentCountDf['fragment_count_per_base'] = rRNAfragmentCountDf.apply(lambda x:
                                                                                       (x['fragment_count'] /
                                                                                        (x['end'] - x['start'])),
                                                                                       axis=1)
        fragmentCountDf = pd.concat([fragmentCountDf, rRNAfragmentCountDf])
    else:
        rRNAfragmentCountDf = None


    dataDf = pd.merge(avgCovDf, fragmentCountDf, on=['ref', 'name', 'start', 'end', 'strand', 'mapq'], how='outer')
    dataDf.rename(columns={'name':'id'}, inplace=True)
    dataDf.sort_index(inplace=True)

    # We will compute also the same measure but removing all the rRNA counts
    dataDf_noRRNA = dataDf[dataDf['id'].map(lambda x: x not in rRNAList)].copy()
    print("dataDf_noRRNA.columns", dataDf_noRRNA.columns)

    # Compute fragment per kilobase per million fragments (FPKM)
    nFragmentTotal = dataDf['fragment_count'].sum()
    dataDf['FPKM'] = dataDf['fragment_count']/((nFragmentTotal/1e6) * ((dataDf['end'] - dataDf['start'])/1e3))

    nFragmentTotal = dataDf_noRRNA['fragment_count'].sum()
    dataDf_noRRNA['FPKM_no_rRNA'] = dataDf_noRRNA['fragment_count']/((nFragmentTotal/1e6) * ((dataDf_noRRNA['end'] - dataDf_noRRNA['start'])/1e3))

    # Compute transcript per million (TMP)
    # TPM is just a way of normalizing by the total amount of transcripts.
    # Here by transcript we mean the CDS regions. We do not consider polycistronic
    # transcripts.

    # The usual way is to use count per base to derive TPM.
    fragmentPerBase = dataDf['fragment_count']/(dataDf['end'] - dataDf['start'])
    fragmentPerBaseTot = fragmentPerBase.sum()
    dataDf['TPM (based on fragment count)'] = (1e6 * fragmentPerBase / fragmentPerBaseTot)

    fragmentPerBase = dataDf_noRRNA['fragment_count']/(dataDf_noRRNA['end'] - dataDf_noRRNA['start'])
    fragmentPerBaseTot = fragmentPerBase.sum()
    dataDf_noRRNA['TPM (based on fragment count) no_rRNA'] = (1e6 * fragmentPerBase / fragmentPerBaseTot)

    # We also use average coverage to derive TPM.
    avgCovTot = dataDf['avg_coverage'].sum()
    dataDf['TPM (based on average coverage)'] = (1e6 * dataDf['avg_coverage'] / avgCovTot)

    avgCovTot = dataDf_noRRNA['avg_coverage'].sum()
    dataDf_noRRNA['TPM (based on average coverage) no_rRNA'] = (1e6 * dataDf_noRRNA['avg_coverage'] / avgCovTot)

    colList = ['ref', 'id', 'start', 'end', 'strand', 'mapq'] + list(set(dataDf_noRRNA.columns) - set(dataDf.columns))
    dataDf = pd.merge(dataDf, dataDf_noRRNA[colList], on=['ref', 'id', 'start', 'end', 'strand', 'mapq'], how='outer')

    # We switch to 1-based index
    dataDf['start'] = dataDf['start'] + 1
    dataDf['end'] = dataDf['end'] + 1
    dataDf.rename(columns={'start':'start_1-based', 'end':'end_1-based'}, inplace=True)

    print("dataDf.columns", dataDf.columns)
    dataDf.to_csv(str(Path(args.output_path) / "{}.CDS_values.csv".format(args.sample_name)))
    dataDf['fragment_count < 100'] = dataDf['fragment_count'] < 100



    # Draw correlation
    plot_filename = str(Path(args.output_path) / "{}.CDS_fragment_count_avg_coverage_corr.png".format(args.sample_name))
    fig = plt.figure(figsize=(14,12))
    ax = fig.add_subplot(111)
    ax.set_title(args.sample_name)
    # Drop the ribosomal RNA which would mask the real correlation beteween CDS measures
    plotDf = dataDf[dataDf['id'].map(lambda x: x not in rRNAList)].copy()
    if len(plotDf) > 5:
        seaborn.lmplot(x="fragment_count_per_base", y="avg_coverage", data=plotDf)
        plt.savefig(plot_filename)

        plot_filename = str(Path(args.output_path) / "{}.CDS_normalization_methods_corr.png".format(args.sample_name))
        fig, ax = plt.subplots(figsize=(22, 22))
        ax.set_title(args.sample_name)

        seaborn.set(font_scale=0.7)
        grid = seaborn.pairplot(data=plotDf,
                                vars=['avg_coverage', 'fragment_count', 'FPKM',
                                      'TPM (based on fragment count)', 'TPM (based on average coverage)'],
                                hue='fragment_count < 100', kind='reg')
        plt.savefig(plot_filename)

    print("Plotting and saving average coverage per base and fragment count... finished")
else:
    print("No CDS annotation, writing empty CDS_values file.")
    # TODO, add the same columns as the normal dataframe
    dataDf = pd.DataFrame(columns=['ref', 'start_1-based', 'end_1-based', 'id', 'mapq', 'strand',
                                   'avg_coverage', 'fragment_count', 'fragment_count_per_base', 'FPKM',
                                   'TPM (based on fragment count)', 'TPM (based on average coverage)',
                                   'FPKM_no_rRNA', 'TPM (based on average coverage) no_rRNA',
                                   'TPM (based on fragment count) no_rRNA'])
    dataDf.to_csv(str(Path(args.output_path) / "{}.CDS_values.csv".format(args.sample_name)))
