#!/usr/bin/env python3
# coding: utf-8
import pandas as pd
import argparse
import re
from pathlib import Path
import textwrap
import subprocess
import shlex

import Bio
from Bio import SeqIO

from mwTools.general import glob_file_list
from mwTools.bio import convert_genbank_to_annotation_df
from mwTools.bio import convert_annotation_df_to_bed
from mwTools.id import extract_fasta_id
from mwTools.pandas import sort_df



class GenomeError(Exception):
    def __init__(self, message):
        self.message = message


def sort_bed_file(filePath, cwdPath):
    # Sort the BED file following same order than bedtools
    cmd = ("/users/lserrano/mweber/local/bin/sort {}".format(str(filePath)) +
           " -k 1,1 -k 2,2n -k 3,3n > {}.sorted".format(str(filePath)) +
           " && sleep 10 && mv {}.sorted {}".format(str(filePath), str(filePath)))
    cmd_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=str(cwdPath), shell=True)
    return cmd_output


def index_genome_files_bowtie2(genbankFileList=None, fastaFileList=None, outputName=None, outputDir='.'):

    outputPath = Path(outputDir).resolve()

    # Parse Genbank files

    annotDfList = []
    fastaFileToBeIndexedList = []
    genomeList = []

    if genbankFileList is not None:
        
        for genbankFilename in genbankFileList:

            genbankBioList = list(Bio.SeqIO.parse(genbankFilename, format="genbank"))
            print("\n### Parsing genbank file", genbankFilename)
            print("Nb of genome entries in Genbank file:", len(genbankBioList))

            for genomeBio in genbankBioList:
                
                print(genomeBio.id, genomeBio.description)
                
                # Convert genome DNA sequence to FASTA nucleotide.

                # Test if the genbank has a well-defined DNA sequence.
                hasWellDefinedSeq = True
                hasWellDefinedSeq = \
                    not (
                        (genomeBio.seq.count('N') / len(genomeBio.seq)) > 0.5 or
                        genomeBio.seq == ''
                    )

                if not hasWellDefinedSeq:
                    print("######\nWARNING: Genbank file", genbankFilename, "contains an undefined sequence.\n######")
                    # We still record this genome annotation but we will require later to have the genome sequence
                    # in fasta format below
                    genomeList.append({'genomeId':genomeBio.id, 'genomeSeq':None, 'genomeLength':0})
                else:
                    print("Converting DNA sequence from genbank to fasta.")
                    fastaId = genomeBio.id
                    fastaDescription = genomeBio.description
                    genomeSeq = textwrap.fill(str(genomeBio.seq), width=70)

                    genomeList.append({'genomeId':genomeBio.id, 'genomeSeq':str(genomeBio.seq),
                                       'genomeLength':len(genomeBio.seq)})

                    genomeFastaFile = str(outputPath / (fastaId + '.fna'))
                    fastaFileToBeIndexedList.append(genomeFastaFile)

                    with open(genomeFastaFile, 'w') as f:
                        # Note: we should not include any description because bowtie will need to capture only the id
                        # such that the id of the fasta file corresponds to the id of the BED files.
                        f.write(">{}\n{}".format(fastaId, genomeSeq))
                    print("Writing fasta file:", genomeFastaFile)


                # Convert genbank annotations to BED file
                annotDf = convert_genbank_to_annotation_df(genomeBio, verbose=1)
                annotDfList.append(annotDf)
                print("annotDf:\n", annotDf.head(20))


    # Parse fasta files

    if fastaFileList is not None:
        fastaFileToBeIndexedList += fastaFileList

        for fastaFile in fastaFileList:
            for fastaBio in SeqIO.parse(fastaFile, format="fasta"):
                fastaId = extract_fasta_id(fastaBio.id)

                # Check if we already have the same genome in Genbank format
                sameGenomeGenbank = [g for g in genomeList if g['genomeId'] == fastaId]
                if len(sameGenomeGenbank) == 0:
                    # Just add the new genome
                    genomeList.append({'genomeId':fastaId, 'genomeSeq':str(fastaBio.seq), 'genomeLength':len(fastaBio.seq)})
                elif len(sameGenomeGenbank) == 1:
                    # If the same genome had no defined sequence, add it
                    if sameGenomeGenbank['genomeSeq'] is None:
                        for i, g in enumerate(genomeList):
                            if g['genomeId'] == fastaId:
                                genomeList[i] = {'genomeId':fastaId, 'genomeSeq':str(fastaBio.seq), 'genomeLength':len(fastaBio.seq)}
                    else:
                        # If the genome already had a sequence, ERROR
                        print("ERROR: the DNA sequence for fasta file", fastaFile, "has the same id as the DNA sequence",
                              "already defined in the Genbank file.")
                        return
                elif len(sameGenomeGenbank) > 1:
                    print("ERROR: multiple DNA sequences between fasta and genbank files that have the same id.")


    # Check that we do not have any undefined DNA sequence
    undefinedGenomeList = [g for g in genomeList if g['genomeSeq'] is None]
    if len(undefinedGenomeList) > 0:
        message = ("ERROR: we have undefined DNA sequence(s):" + undefinedGenomeList.__repr__() +
                   "Please provide an additional fasta file containing the DNA sequence with the same id as the Genbank file,"
                   " or an updated Genbank file that contains the DNA sequence.")
        raise GenomeError(message)
        return

    # Choose output name

    if outputName is None:
        if len(genomeList) > 0:
            print("genomeList", [g['genomeId'] for g in genomeList])
            genomeIdComposed = "_".join([g['genomeId'] for g in genomeList])
            outputName = genomeIdComposed
        else:
            print("ERROR: genome id could not be extracted from the genbank or fasta files,",
                  "please use the -o option to set the output name.")
    outputPath = outputPath / outputName
    outputPath.mkdir(exist_ok=True)
    print("Output path:", str(outputPath))


    # Parse CDS and rRNA list

    if annotDfList != []:

        annotDf = pd.concat(annotDfList)
        # In fact, sorting here is useless because it does not follow the same sorting criteria as bedtool.
        annotDf.sort_values(['start'], inplace=True)
        annotDf = sort_df(annotDf, 'chromosome', key=lambda x: (x.upper(), x[0].islower()), reverse=False)

        CDSDf = annotDf[annotDf['feature'].str.contains('CDS')]
        filePath = outputPath / (outputName + '_CDS.bed')
        with filePath.open('w') as f:
            f.write(convert_annotation_df_to_bed(CDSDf))
        print("Writing CDS BED file", (outputName + '_CDS.bed'))
        # See Note 1 below
        sort_bed_file(filePath, outputPath)

        rRNADf = annotDf[annotDf['feature'] == 'rRNA']
        filePath = outputPath / (outputName + '_rRNA.bed')
        with filePath.open('w') as f:
            f.write(convert_annotation_df_to_bed(rRNADf))
        print("Writing rRNA BED file", (outputName + '_rRNA.bed'))
        # See Note 1 below
        sort_bed_file(filePath, outputPath)

        tRNADf = annotDf[annotDf['feature'] == 'tRNA']
        filePath = outputPath / (outputName + '_tRNA.bed')
        with filePath.open('w') as f:
            f.write(convert_annotation_df_to_bed(tRNADf))
        print("Writing tRNA BED file", (outputName + '_tRNA.bed'))
        # See Note 1 below
        sort_bed_file(filePath, outputPath)

        filePath = outputPath / (outputName + '_rRNA_tRNA.bed')
        with filePath.open('w') as f:
            f.write(convert_annotation_df_to_bed(pd.concat([tRNADf, rRNADf])))
        print("Writing rRNA_tRNA BED file", (outputName + '_rRNA_tRNA.bed'))
        # See Note 1 below
        sort_bed_file(filePath, outputPath)

    # Write genome list BED file
    # Note 1: we **need** to sort the chromosome names in alplhabetical order and then the
    # reads in start coordinates. Failing to do so will result in an error when computing the genome coverage using
    # bedtools.
    with (outputPath / (outputName + '_genome.bed')).open('w') as f:
        genomeList = sorted(genomeList, key=lambda x: (x['genomeId'].upper(), x['genomeId'][0].islower()))
        print("genomeList", [g['genomeId'] for g in genomeList])
        for g in genomeList:
            f.write("{}\t{}\n".format(g['genomeId'], g['genomeLength']))

    # Index reference genome

    print("We build the bowtie2 index from the following set of fasta files:\n", "\n".join(fastaFileToBeIndexedList))
    if len(fastaFileToBeIndexedList) == 0:
        print("ERROR: len(fastaFileToBeIndexedList) == 0")
        return
    # bowtie2Build = '/users/lserrano/mweber/Software/bowtie2-2.2.9/bowtie2-build'
    bowtie2Build = 'bowtie2-build'
    print("### Running bowtie2-build...")
    cmd = bowtie2Build + ' ' + ",".join(fastaFileToBeIndexedList) + " " + outputName
    print(cmd)
    cmd = shlex.split(cmd)
    cmd_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=str(outputPath))
    cmd_output = re.sub(r'\\n','\n',str(cmd_output))
    # print(cmd_output)
    print("### bowtie2-build finished")

    # Return the path to the basename of the indexed genome
    # This path is what bowtie2 aligner expects as input
    return str(outputPath / outputName)



if __name__ == "__main__":

    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    argParser = argparse.ArgumentParser(formatter_class=formatter_class)
    argParser.add_argument('-g', '--genbank', type=str, nargs='+',
                           help="List of filenames for reference DNA sequences in Genbank format. Wildcards such as ? * and [] can be used.")
    argParser.add_argument('-f', '--fasta', type=str, nargs='+',
                           help="List of filenames for reference DNA sequences in fasta format. Wildcards such as ? * and [] can be used.")
    argParser.add_argument('-o', '--output-name', type=str, dest='outputName',
                           help="Base name of the indexed genome and BED files. By default, the base name is set to the first id found in the Genbank or fasta files.")
    argParser.add_argument('-d', '--output-dir', type=str, dest='outputDir', default='.',
                           help="Directory of the indexed genome and BED files.")
    args = argParser.parse_args()


    # Example
    # args = argParser.parse_args('-g /users/lserrano/mweber/Databases/RefSeq/Release_2016_10/Genomes_extracted/GCF_000183385.1_ASM18338v1_genomic.gbff /users/lserrano/mweber/Databases/RefSeq/Release_2016_10/Genomes_extracted/GCF_000092585.1_ASM9258v1_genomic.gbff -f /users/lserrano/mweber/Translation_model/Ribosome_profiling_data/bowtie2_indexed_genome/NC_013948.1/example.fasta -d /users/lserrano/mweber/Translation_model/Ribosome_profiling_data/bowtie2_indexed_genome'.split())

    genbankFileList = glob_file_list(args.genbank)
    fastaFileList = glob_file_list(args.fasta)

    genbankFileListPrint = genbankFileList if genbankFileList is not None else [""]
    fastaFileListPrint = fastaFileList if fastaFileList is not None else [""]
    print("genbankFileList:\n", "\n".join(genbankFileListPrint))
    print("fastaFiles:\n", "\n".join(fastaFileListPrint))
    
    indexedGenomePath = index_genome_files_bowtie2(genbankFileList=genbankFileList, fastaFileList=fastaFileList,
                                                   outputName=args.outputName, outputDir=args.outputDir)
    print("### Indexed genome path:", indexedGenomePath)
