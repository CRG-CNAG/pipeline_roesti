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

# from mwTools.general import glob_file_list
# from mwTools.bio import convert_genbank_to_annotation_df
# from mwTools.bio import convert_annotation_df_to_bed
# from mwTools.id import extract_fasta_id
# from mwTools.pandas import sort_df
import random
from glob import glob



def glob_file_list(fileList):
    "Apply glob on the list of file pattern, filter for files only, remove duplicated files."
    if fileList is None or fileList == '':
        fileList = None
    else:
        fileList = [fn for pattern in fileList for fn in glob(pattern) if Path(fn).is_file()]
        # Remove duplicated files from list (can happen if multiple patterns match same files)
        fileList = sorted(list(set(fileList)))
    return fileList


def convert_genbank_to_annotation_df(genomeBio, extractDNASeq=False, extractProteinSeq=False, verbose=1):
    
    species_name = genomeBio.annotations.get('organism')
    species_name = species_name if species_name is not None else ''
    
    def get_proteinGI(feature):
        proteinGI = feature.qualifiers.get('db_xref')
        proteinGI = int(re.sub(r'GI:', '', proteinGI[0])) if proteinGI is not None else None
        return proteinGI
    
    def get_attribute(feature, attName):
        if attName == 'proteinGI':
            att = get_proteinGI(feature)
        else:
            att = feature.qualifiers.get(attName)
            att = att[0] if att is not None else None
        return att
      
#     def get_attribute_dict(feature):
#         # Extract the first value of each qualifier
#         return {q:v[0] for q, v in feature.qualifiers.items()}

    def get_attribute_dict(feature):
        return {att:get_attribute(feature, att) for att in feature.qualifiers.keys()}
    
    hasWellDefinedGenomeSeq = \
        not (
            (genomeBio.seq.count('N') / len(genomeBio.seq)) > 0.5 or
            genomeBio.seq == ''
        )

    geneIdentifierPriorityList = ['locus_tag', 'gene', 'label', 'protein_id']
    
    CDSList = []
    if verbose >= 1: print("len(genomeBio.features):", len(genomeBio.features))
    for feature in genomeBio.features:
        attDict = get_attribute_dict(feature)
        location = convert_location_bio_to_dict(feature)
        locationBio = feature.location

        # Use the first of the attributes of the priority list found in the feature as main identifier
        featId = None
        for identifier in geneIdentifierPriorityList:
            if attDict.get(identifier) is not None:
                featId = attDict[identifier]
                break

        CDSDict = {'chromosome':genomeBio.id,
                   'id':featId,
                   'feature':feature.type,
                   'strand':location['strand'],
                   'start':location['start'],
                   'end':location['end'],
                   'locus_tag':attDict.get('locus_tag'),
                   'gene':attDict.get('gene'),
                   'protein_id':attDict.get('protein_id')
                   }
        if extractDNASeq and hasWellDefinedGenomeSeq:
            dnaSeqBio = SeqFeature(location=locationBio).extract(genomeBio)
            dnaSeq = str(dnaSeqBio.seq)
            CDSDict['DNA_seq'] = dnaSeq
        else:
            dnaSeqBio = None
            dnaSeq = None
        
        if extractProteinSeq and feature.type == 'CDS' and dnaSeqBio is not None:
            codonTableId = attDict.get('transl_table')
            codonStartPos = attDict.get('codon_start')    # in 1-based index
            codonStartPos = int(codonStartPos) if codonStartPos is not None else 1
            if codonTableId is not None:
                try:
                    proteinSeq = str(dnaSeqBio.seq[3*(codonStartPos - 1):].translate(table=codonTableId, cds=True))
                except TranslationError:
                    try:
                        proteinSeq = str(dnaSeqBio.seq[3*(codonStartPos - 1):].translate(table=codonTableId, cds=False))
                    except TranslationError:
                        proteinSeq = None

            CDSDict['protein_seq'] = proteinSeq

        CDSList.append(CDSDict)
        if verbose >= 2: print("\n\n")

    CDSDf = pd.DataFrame(CDSList)
    CDSDf = CDSDf.sort_values(by=['start', 'end', 'strand'])
    
    return CDSDf


def convert_annotation_df_to_bed(annotDf, chromosomeCol='chromosome', startCol='start', endCol='end',
                                 strandCol='strand', idCol='id', featureCol='feature', combineFeatureAndId=False,
                                 sort=True, sortBy=None, sortAscending=None
                                 ):
    """We assume that annotation features dataframe uses 0-based start-inclusive end-exclusive
    counting with start <= end."""

    df = annotDf
    if type(df) is pd.Series:
        df = df.to_frame().T
    sortColList = [startCol, endCol, strandCol]
    sortAscendingList = [True, True, True]
    if sort:
        if sortBy is not None:
            if 'chromosome' in sortBy:
                raise ValueError("chromosome is included by default and has to be sorted.")
            sortColList = sortBy + sortColList
            sortAscendingList = sortAscending + sortAscendingList
        df = df.sort_values(by=sortColList, ascending=sortAscendingList)
        df = sort_df(df, chromosomeCol, key=lambda x: (x.upper(), x[0].islower()), reverse=False)
    print(df)
    df.index.name = idCol

    if idCol not in df.columns:
        df = df.reset_index()

    bed_string = ""
    for index, annot in df.iterrows():
        # BED format uses 0-based start-inclusive end-exclusive counting (as in Python)
        # BED format start < end
        chromosomeName = annot[chromosomeCol]
        if combineFeatureAndId:
            bed_id = "{};{}".format(annot[featureCol], annot[idCol])
        else:
            bed_id = annot[idCol]
        bed_start = int(annot[startCol])
        bed_end = int(annot[endCol])
        bed_strand = annot[strandCol]
        if bed_start > bed_end:
            bed_start, bed_end = bed_end, bed_start
        bed_string += "{}\t{}\t{}\t{}\t{:d}\t{}\n".format(chromosomeName, str(bed_start), str(bed_end),
                                                          bed_id, 0, bed_strand)
    return bed_string


def extract_fasta_id(seqIdString):
    """
    Extract the id in the fasta description string, taking only the characters till the '|'.
    """
    seqId = None
    regex1 = re.compile(r'^\s*?([^|]+?)(\s|\||$)')
    seqIdSearch = re.search(regex1, seqIdString)
    if seqIdSearch:
        seqId = seqIdSearch.group(1)
    return seqId


def sort_df_by_value_list(df, col, valueList):
    sorterIndex = dict(zip(valueList, range(len(valueList))))
    ranInt = random.random_integers(1e6)
    rank = df[col].map(sorterIndex)
    rank.name = 'rank{:d}'.format(ranInt)
    return df.join(rank).sort_values(rank.name).drop(rank.name, axis=1)


def sort_df(df, col, key=None, reverse=False):
    """
    Take into account duplicated entries by using a rank column and applying the sort_df_by_value_list
    method.

    Example, sorting by alphabetical order grouping together capital and small case, with capital case first:
    sort_df(df, 'col', key=lambda x: (x.upper(), x[0].islower()))
    """

    def sorter(key=None, reverse=False):
        """
        Using python built-in sorted function.
        """
        def sorter_(series):
            series_list = list(series)
            return [series_list.index(i)
                    for i in sorted(series_list, key=key, reverse=reverse)]
        return sorter_

    df2 = df.reset_index(drop=True)
    sortedIndex = sorter(key=key, reverse=reverse)(df2[col])
    sortedValuesUnique = list(OrderedDict.fromkeys(np.array([df2[col].iloc[i] for i in sortedIndex])))
    return sort_df_by_value_list(df2, col=col, valueList=sortedValuesUnique)



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

        CDSDf = annotDf[annotDf['feature'] == 'CDS']
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
