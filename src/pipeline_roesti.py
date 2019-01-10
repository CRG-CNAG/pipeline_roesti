#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("darkgrid")
seaborn.set_style("whitegrid")

from ruffus import *
import ruffus.cmdline as cmdline
import time
from pathlib import Path
import os
import shutil
import re
import argparse
from glob import glob
from itertools import islice
import datetime
import subprocess
import shlex
import pandas as pd
import json
from socketIO_client import SocketIO, LoggingNamespace
from send_socket_message import send_socket_message

from index_genome_files_bowtie2 import index_genome_files_bowtie2
from mwTools.general import glob_file_list
from mwTools.general import open_by_suffix

# --------------------------------------------------

pipeline_name = 'pipeline_roesti'
pipelineDoc = pipeline_name + """

Pipeline to analyze RNA-seq data from RNA-seq and Ribo-seq (ribosome profiling) experiments.
Python 3.5 script. |

Author: Marc Weber |
Last updated: 2017.11.23 |
Affiliation: Center for Genomic Regulation, Luis Serrano's lab |
email: webermarcolivier@gmail.com; marc.weber@crg.eu |

Remark: we refer to insert as the original RNA fragment flanked by adapters and read by sequencing (as opposed to "insert" as the inner interval between two paired reads).

"""

# Pipeline test options
run_on_cluster = True

# Command-line arguments
formatter_class = argparse.RawDescriptionHelpFormatter
formatter_class = argparse.ArgumentDefaultsHelpFormatter
parser = cmdline.get_argparse(description=pipelineDoc, formatter_class=formatter_class)
parser.add_argument('-g', '--ref-genbank', dest='refGenbank', type=str, nargs='+',
                    help="List of filenames for reference DNA sequences in Genbank format. Wildcards such as ? * and [] can be used. Annotations in genbank format will be used to create the CDS and rRNA BED files.")
parser.add_argument('-f', '--ref-fasta', dest='refFasta', type=str, nargs='+',
                    help="List of filenames for reference DNA sequences in fasta format. Wildcards such as ? * and [] can be used.")
parser.add_argument('--ref-output-name', dest='refOutputName', type=str,
                    help="Base name of the indexed genome and BED files. By default, the base name is set to the concatenated id's found in the Genbank or fasta files.")
parser.add_argument('--ref-output-dir', dest='refOutputDir', type=str, default='.',
                    help="Directory of the indexed genome and BED files.")
parser.add_argument('-i', '--fastq-files', dest='input_fastq_files', default=['*.fastq*'], nargs='+',
                    help="List of input fastq files separated by white spaces. Wildcards such as ? * and [] can be used.")
parser.add_argument('-o', '--output-dir', dest='output_dir', default='.',
                    help='Output directory.')
parser.add_argument('--test', dest='run_test', action='store_true',
                    help='Run a test by taking as input the first 25k reads of the first fasta file.')
parser.add_argument('--run-locally', dest='run_locally', action='store_true',
                    help='Run the pipeline on local machine (as opposed to submitting jobs to the cluster). Note that multiple threads could be used.')
parser.add_argument('--pipeline-name', dest='pipeline_name', default=pipeline_name,
                    help='Name of the pipeline. Important: the history of up-to-date files is kept in a database of the same name in the output folder, .pipeline_roesti.ruffus_history.sqlite.')
parser.add_argument('--library-type', dest='library_type', default='rna-seq', choices=['ribo-seq', 'rna-seq'],
                    help="Type of RNA-seq library. In ribosome profiling data analysis, additional fragment filtering is applied in order to select ribosome footprints.")
parser.add_argument('--seq-end', dest='seq_end', default='paired-end', choices=['single-end', 'paired-end'],
                    help="Single-end or paired-end sequencing data.")
parser.add_argument('--adapter-seq-fw', dest='trim_adapter_seq_forward', default='',
                    help='Adapter sequence forward. Default will be set to TruSeq Universal Adapter 5’ AGATCGGAAGAGCACACGTCT')
parser.add_argument('--adapter-seq-rv', dest='trim_adapter_seq_reverse', default='',
                    help='Adapter sequence reverse. Default will be set to TruSeq Universal Adapter 5’ AGATCGGAAGAGCACACGTCT for rna-seq library, and to Illumina RNA PCR Primer GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTA for ribo-seq library.')
parser.add_argument('--no-trimming', dest='no_trimming', action='store_true',
                    help='Do not trim reads before alignment. By default bad quality bases are trimmmed from read ends.')
parser.add_argument('--remove-rRNA', dest='remove_rRNA', action='store_true',
                    help='Remove reads that align to the rRNA.')
parser.add_argument('--align-quality-threshold', dest='filter_alignments_quality_threshold', default=15, type=int,
                    help='Filter out reads with alignment quality score MAPQ smaller than the threshold.')
parser.add_argument('--indexed-ref-genome', dest='align_indexed_ref_genome_path', default='/users/lserrano/mweber/RNA-seq_data/bowtie2_indexed_genome/Mpn/NC_000912',
                    help="Path to the basename of the index for the reference genome built with bowtie2-build.")
parser.add_argument('--rRNA-bedfile', dest='rRNA_bedfile',
                    default="/users/lserrano/mweber/Research_cloud/Mycoplasma_pneumoniae_experimental_data/Annotation/mpn_rRNA.bed",
                    help="Path to the BED file of rRNA regions. Reads aligning in the first rRNA region will be used to determine the strandness.")
parser.add_argument('--rRNA-tRNA-bedfile', dest='rRNA_tRNA_bedfile',
                    default="/users/lserrano/mweber/Research_cloud/Mycoplasma_pneumoniae_experimental_data/Annotation/mpn_rRNA_tRNA.bed",
                    help="Path to the BED file of rRNAs and tRNAs regions of the genome. If the option remove-rRNA is set, all reads aligning in these regions will be filtered out.")
parser.add_argument('--genome-bedfile', dest='genomeBedFile',
                    default="/users/lserrano/mweber/Research_cloud/Mycoplasma_pneumoniae_experimental_data/Genome/NC_000912.1.genome",
                    help="Path to the BED file genome. Simple BED file that lists the names of the chromosomes (or scaffolds, etc.) and their size (in basepairs).")
parser.add_argument('--genome-CDS-bedfile', dest='genomeCDSBedFile',
                    default="/users/lserrano/mweber/Research_cloud/Mycoplasma_pneumoniae_experimental_data/Annotation/mpn_CDS.bed",
                    help="Path to the BED file for all CDS. Will be used to count mRNA fragments for each gene. If set to empty string \"\", the computation of fragment count per CDS will be skipped.")
parser.add_argument('--nthreads', dest='nThreads', default=12, type=int,
                    help='Number of threads to use in each cluster node (shared memory). This will reduce computational time, in particular for bowtie2 (alignment).')
parser.add_argument('--njobsmax', dest='njobs', default=200, type=int,
                    help='Number of concurrent jobs to launch on the cluster (each job will use nThreads nodes in shared memory).')
parser.add_argument('--bash-profile', dest='bash_profile', default='',
                    help='Bash profile is executed before each job on the cluster in order to load the dependencies. By default bash profile path is automatically detected in user\'s home directory, this option sets the path manually.')
parser.add_argument('--analysisId', default=None, type=str,
                    help="String that identifies the overall pipeline run. It is independent from the jobids of the job submissions on the cluster grid engine.")
parser.add_argument('--deleteIntermediateFiles', dest='delete_intermediate_files', action='store_true', default=False)
parser.add_argument('--sendMessageToWebServer', action='store_true',
                    help="Send a message to the webserver dbspipe when the pipeline has finished. Only for pipeline launched by the web server application.")
parser.add_argument('--host', dest='host', default='cluster', type=str)
options = parser.parse_args()


#################################
# TODO

# change the name of the tRNA_rRNA file to something more consistent like
#     the file "regions_to_remove"

# INCLUDE RRNA COUNTS IN THE FINAL VALUES


## Global

pipeline_name = options.pipeline_name
run_locally = options.run_locally

# Glob fastq files following list of patterns
fastqFiles = [fn for pattern in options.input_fastq_files for fn in glob(pattern) if Path(fn).is_file()]
if len(fastqFiles) == 0:
    raise ValueError("ERROR: no input file exists with path(s) {}".format(options.input_fastq_files))
    raise SystemExit
# Remove duplicated files from list (can happen if multiple patterns match same files)
fastqFiles = sorted(list(set(fastqFiles)))
# Filter out files that start with "test"
# Note: This is not strictly necessary. Maybe we want to test the pipeline with already existing
# test fastq files.
# fastqFiles = [fn for fn in fastqFiles if not re.match(r'test.*', Path(fn).name)]
fastqPath = Path(fastqFiles[0]).parent
print("fastqFiles:", fastqFiles)


def group_paired_end_fastq_files(fastqFiles):
    """Group paired-end fastq files together"""
    fastqFilesPE = []
    for f1 in fastqFiles:
        f1match = re.match(r"^(.+)(read1|r1)\.fastq(\.gz)?$", Path(f1).name)
        if f1match:
            read1Filename = f1
            read1FastaBasename = f1match.group(1)
            for f2 in fastqFiles:
                if re.match(r'^' + read1FastaBasename + r"(read2|r2)\.fastq(\.gz)?$", Path(f2).name):
                    read2Filename = f2
                    fastqFilesPE.append((read1Filename, read2Filename))
    return fastqFilesPE


# Paths
rootPath = Path('/users/lserrano/mweber')
scriptPath = rootPath / 'Research_cloud' / 'RNA-seq_data_analysis' / 'src'
# rnaSeq_data_path = rootPath / 'Translation_model' / 'Ribosome_profiling_data'
# Note: here we must use the full pathname!!!
outputPath = Path(options.output_dir).resolve()
os.chdir(str(outputPath))
# Make pipeline folder
pipeline_path = outputPath
pipeline_path.mkdir(exist_ok=True)
pipeline_folder = pipeline_path.name


# The isis filesystem on the cluster can be very slow to propagate file changes, which can be a problem in the pipeline
# because we frequently write and read files. We introduce a waiting loop in between file creation and file opening in
# order to let the filesystem update the file changes. The loop will keep checking for the existence of the file every X
# seconds.
sleepTimeFilesystem = 2    # seconds


# Cluster queues
short_queue = 'short-sl7'
long_queue = 'long-sl7'

# The bash profile together with the load_dependencies script will be executed before each command in the pipeline.
# This is necessary when running computation on the nodes of the cluster.
if options.bash_profile == '':
    if (Path.home() / '.bash_profile').exists():
        options.bash_profile = Path.home() / '.bash_profile'
    elif (Path.home() / '.bashrc').exists():
        options.bash_profile = Path.home() / '.bashrc'
    else:
        print("ERROR: bash profile could not be found.")
        raise SystemError
loadDependenciesScriptPath = scriptPath / "load_dependencies.sh"
if options.host == 'cluster':
    cmd_source_bash_profile = ". {} && . {} && cd {} &&".format(options.bash_profile,
                                                                str(loadDependenciesScriptPath),
                                                                str(outputPath))
else:
    cmd_source_bash_profile = "cd {} &&".format(str(outputPath))

# Global bowtie2 options
options.phredEncoding = 'phred33'
# If running locally, only run 1 job
if run_locally:
    options.njobs = 1

## trim adapter
options.no_trimming
options.trim_adapter_min_length = 12
options.trim_adapter_trim_end_quality = 10
options.trim_adapter_nthreads = min(4, options.nThreads)
if options.library_type == 'ribo-seq':
    if options.trim_adapter_seq_forward == '':
        options.trim_adapter_seq_forward = 'AGATCGGAAGAGCACACGTCT'
    if options.trim_adapter_seq_reverse == '':
        options.trim_adapter_seq_reverse = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTA'
elif options.library_type == 'rna-seq':
    if options.trim_adapter_seq_forward == '':
        options.trim_adapter_seq_forward = 'AGATCGGAAGAGCACACGTCT'
    if options.trim_adapter_seq_reverse == '':
        options.trim_adapter_seq_reverse = 'AGATCGGAAGAGCACACGTCT'

## align
options.align_alignmentMode = 'end-to-end'

# As a reference, the default options for bowtie2 in the mode --very-sensitive are:
# -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
if options.library_type == 'ribo-seq':
    options.align_nMismatches = 1
    options.align_seedLength = 12
elif options.library_type == 'rna-seq':
    options.align_nMismatches = 0
    options.align_seedLength = 20
options.align_seedInterval = 'S,1,1.15'    # mode sensitive. For read length of 50, seed interval is 1 + 1.15*sqrt(50) = 9.13
options.align_seedInterval = 'S,1,0.50'    # mode very sensitive. For read length of 50, seed interval is 1 + 0.5*sqrt(50) = 4.53

options.align_maxAlignAttempts = 20    # very sensitive: -D 20
options.align_maxReSeed = 3    # very sensitive: -R 3

if options.library_type == 'ribo-seq':
    options.align_maxInsertLength = 400
elif options.library_type == 'rna-seq':
    options.align_maxInsertLength = 1000     # This is the theoretical maximum fragment length in the library preparation
options.align_max_reported_alignments = 0    # set to 0 to deactivate the option

## convert_sam_to_bam
if run_locally:
    options.samtools_sort_tmp_dir = '.'
else:
    options.samtools_sort_tmp_dir = '$TMPDIR'
options.samtools_sort_max_mem = 24000    # M
options.samtools_sort_nthread = 4
# Note: maximum memory per thread **has to be an integer**, otherwise it is interpreted as bytes
# Note: we have to allocate some free memory for the main samtools thread, probably for the
# file merging, otherwise the job will reach memory limit.
options.samtools_sort_max_mem_per_thread = int((options.samtools_sort_max_mem - 6000) / options.samtools_sort_nthread)

## filter_alignments
# ...

## Write information about the pipeline in a text file
# pipelineDocFile = pipeline_path / (pipeline_name + '.{:d}.readme'.format(options.jobid))
pipelineDocFile = pipeline_path / (pipeline_name + '.readme')

# Standard python logger which can be synchronised across concurrent Ruffus tasks
# options.log_file = str(pipeline_path / (pipeline_name + '.{:d}.log'.format(options.jobid)))
options.log_file = str(pipeline_path / (pipeline_name + '.log'))
options.verbose = 2
logger, logger_mutex = cmdline.setup_logging(pipeline_name, options.log_file, options.verbose)
with logger_mutex:
    logger.debug("Pipeline options:\n" + str(options))
pipelineDoc += "Type of RNA-seq data:" + options.library_type + "\n\n"
print("\n\n########## TYPE OF RNA-SEQ DATA: ",options.library_type, "\n\n")
pipelineDoc += "Pipeline options:\n\n"
for key, value in vars(options).items():
    pipelineDoc += '{}:{}\n'.format(key, value)


# Generate test files
writeTestFiles = options.run_test
# takes 25k first reads in the first pair of fastq files
if writeTestFiles:
    nReads = int(50e3)
    if options.seq_end == 'single-end':
        exampleFiles = [fastqFiles[0]]
    elif options.seq_end == 'paired-end':
        fastqFilesPE = group_paired_end_fastq_files(fastqFiles)
        print("fastqFilesPE", fastqFilesPE)
        if len(fastqFilesPE) == 0:
            raise ValueError("No fastq file pair has been found within the list of input files. Check that the filenames are identical for read1 and read2.")
        exampleFiles = fastqFilesPE[0]
    
    for i, fastqFilename in enumerate(exampleFiles):
        print("Writing test files from head of file: ",fastqFilename)
        with open_by_suffix(str(fastqFilename)) as fastqFile:
            head = list(islice(fastqFile, 4*nReads))  # Note: each read is 4 lines in fastq format
        if options.seq_end == 'single-end':
            testFilename = str(outputPath / "test.fastq")
        elif options.seq_end == 'paired-end':
            testFilename = str(outputPath / ("test_read{:d}.fastq".format(i + 1)))
        with open(testFilename, "w") as testFile:
            testFile.writelines(head)
        fastqPath = outputPath

if options.run_test:
    # Test fastq files with 25k first reads
    print("Taking as a test the 25k first reads of the first (pair of) fastq file(s).")
    fastqFiles = [str(filepath.resolve()) for filepath in fastqPath.glob('test*.fastq') if filepath.is_file()]

print("fastqFiles: ", fastqFiles)


#############################################################################


# Tools
def estimateNbLines(filename, learn_size_bytes=1024*1024):
    """ Estimate the number of lines in the given file without reading the whole file."""

    file_size = os.path.getsize(filename)
    learn_size_bytes = min(learn_size_bytes, file_size)

    if Path(filename).suffix == '.gz':
        # Rough approximation for gzipped fasta file of RNA-seq reads,
        # ~40k reads per 1MB of gz fasta file, each read is 4 lines
        numLines = (file_size/(1024*1024)) * 40e3 * 4
    else:
        with open(filename, 'rb') as file:
            buf = file.read(learn_size_bytes)
            numLines = file_size / (len(buf) // buf.count(b'\n'))

    return numLines


def printTimeDelta(time_delta):
    s = time_delta.seconds + time_delta.days*24*3600
    hours, remainder = divmod(s, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:d}:{:d}:{:f}".format(hours, minutes, seconds)


def plot_trimmed_reads_length_dist(summaryFilePath):
    histData = []
    with open(summaryFilePath,'r') as summaryFile:
        for line in summaryFile:
            if line.strip() == 'Read length distribution after trimming:' \
                    or line.strip() == 'Length distribution of reads after trimming:':
                break
        for line in summaryFile:
            regexSearch = re.search(r'([0-9]+):?\s*([0-9]+)', line.strip())
            if regexSearch:
                length = int(regexSearch.group(1))
                counts = int(regexSearch.group(2))
                histData.append((length, counts))
    
    x, y = zip(*histData)
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel('length trimmed read (insert)')
    fig.savefig(summaryFilePath + '.length_dist.png', dpi=200)
    
    return histData


def wait_for_file(filename, sleepTime=2):
    filePath = Path(filename)
    while not filePath.exists():
        time.sleep(sleepTime)
    return


def wait_for_any_of_files(filenameList, sleepTime=2):
    filePathList = [Path(fn) for fn in filenameList]
    while not any([fp.exists() for fp in filePathList]):
        print("sleeping...")
        time.sleep(sleepTime)
    return




#############################################################################

# Start DRMAA session
if run_on_cluster:
    # Start shared drmaa session for all jobs / tasks in pipeline
    # drmaa Open Grid Forum API
    import drmaa
    drmaa_session = drmaa.Session()
    drmaa_session.initialize()
    from ruffus.drmaa_wrapper import run_job, error_drmaa_job

iTask = 0

#############################################################################

cmdSendSocketMessageJobStarted = '{} --analysisId {} --status 1 {}'.format(str(scriptPath / 'send_socket_message.py'),
                                                                           options.analysisId,
                                                                           '--sendMessageToWebServer' if options.sendMessageToWebServer else '')

#############################################################################

infoStr = "\n\n#############################\n"
infoStr += """
Task: index_genome_files
Indexing reference DNA sequences for bowtie2 aligner.

Software: bowtie2-build version 2.2.9

Parameters used: all default\n\n"""
pipelineDoc += infoStr


computeIndexedGenome = (options.refGenbank is not None) or (options.refFasta is not None)
task_name = "index_genome_files"

if computeIndexedGenome:

    genbankFileList = glob_file_list(options.refGenbank)
    fastaFileList = glob_file_list(options.refFasta)
    
    with logger_mutex:
        genbankFileListPrint = genbankFileList if genbankFileList is not None else [""]
        fastaFileListPrint = fastaFileList if fastaFileList is not None else [""]
        logger.debug("genbankFileList:\n" + "\n".join(genbankFileListPrint))
        logger.debug("fastaFiles:\n" + "\n".join(fastaFileListPrint))

    refOutputDir = options.refOutputDir
    
    indexedGenomePath = index_genome_files_bowtie2(genbankFileList=genbankFileList, fastaFileList=fastaFileList,
                                                   outputName=options.refOutputName, outputDir=refOutputDir)
    with logger_mutex:
        logger.debug("### Indexed genome path:\n" + indexedGenomePath)


    options.align_indexed_ref_genome_path = indexedGenomePath
    options.rRNA_bedfile = indexedGenomePath + '_rRNA.bed'
    options.rRNA_tRNA_bedfile = indexedGenomePath + '_rRNA_tRNA.bed'
    options.genomeBedFile = indexedGenomePath + '_genome.bed'
    options.genomeCDSBedFile = indexedGenomePath + '_CDS.bed'


#############################################################################


infoStr = "\n\n#############################\n"
infoStr += """
Task: trim_adapter_PE_reads.
Trimming adapters from PE reads and extracting insert.

Software: SeqPurge 0.1-478-g3c8651b

Parameters used:\n\n"""
infoStr += "-min_len <int>      Minimum read length after adapter trimming. Shorter reads are discarded. Default value: '15'\n"
infoStr += "min_len " + str(options.trim_adapter_min_length) + "\n\n"
infoStr += "-a1 <string>        Forward adapter sequence (at least 15 bases). Default value: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'\n"
infoStr += "a1 " + options.trim_adapter_seq_forward + "\n\n"
infoStr += "-a2 <string>        Reverse adapter sequence (at least 15 bases). Default value: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'\n"
infoStr += "a2 " + options.trim_adapter_seq_reverse + "\n\n"
infoStr += "-threads " + str(options.nThreads) + "\n\n"
pipelineDoc += infoStr

taskPathList = []
iTask += 1
task_name = "trim_adapter_PE_reads"
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
taskPathList.append(task_path)
@follows(mkdir(str(task_path)))
# Only run this analysis for paired-end library type
@active_if(options.seq_end == 'paired-end')
# Group together file pairs
@collate(
    fastqFiles,

    formatter(r'^/?(.+/)*(?P<SAMPLENAME>.+)(read[12]|r[12])\.fastq(\.gz)?$'),

    # Create output parameter supplied to next task
    [str(task_path) + "/{SAMPLENAME[0]}read1.trimmed.fastq.gz",   # paired file 1
     str(task_path) + "/{SAMPLENAME[0]}read2.trimmed.fastq.gz",   # paired file 2
     str(task_path) + "/{SAMPLENAME[0]}.trimmed.nreads"],

    # Extra parameters for our own convenience and use
    [str(task_path) + "/{SAMPLENAME[0]}unpaired.read1.fastq.gz",    # unpaired file 1
     str(task_path) + "/{SAMPLENAME[0]}unpaired.read2.fastq.gz"],   # unpaired file 2

    task_name, logger, logger_mutex)
def trim_adapter_PE_reads(input_files,
                          output_paired_files,
                          discarded_unpaired_files,
                          task_name, logger, logger_mutex):
    if len(input_files) != 2:
        raise Exception("One of read pairs %s missing" % (input_files,))

    if options.no_trimming:
        with logger_mutex:
            logger.debug("Input files:")
            logger.debug(input_files[0])
            logger.debug(input_files[1])

        # Create symbolic links to the original read files
        Path(output_paired_files[0]).symlink_to(Path(input_files[0]).resolve())
        Path(output_paired_files[1]).symlink_to(Path(input_files[1]).resolve())

        # Compute nb of valid reads and write in file
        cmd = "wc -l < " + output_paired_files[0]
        cmd_output = subprocess.check_output(cmd, shell=True)
        cmd_output = cmd_output.decode().strip()
        print(cmd_output)
        nValidReads = int(cmd_output)

        metadata_filename = output_paired_files[2]
        with open(metadata_filename,'w') as metadata_file:
            metadata_file.write(str(nValidReads) + '\n')

    else:

        nReadsApprox = max([estimateNbLines(input_file, 100*1024*1024) / 4. for input_file in input_files])
        # Approximate computation time in seconds per reads per thread measured for test run
        comp_time = (nReadsApprox*(((2*60+59)*16)/12.5e6))/options.trim_adapter_nthreads
        with logger_mutex:
            logger.debug("trim_adapter_PE_reads, nReadsApprox " + str(nReadsApprox))
            logger.debug("trim_adapter_PE_reads, comp_time " + str(comp_time))

        summaryFilename = re.search(r'(.+)\.fastq.*', output_paired_files[0]).group(1) + ".summary"

        # We send the message to the web server from the first node of computation when the job has started.
        cmd = (cmd_source_bash_profile +
               cmdSendSocketMessageJobStarted + " && " +
               " SeqPurge " +
               " -in1 {} -in2 {} ".format(input_files[0], input_files[1]) +
               " -out1 {} -out2 {} ".format(output_paired_files[0], output_paired_files[1]) +
               " -a1 " + options.trim_adapter_seq_forward + " -a2 " + options.trim_adapter_seq_reverse + " " +
               " -min_len " + str(options.trim_adapter_min_length) +
               " -threads " + str(options.trim_adapter_nthreads) +
               " -summary " + summaryFilename)
        with logger_mutex:
            logger.debug(cmd)
            logger.debug("estimated computation time " + str(datetime.timedelta(seconds=comp_time)))

        try:
            stdout_res, stderr_res = "",""
            walltime = datetime.timedelta(seconds=max(4*comp_time, 1*60*60))
            if walltime < datetime.timedelta(hours=6):
                job_queue_name = short_queue
            else:
                job_queue_name = long_queue
            job_other_options = " -pe smp " + str(options.trim_adapter_nthreads) +\
                                " -q " + job_queue_name +\
                                " -l h_rt=" + printTimeDelta(walltime) +\
                                " -l h_vmem=3G,virtual_free=3G" +\
                                " -cwd"

            # ruffus.drmaa_wrapper.run_job
            stdout_res, stderr_res  = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                              drmaa_session=drmaa_session, run_locally=run_locally,
                                              job_other_options=job_other_options, retain_job_scripts=False)

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            if options.sendMessageToWebServer:
                send_socket_message(options.analysisId, -1)
            raise Exception("\n".join(map(str,["Failed to run:", cmd, err, stdout_res, stderr_res])))

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

        wait_for_file(summaryFilename, sleepTimeFilesystem)

        # Extract nb of valid reads from the summary file and write to a small metadata file
        with open(summaryFilename,'r') as summary_file:
            nValidReads = None
            nRemovedReads = None
            nRawReads = None
            for line in summary_file:
                regex = re.search(r'Reads \(forward \+ reverse\): ([0-9]+).*', line)
                if regex:
                    nRawReads = int(regex.group(1))
                    print('nRawReads', nRawReads)
                regex = re.search(r'Removed reads: ([0-9]+) of .*', line)
                if regex:
                    nRemovedReads = int(regex.group(1))
                    print('nRemovedReads', nRemovedReads)
            if nRawReads is not None and nRemovedReads is not None:
                nValidReads = nRawReads - nRemovedReads

        metadata_filename = output_paired_files[2]
        with open(metadata_filename,'w') as metadata_file:
            metadata_file.write(str(nValidReads) + '\n')

        # Small analysis that will be run on the login node
        plot_trimmed_reads_length_dist(summaryFilename)

    with logger_mutex:
        logger.debug("Trimming of pair-end reads finished.")



#############################################################################


infoStr = "\n\n#############################\n"
infoStr += """
Task: trim_adapter_SE_reads
Trimming adapters from SE reads and extracting insert.

Software: skewer-0.2.2

Parameters used:\n\n"""
infoStr += "-x <str> Adapter sequence/file (AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)\n"
infoStr += "-x " + options.trim_adapter_seq_forward + "\n\n"

infoStr += "-l, --min <int> The minimum read length allowed after trimming; (18)\n"
infoStr += "-l " + str(options.trim_adapter_min_length) + "\n\n"

infoStr += "-q, --end-quality  <int> Trim 3' end until specified or higher quality reached; (0)\n"
infoStr += "-q " + str(options.trim_adapter_trim_end_quality) + "\n\n"

infoStr += "-t, --threads <int>   Number of concurrent threads [1, 32]; (1)\n"
infoStr += "-t " + str(options.trim_adapter_nthreads) + "\n\n"

pipelineDoc += infoStr


iTask += 1
task_name = "trim_adapter_SE_reads"
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
taskPathList.append(task_path)
@follows(mkdir(str(task_path)))
# Only run this analysis for single-end library type
@active_if(options.seq_end == 'single-end')
@transform(fastqFiles,

    formatter(r'^/?(.+/)*(?P<SAMPLENAME>.+)\.fastq(\.gz)?$'),

    # Create output parameter supplied to next task
    [str(task_path) + "/{SAMPLENAME[0]}.trimmed.fastq",
     str(task_path) + "/{SAMPLENAME[0]}.trimmed.nreads"],

    task_name, task_path, logger, logger_mutex)
def trim_adapter_SE_reads(input_file,
                          output_files,
                          task_name, task_path, logger, logger_mutex):

    if type(input_file) is not str:
        raise Exception("Input file should be unique, input_file: %s" % (input_file,))

    if options.no_trimming:
        with logger_mutex:
            logger.debug("Input files:")
            logger.debug(input_file)

        # Create symbolic links to the original read files
        Path(output_files[0]).symlink_to(Path(input_file).resolve())

        # Compute nb of valid reads and write in file
        cmd = "wc -l < " + output_files[0]
        cmd_output = subprocess.check_output(cmd, shell=True)
        cmd_output = cmd_output.decode().strip()
        print(cmd_output)
        nValidReads = int(cmd_output)

        metadata_filename = output_files[1]
        with open(metadata_filename,'w') as metadata_file:
            metadata_file.write(str(nValidReads) + '\n')

    else:


        nReadsApprox = estimateNbLines(input_file, 100*1024*1024) / 4.
        # Approximate computation time in seconds per reads per thread measured for test run
        comp_time = (nReadsApprox * (((2*60) * 1) / 5e4)) / options.trim_adapter_nthreads

        # We send the message to the web server from the first node of computation when the job has started.
        cmd = (cmd_source_bash_profile +
               cmdSendSocketMessageJobStarted + " && " +
               " skewer-0.2.2-linux-x86_64 " +
               " " + input_file +
               " -x " + options.trim_adapter_seq_forward +
               " -l " + str(options.trim_adapter_min_length) +
               " -q " + str(options.trim_adapter_trim_end_quality) +
               " -t " + str(options.trim_adapter_nthreads) +
               " --quiet" +
               " -1 > " + output_files[0])

        with logger_mutex:
            logger.debug(cmd)
            logger.debug("estimated computation time " + str(datetime.timedelta(seconds=comp_time)))

        try:
            stdout_res, stderr_res = "",""
            walltime = datetime.timedelta(seconds=max(4*comp_time, 1*60*60))
            if walltime < datetime.timedelta(hours=6):
                job_queue_name = short_queue
            else:
                job_queue_name = long_queue
            job_other_options = " -pe smp " + str(options.trim_adapter_nthreads) +\
                                " -q " + job_queue_name +\
                                " -l h_rt=" + printTimeDelta(walltime) +\
                                " -l h_vmem=3G,virtual_free=3G" +\
                                " -cwd"

            # ruffus.drmaa_wrapper.run_job
            stdout_res, stderr_res  = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                              drmaa_session=drmaa_session, run_locally=run_locally,
                                              job_other_options=job_other_options, retain_job_scripts=False)

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            if options.sendMessageToWebServer:
                send_socket_message(options.analysisId, -1)
            raise Exception("\n".join(map(str,["Failed to run:", cmd, err, stdout_res, stderr_res])))

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

        trimmedFilename = output_files[0]

        # Compute nb of valid reads and write in file
        cmd = "wc -l < " + trimmedFilename
        cmd_output = subprocess.check_output(cmd, shell=True)
        cmd_output = cmd_output.decode().strip()
        print(cmd_output)
        nValidReads = int(cmd_output)

        metadata_filename = output_files[1]
        with open(metadata_filename,'w') as metadataFile:
            metadataFile.write(str(nValidReads) + '\n')

    with logger_mutex:
        logger.debug(cmd_output)
        logger.debug("Trimming of single-end reads finished.")



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """Align reads to genome.

Aligning reads to the reference genome.

We align all the trimmed raw reads to the genome, including rRNA and tRNAs, which we will filter in the next stage. Note that we get tons of very short reads (<20 bp) that will virtually map everywhere on the genome but with a very low mapping quality (MAPQ). These will give raise to paired aligned reads that align discordantly, most often because they get assigned random position on the genome and does not satisfy the excpected insert size or relative orientation (Forward-Reverse). We will only analyze the reads that align concordantly in the next stage.

In the case of ribosome profiling, footprints are very short inserts. Sensitive alignment tool might be better to align those. We use a small seed length and a smaller seed step, in order to improve alignment of shorter reads and allow more mismatches.

Software: bowtie2-align-s version 2.2.9

Parameters used: 

"""

infoStr += "options.align_indexed_ref_genome_path: " + str(options.align_indexed_ref_genome_path) + "\n\n"

infoStr += "--phred33          Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the \"Phred+33\" encoding, which is used by the very latest Illumina pipelines.\n"
infoStr += "Phred quality encoding: " + options.phredEncoding + "\n\n"

infoStr += "--end-to-end       In this mode, Bowtie 2 requires that the entire read align from one end to the other, without any trimming (or \"soft clipping\") of characters from either end. The match bonus --ma always equals 0 in this mode, so all alignment scores are less than or equal to 0, and the greatest possible alignment score is 0. This is mutually exclusive with --local. --end-to-end is the default mode.\n"
infoStr += "--local            In this mode, Bowtie 2 does not require that the entire read align from one end to the other. Rather, some characters may be omitted (\"soft clipped\") from the ends in order to achieve the greatest possible alignment score. The match bonus --ma is used in this mode, and the best possible alignment score is equal to the match bonus (--ma) times the length of the read. Specifying --local and one of the presets (e.g. --local --very-fast) is equivalent to specifying the local version of the preset (--very-fast-local). This is mutually exclusive with --end-to-end. --end-to-end is the default mode.\n"
infoStr += "options.align_alignmentMode: " + options.align_alignmentMode + "\n\n"

infoStr += """-N <int>
Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.\n"""
infoStr += "-N " + str(options.align_nMismatches) + "\n\n"

infoStr += "-L <int>            Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 20 both in --end-to-end mode and in --local mode.\n"
infoStr += "-L " + str(options.align_seedLength) + "\n\n"

infoStr += "-L <int>            Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 20 both in --end-to-end mode and in --local mode.\n"
infoStr += "-L " + str(options.align_seedLength) + "\n\n"

infoStr += "-D <int>            Up to <int> consecutive seed extension attempts can \"fail\" before Bowtie 2 moves on, using the alignments found so far. A seed extension \"fails\" if it does not yield a new best or a new second-best alignment. This limit is automatically adjusted up when -k or -a are specified. Default: 15.\n"
infoStr += "Default value in the mode very sensitive: 20\n"
infoStr += "-D " + str(options.align_maxAlignAttempts) + "\n\n"

infoStr += '-R <int>            <int> is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds. When "re-seeding," Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments. A read is considered to have repetitive seeds if the total number of seed hits divided by the number of seeds that aligned at least once is greater than 300. Default: 2.\n'
infoStr += "Default value in the mode very sensitive: 3\n"
infoStr += "-R " + str(options.align_maxReSeed) + "\n\n"

infoStr += """-i <func>
Sets a function governing the interval between seed substrings to use during multiseed alignment. For instance, if the read has 30 characters, and seed length is 10, and the seed interval is 6, the seeds extracted will be:

Read:      TAGCTACGCTCTACGCTATCATGCATAAAC
Seed 1 fw: TAGCTACGCT
Seed 1 rc: AGCGTAGCTA
Seed 2 fw:       CGCTCTACGC
Seed 2 rc:       GCGTAGAGCG
Seed 3 fw:             ACGCTATCAT
Seed 3 rc:             ATGATAGCGT
Seed 4 fw:                   TCATGCATAA
Seed 4 rc:                   TTATGCATGA
Since it's best to use longer intervals for longer reads, this parameter sets the interval as a function of the read length, rather than a single one-size-fits-all number. For instance, specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. See also: setting function options. If the function returns a result less than 1, it is rounded up to 1. Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 in --end-to-end mode to -i S,1,0.75 in --local mode.\n"""
infoStr += "-i " + options.align_seedInterval + "     Note: we use the same function as in the --very-sensitive preset, f(x) = 1 + 0.5 * sqrt(x)\n\n"

infoStr += """-X/--maxins <int>
The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the -X constraint is applied with respect to the untrimmed mates, not the trimmed mates.

The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences bewteen -I and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.

Default: 500.\n"""
infoStr += "-X " + str(options.align_maxInsertLength) + "\n\n"

infoStr += """-p/--threads NTHREADS
Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).\n"""
infoStr += "-p " + str(options.nThreads) + "\n\n"

infoStr += """-k <int>
By default, bowtie2 searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better. The best alignment found is reported (randomly selected from among best if tied). Information about the best alignments is used to estimate mapping quality and to set SAM optional fields, such as AS:i and XS:i.

When -k is specified, however, bowtie2 behaves differently. Instead, it searches for at most <int> distinct, valid alignments for each read. The search terminates when it can't find more distinct valid alignments, or when it finds <int>, whichever happens first. All alignments found are reported in descending order by alignment score. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM 'secondary' bit (which equals 256) set in its FLAGS field. For reads that have more than <int> distinct, valid alignments, bowtie2 does not guarantee that the <int> alignments reported are the best possible in terms of alignment score. -k is mutually exclusive with -a.

Note: Bowtie 2 is not designed with large values for -k in mind, and when aligning reads to long, repetitive genomes large -k can be very, very slow.\n"""
infoStr += "Note: if k=0, we do not include this option in the command.\n"
infoStr += "-k " + str(options.align_max_reported_alignments) + "\n\n"

pipelineDoc += infoStr


iTask += 1
task_name = 'align_seq'
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
taskPathList.append(task_path)
if options.seq_end == 'single-end':
    regexInputFiles = r'^(.+/)*(?P<SAMPLENAME>.+)\.trimmed\.fastq(\.gz)?.*'
elif options.seq_end == 'paired-end':
    regexInputFiles = r'^(.+/)*(?P<SAMPLENAME>.+?)_?read[12]\.trimmed\.fastq(\.gz)?$'
@follows(trim_adapter_PE_reads, trim_adapter_SE_reads, mkdir(str(task_path)))
@transform([trim_adapter_PE_reads, trim_adapter_SE_reads],

           # Match any of SE or PE trimmed reads fastq file
           formatter(regexInputFiles),

           # Output parameter supplied to next task
           str(task_path) + "/{SAMPLENAME[0]}.sam",

           # Extra parameters
           str(task_path) + "/{SAMPLENAME[0]}.summary",
           task_name, logger, logger_mutex)
def align_seq(input_files,
              sam_file,
              summary_file,
              task_name, logger, logger_mutex):

    if options.seq_end == 'single-end' and len(input_files) != 2:
        raise Exception("Number of input fastq files incorrect for single-end RNA-seq alignment, input_files: %s" % (input_files,))
    elif options.seq_end == 'paired-end' and len(input_files) != 3:
        raise Exception("Number of input fastq files incorrect for paired-end RNA-seq alignment, input_files: %s" % (input_files,))

    with logger_mutex:
        logger.debug("align_seq, input_files: " + ", ".join(input_files))

    if options.seq_end == 'single-end':
        readU = input_files[0]
        nValidReadsFilename = input_files[1]
    elif options.seq_end == 'paired-end':
        read1 = input_files[0]
        read2 = input_files[1]
        nValidReadsFilename = input_files[2]
        
    with open(nValidReadsFilename, 'r') as nreadsFile:
        nreads = int(re.match(r'^([0-9.]+).*', nreadsFile.readline()).group(1))
    with logger_mutex:
        logger.debug("align_seq, nreads " + str(nreads))

    ## Decompress the input file
    #with gzip.GzipFile(read1_file, mode='r') as gzip_file:
        #input_file_unzipped = re.search(r'(.+)\.gz', read1_file).group(1)
        #print(input_file_unzipped)
        #with open(input_file_unzipped, 'wb') as gunzip_file:
            #gunzip_file.write(gzip_file.read())
    #nReadsApprox = estimateNbLines(input_file_unzipped, 10*1024*1024) / 4.
    ## Unzipped fastq file can be deleted to save storage
    #os.remove(input_file_unzipped)
    # Approximate computation time in seconds per reads per thread measured for test run
    comp_time = ( (3.3 if options.align_nMismatches == 1 else 1)*
                  (2.7 if options.align_max_reported_alignments > 0 else 1)*
                  ((nreads*(((30.1)*16)/(2*2e5)))/options.nThreads) )

    if options.seq_end == 'single-end':
        bowtieInputOptions = '-U ' + readU
    elif options.seq_end == 'paired-end':
        bowtieInputOptions = " -1 {} -2 {} ".format(read1, read2)

    cmd = cmd_source_bash_profile +\
          "bowtie2 " +\
          " -x {} ".format(str(options.align_indexed_ref_genome_path)) +\
          bowtieInputOptions +\
          " -S {} ".format(sam_file) +\
          " --{} ".format(options.phredEncoding) +\
          " --{} ".format(options.align_alignmentMode) +\
          " -N " + str(options.align_nMismatches) +\
          " -L " + str(options.align_seedLength) +\
          " -i {} ".format(options.align_seedInterval) +\
          " -D {} ".format(options.align_maxAlignAttempts) +\
          " -R {} ".format(options.align_maxReSeed) +\
          " -X " + str(options.align_maxInsertLength) +\
          " -p " + str(options.nThreads)
    if options.align_max_reported_alignments > 0:
          cmd += " -k " + str(options.align_max_reported_alignments)

    with logger_mutex:
        logger.debug(cmd)
        logger.debug("estimated computation time " + str(datetime.timedelta(seconds=comp_time)))

    try:
        stdout_res, stderr_res = "",""
        walltime = datetime.timedelta(seconds=max(2.5*comp_time, 60))
        if walltime < datetime.timedelta(hours=6):
            job_queue_name = short_queue
        else:
            job_queue_name = long_queue
        job_other_options = " -pe smp " + str(options.nThreads) +\
                            " -q " + job_queue_name +\
                            " -l h_rt=" + printTimeDelta(walltime) +\
                            " -l h_vmem=8G,virtual_free=8G" +\
                            " -cwd"
        with logger_mutex:
            logger.debug("Submitting job, cmd:\n" + cmd + "\njob options:" + job_other_options)
        # ruffus.drmaa_wrapper.run_job
        stdout_res, stderr_res = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                         drmaa_session=drmaa_session, run_locally=run_locally,
                                         job_other_options=job_other_options, retain_job_scripts=False)

    # relay all the stdout, stderr, drmaa output to diagnose failures
    except error_drmaa_job as err:
        if options.sendMessageToWebServer:
            send_socket_message(options.analysisId, -1)
        raise Exception("\n".join(map(str,["Failed to run:", cmd, err, stdout_res, stderr_res])))

    std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
    with logger_mutex:
        logger.debug(std_err_string)

    with open(summary_file, 'w') as summaryFile:
        summaryFile.write(str(std_err_string))

    with logger_mutex:
        logger.debug("Alignment of reads finished.")



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """Convert alignment file SAM to sorted BAM.

Software: samtools Version: 1.3.1 (using htslib 1.3.1)

"""
pipelineDoc += infoStr


iTask += 1
task_name = 'convert_sam_to_bam'
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
taskPathList.append(task_path)
@follows(align_seq, mkdir(str(task_path)))
@transform(align_seq,

           # sam file
           formatter(r'^(.+/)*(?P<SAMPLENAME>.+)\.sam$'),

           # sorted bam file
           str(task_path) + "/{SAMPLENAME[0]}_sorted.bam",

           # Extra parameters
           "{SAMPLENAME[0]}",
           task_name, logger, logger_mutex)
def convert_sam_to_bam(sam_file,
                       sorted_bam_file,
                       sample_name,
                       task_name, logger, logger_mutex):

    # Remark: beware of the option -m XXXG that sets memory limit for samtools sort,
    # if a float is passed, it seems that samtools takes the value as bytes and will create
    # hundred of thousands of temporary files, potentially collapsing the filesystem.
    cmd = cmd_source_bash_profile +\
          " samtools view -b -h -u " + sam_file +\
          " | samtools sort -@ {:d} -m {:d}M -T {} -o {}".format(options.samtools_sort_nthread,
                                                                 options.samtools_sort_max_mem_per_thread,
                                                                 options.samtools_sort_tmp_dir,
                                                                 sorted_bam_file) +\
          " && samtools index " + sorted_bam_file
    with logger_mutex:
        logger.debug(cmd)

    try:
        stdout_res, stderr_res = "",""
        walltime = datetime.timedelta(hours=5)
        if walltime < datetime.timedelta(hours=6):
            job_queue_name = short_queue
        else:
            job_queue_name = long_queue
        job_other_options = " -pe smp " + str(options.samtools_sort_nthread) +\
                            " -q " + job_queue_name +\
                            " -l h_rt=" + printTimeDelta(walltime) +\
                            " -l virtual_free={}M".format(options.samtools_sort_max_mem)

        # ruffus.drmaa_wrapper.run_job
        stdout_res, stderr_res  = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                          drmaa_session=drmaa_session, run_locally=run_locally,
                                          job_other_options=job_other_options, retain_job_scripts=False)

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

        # relay all the stdout, stderr, drmaa output to diagnose failures
    except error_drmaa_job as err:
        if options.sendMessageToWebServer:
            send_socket_message(options.analysisId, -1)
        raise Exception("\n".join(map(str,["Failed to run:" + cmd, err, stdout_res, stderr_res])))

    with logger_mutex:
        logger.debug(task_name + " finished.")



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """

Filter alignments by quality and size and report statistics.

We filter reads that do not align concordantly on the genome (mostly very short reads that align almost randomly).

In order to detect properly the strand of the RNA reads, we compare the number of Forward reads and Reverse reads in the first rRNA region. The case with higher nb of reads is the + strand. We convert the BAM file to BEDPE file (BED with paired-end information with two read mates in the same record) or directly to BED file. Using the correct strand information, we run the python script pipeline_roesti_bedpe2bed.py to convert the paired-end reads to a single BED reads corresponding to the physical RNA fragment with strand correctly assigned. Then, we remove aligned reads that overlap with the rRNA or tRNA regions (optional).

Software: samtools Version: 1.3.1 (using htslib 1.3.1)
Software: bedtools v2.26.0
Software: sort (GNU coreutils) 8.4

"""
pipelineDoc += infoStr


iTask += 1
task_name = 'filter_alignments'
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
taskPathList.append(task_path)
@follows(convert_sam_to_bam, mkdir(str(task_path)))
@transform(convert_sam_to_bam,

           # sorted bam file
           formatter(r'^(.+/)*(?P<SAMPLENAME>.+)_sorted\.bam$'),

           # output files
           [str(task_path) + '/{SAMPLENAME[0]}.filtered.bed',
            str(task_path) + '/{SAMPLENAME[0]}.filtered.bed.nreads'],

           # Sample name
           "{SAMPLENAME[0]}",
           # Input path
           "{path[0]}",
           # Output path
           str(task_path),
           task_name, logger, logger_mutex)
def filter_alignments(sorted_bam_file,
                      output_files,
                      sample_name,
                      input_path,
                      output_path,
                      task_name, logger, logger_mutex):
       

    filter_script_filename = str(scriptPath / 'pipeline_roesti_filter_script.sh')

    cmd = cmd_source_bash_profile +\
          " " + filter_script_filename +\
          " " + sample_name +\
          " " + input_path +\
          " " + output_path +\
          " " + str(options.filter_alignments_quality_threshold) +\
          " false" +\
          " " + options.rRNA_bedfile +\
          " " + options.rRNA_tRNA_bedfile +\
          " " + str(scriptPath) +\
          " " + ("true" if options.remove_rRNA else "false") +\
          " " + options.seq_end
    with logger_mutex:
        logger.debug(cmd)

    try:
        stdout_res, stderr_res = "",""
        walltime = datetime.timedelta(hours=20)
        if walltime < datetime.timedelta(hours=6):
            job_queue_name = short_queue
        else:
            job_queue_name = long_queue
        job_other_options = " -pe smp " + str(8) +\
                            " -q " + job_queue_name +\
                            " -l h_rt=" + printTimeDelta(walltime) +\
                            " -l h_vmem=24G,virtual_free=24G"

        # ruffus.drmaa_wrapper.run_job
        stdout_res, stderr_res = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                         drmaa_session=drmaa_session, run_locally=run_locally,
                                         job_other_options=job_other_options, retain_job_scripts=False)

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

    # relay all the stdout, stderr, drmaa output to diagnose failures
    except error_drmaa_job as err:
        if options.sendMessageToWebServer:
            send_socket_message(options.analysisId, -1)
        raise Exception("\n".join(map(str,["Failed to run:" + cmd, err, stdout_res, stderr_res])))

    time.sleep(sleepTimeFilesystem)

    # Read results from samtools read counts
    nreads_filename = output_path + '/' + sample_name + '.bed.nreads'
    wait_for_file(nreads_filename, sleepTimeFilesystem)
    with open(nreads_filename) as nreads_file:
        nreads_bed = int(next(nreads_file).split()[0])
    os.remove(nreads_filename)

    # nreads_filename = output_path + '/' + sample_name + '.filtered.bed.nreads'
    nreads_filename = output_files[1]
    wait_for_file(nreads_filename, sleepTimeFilesystem)
    with open(nreads_filename) as nreads_file:
        nreads_bed_filtered = int(next(nreads_file).split()[0])

    # total nb of raw reads
    if options.seq_end == 'paired-end':
        nreads_raw_filename = output_path + '/../Task01_trim_adapter_PE_reads/' + sample_name + '_.trimmed.nreads'
        if not Path(nreads_raw_filename).exists():
            nreads_raw_filename = output_path + '/../Task01_trim_adapter_PE_reads/' + sample_name + '.trimmed.nreads'
    elif options.seq_end == 'single-end':
        nreads_raw_filename = output_path + '/../Task02_trim_adapter_SE_reads/' + sample_name + '_.trimmed.nreads'
        if not Path(nreads_raw_filename).exists():
            nreads_raw_filename = output_path + '/../Task02_trim_adapter_SE_reads/' + sample_name + '.trimmed.nreads'
    wait_for_file(nreads_raw_filename, sleepTimeFilesystem)
    with open(nreads_raw_filename) as nreads_file:
        nreads_raw = int(next(nreads_file).split()[0])

    # Write percentage of reads passing the filter in CSV file
    with open(output_path + '/' + sample_name + '.nreads','w') as nreads_file:
        nreads_file.write('nreads_raw,' + str(nreads_raw) + '\n')
        nreads_file.write('nreads_bed,' + str(nreads_bed) + '\n')
        nreads_file.write('nreads_bed_filtered,' + str(nreads_bed_filtered) + '\n')
        nreads_file.write('nreads_bed_filtered/nreads_bed' + str(nreads_bed_filtered/nreads_bed) + '\n')

    with logger_mutex:
        logger.debug(task_name + " finished.")



#############################################################################
#############################################################################
# FROM HERE WE FORK THE PIPELINE INTO TWO DIRECTIONS FOR RIBO-SEQ OR RNA-SEQ
#############################################################################
#############################################################################



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """
Extract ribosome footprints inserts by size.

Software: samtools Version: 1.3.1 (using htslib 1.3.1)
Software: bedtools v2.26.0

"""
pipelineDoc += infoStr


iTask += 1
task_name = 'extract_footprints'
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
# Only run this analysis for ribo-seq data type
if options.library_type == 'ribo-seq':
    task_path.mkdir(exist_ok=True)
    finalResultsPath = task_path
@active_if(options.library_type == 'ribo-seq')
# @follows(filter_alignments, mkdir(str(task_path)))
@transform(filter_alignments,

           # BED file of paired-end reads (inserts), nreads file
           formatter(r'^(.+/)*(?P<SAMPLENAME>.+)\.filtered\.bed.*$'),

           # output file
           str(task_path) + "/{SAMPLENAME[0]}.footprints_narrow.bed",

           # Sample name
           "{SAMPLENAME[0]}",
           # Input path
           "{path[0]}",
           # Output path
           str(task_path),
           task_name, logger, logger_mutex)
def extract_footprints(input_files,
                       footprints_file,
                       sample_name,
                       input_path,
                       output_path,
                       task_name, logger, logger_mutex):

    if re.match(r'^(.+/)*(.+)\.bed$', input_files[0]):
        reads_bed_file = input_files[0]
    if re.match(r'^(.+/)*(.+)\.bed$', input_files[1]):
        reads_bed_file = input_files[1]

    filter_script_filename = str(scriptPath / 'pipeline_roesti_extract_footprints.sh')

    cmd = cmd_source_bash_profile +\
          " " + filter_script_filename +\
          " " + sample_name +\
          " " + input_path +\
          " " + output_path +\
          " false" +\
          " " + str(scriptPath) +\
          " " + str(options.genomeBedFile)
    with logger_mutex:
        logger.debug(cmd)

    try:
        stdout_res, stderr_res = "",""
        walltime = datetime.timedelta(hours=20)
        if walltime < datetime.timedelta(hours=6):
            job_queue_name = short_queue
        else:
            job_queue_name = long_queue
        job_other_options = " -pe smp " + str(4) +\
                            " -q " + job_queue_name +\
                            " -l h_rt=" + printTimeDelta(walltime) +\
                            " -l h_vmem=16G,virtual_free=16G"

        # ruffus.drmaa_wrapper.run_job
        stdout_res, stderr_res = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                         drmaa_session=drmaa_session, run_locally=run_locally,
                                         job_other_options=job_other_options, retain_job_scripts=False)

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

    # relay all the stdout, stderr, drmaa output to diagnose failures
    except error_drmaa_job as err:
        if options.sendMessageToWebServer:
            send_socket_message(options.analysisId, -1)
        raise Exception("\n".join(map(str,["Failed to run:" + cmd, err, stdout_res, stderr_res])))


    nreads_filename = output_path + '/' + sample_name + '.footprints_wide.nreads'
    wait_for_file(nreads_filename, sleepTimeFilesystem)
    with open(nreads_filename) as nreads_file:
        nreads_footprints_wide = int(next(nreads_file).split()[0])
    os.remove(nreads_filename)

    nreads_filename = output_path + '/' + sample_name + '.footprints_narrow.nreads'
    wait_for_file(nreads_filename, sleepTimeFilesystem)
    with open(nreads_filename) as nreads_file:
        nreads_footprints_narrow = int(next(nreads_file).split()[0])
    os.remove(nreads_filename)

    wait_for_file(nreads_bed_file, sleepTimeFilesystem)
    with open(nreads_bed_file) as nreads_file:
        nreads_bed_filtered = int(next(nreads_file).split()[0])

    # Write percentage of reads passing the filter in CSV file
    with open(output_path + '/' + sample_name + '.footprints_nreads','w') as nreads_file:
        nreads_file.write('nreads_footprints_wide,' + str(nreads_footprints_wide) + '\n')
        nreads_file.write('nreads_footprints_wide/nreads_bed_filtered,' + str(nreads_footprints_wide/nreads_bed_filtered) + '\n')
        nreads_file.write('nreads_footprints_narrow,' + str(nreads_footprints_narrow) + '\n')
        nreads_file.write('nreads_footprints_narrow/nreads_bed_filtered,' + str(nreads_footprints_narrow/nreads_bed_filtered) + '\n')

    with logger_mutex:
        logger.debug(task_name + " finished.")



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """
Compute genome coverage and fragment count of reads (mRNA).

[script pipeline_roesti_genome_coverage.sh] We compute the fragment count per CDS using `bedtools intersect`, with strand specific option and setting the fractional overlap threshold to 0.5, meaning that at least half of the fragment has to overlap with the feature to be counted. We compute the per base coverage for plus and minus strand using `bedtools genomecov`.

[script pipeline_roesti_mean_coverage_per_CDS.py] We compute the average coverage per base for all CDS by integrating the per-base coverage along the feature and dividing by the length of the feature. This is in theory more precise than just counting the fragments, since fragment size distribution may vary between genes.
"""
pipelineDoc += infoStr


# Note that we need to sort the BED file, otherwise bedtools intersect will use huge memory (70G)
# when computing the number of reads per gene and is much slower.
# We also have to make sure that the CDS BED file is also sorted.
# we sort the CDS BED file from here the python ruffus script and not in the task function below
# because we only need to do it once, otherwise each parallel job will try to sort the CDS BED file,
# resulting in error at runtime.
script_filename = str(scriptPath / 'pipeline_roesti_sort_CDS_BED.sh')
cmd = cmd_source_bash_profile +\
      " " + script_filename +\
      " " + str(options.genomeCDSBedFile)
print(cmd)
cmd_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
cmd_output = re.sub(r'\\n','\n', str(cmd_output))
print(cmd_output)


iTask += 1
task_name = 'genome_coverage_fragment_count'
task_path = pipeline_path / "Task{:02d}_{}".format(iTask, task_name)
# Only run this analysis for rna data type
if options.library_type == 'rna-seq':
    task_path.mkdir(exist_ok=True)
    finalResultsPath = task_path
@active_if(options.library_type == 'rna-seq')
@follows(filter_alignments, mkdir(str(task_path)))
@transform(filter_alignments,

           # BED file of reads (inserts)
           formatter(r'^(.+/)*(?P<SAMPLENAME>.+)\.filtered\.bed$'),

           # output file
           str(task_path) + "/{SAMPLENAME[0]}.strandp_coverage.bed",

           # Sample name
           "{SAMPLENAME[0]}",
           # Input path
           "{path[0]}",
           # Output path
           str(task_path),
           task_name, logger, logger_mutex)
def genome_coverage_fragment_count(reads_bed_file,
                                   genome_coverage_file,
                                   sample_name,
                                   input_path,
                                   output_path,
                                   task_name, logger, logger_mutex):

    # For some reason the formatter filtering does not work and we have to apply a regex again to choose the .bed file only.
    if len(reads_bed_file) != 1:
        reads_bed_file = [fn for fn in reads_bed_file if re.search(r'^(.+/)*(?P<SAMPLENAME>.+)\.filtered\.bed$', fn)][0]
    # with logger_mutex:
    #     logger.debug(task_name + ", filtered input file: " + reads_bed_file)

    nreads_filename = sample_name + '.filtered.bed.nreads'
    nreads_filepath = Path(input_path) / nreads_filename
    if nreads_filepath.is_file():
        with nreads_filepath.open() as f:
            nreads = int(next(f).split()[0])
    else:
        nreads = 10e6
    print("sample", sample_name, "nreads", nreads)

    filter_script_filename = str(scriptPath / 'pipeline_roesti_genome_coverage.sh')

    cmd = cmd_source_bash_profile +\
          " " + filter_script_filename +\
          " " + reads_bed_file +\
          " " + sample_name +\
          " " + input_path +\
          " " + output_path +\
          " false" +\
          " " + str(scriptPath) +\
          " " + str(options.genomeBedFile) +\
          " " + str(options.genomeCDSBedFile) +\
          " " + str(options.rRNA_bedfile)
    with logger_mutex:
        logger.debug(cmd)

    try:
        stdout_res, stderr_res = "",""
        walltime = datetime.timedelta(hours=20)
        if walltime < datetime.timedelta(hours=6):
            job_queue_name = short_queue
        else:
            job_queue_name = long_queue
        # Large number of reads might require a fairly large memory for the bedtools intersect,
        # even if the reads are sorted. Last example case was 55M reads and required 16G of memory.
        if nreads > 30e6:
            memory = '24G'
        else:
            memory = '16G'
        job_other_options = " -pe smp " + str(1) +\
                            " -q " + job_queue_name +\
                            " -l h_rt=" + printTimeDelta(walltime) +\
                            " -l virtual_free={}".format(memory)

        # ruffus.drmaa_wrapper.run_job
        stdout_res, stderr_res = run_job(cmd_str=cmd, job_name=task_name, logger=logger,
                                         drmaa_session=drmaa_session, run_locally=run_locally,
                                         job_other_options=job_other_options, retain_job_scripts=False)

        std_err_string = "".join([line.decode() if isinstance(line, bytes) else line for line in stderr_res])
        with logger_mutex:
            logger.debug(std_err_string)

    except error_drmaa_job as err:
        if options.sendMessageToWebServer:
            send_socket_message(options.analysisId, -1)
        raise Exception("\n".join(map(str,["Failed to run:" + cmd, err, stdout_res, stderr_res])))

    with logger_mutex:
        logger.debug(task_name + " finished.")

    finalResultsPath = task_path



#############################################################################
infoStr = "\n\n#############################\n"
infoStr += """

PCR duplicates.
Note: In the ribosome profiling data analysis, we do not try to remove PCR duplicates. PCR duplicates can only be detected with high confidence when the probability of having in the library RNA inserts (fragments) of same length mapping to the same exact position in the genome is very small. In the ribosome footprints library, fragments length is distributed around 30 nt, and we expect many footprints to map to the same position. Alternative approach to detect PCR duplicates include adding short randomized nucleotides at both ends of fragments during library preparation. See Lecanda, A., Nilges, B. S., Sharma, P., Nedialkova, D. D., Schwarz, J., Vaquerizas, J. M., & Leidel, S. A. (2016). Dual randomization of oligonucleotides to reduce the bias in ribosome-profiling libraries. Methods, 10–12.

Note: In the RNA-seq data analysis, we do not remove PCR duplicates neither.


"""
pipelineDoc += infoStr


#############################################################################
iTask += 1
task_name = 'delete_intermediate_files'


@follows(genome_coverage_fragment_count)
@active_if(options.delete_intermediate_files)
def delete_intermediate_files():

    print("Delete intermediate files.")
    for p in taskPathList:
        for f in p.glob('*'):
            print("deleting file:", f)
            f.unlink()
        Path(p).rmdir()

    with logger_mutex:
        logger.debug(task_name + " finished.")



#############################################################################
iTask += 1
task_name = 'write_jobid_files'


@follows(delete_intermediate_files)
def write_jobid_files():

    # Check that final output files exist
    # task_path is the last task's output directory
    fileList0 = [f for f in Path(task_path).iterdir() if f.is_file()]
    comments = ''
    status = 2

    # We just check that the coverage file for plus and minus strands exist
    # and sum up at least one read.
    for strand in ['plus', 'minus']:
        if strand == 'plus':
            regex = re.compile(r'^(.+/)*(?P<SAMPLENAME>.+)\.strandp_coverage\.bed$')
        elif strand == 'minus':
            regex = re.compile(r'^(.+/)*(?P<SAMPLENAME>.+)\.strandm_coverage\.bed$')
        fileList = [f for f in fileList0 if regex.search(f.name)]
        if len(fileList) > 0:
            file = fileList[0]
            df = pd.read_table(str(file))
            nReads = df.iloc[:, 2].sum()
            if nReads < 1:
                status = -1
                comments = comments + 'strand {} coverage file sums up 0 reads.\n'.format(strand)
        else:
            status = -1
            comments = comments + 'strand {} coverage file not found.\n'.format(strand)

    # We also test that the CDS_values file exists, though we do not consider
    # an error if it does not exists. User could only want to compute the coverage,
    # with no annotation count values.
    regex = re.compile(r'^(.+/)*(?P<SAMPLENAME>.+)\.CDS_values\.csv$')
    fileList = [f for f in fileList0 if regex.search(f.name)]
    if len(fileList) == 0:
        status = 2
        comments = comments + 'CDS_values.csv file not found.\n'

    # filename = 'pipeline_done.{:d}.txt'.format(options.jobid)
    # with (outputPath / filename).open('w') as f:
    #     f.write('')

    if options.sendMessageToWebServer:
        send_socket_message(options.analysisId, status)

#############################################################################


pipelineDocFile.write_text(pipelineDoc)

# Change the history file for Ruffus in order to use several different pipelines in the same root folder
history_file = "." + pipeline_name + ".ruffus_history.sqlite"
options.history_file = history_file
if run_on_cluster:
    pipeline_printout(history_file=history_file)
    cmdline.run(options, multithread=options.njobs, logger=logger, verbose=options.verbose)
    drmaa_session.exit()
