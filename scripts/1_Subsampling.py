#!/usr/bin/env python
#-*-coding:utf-8-*-d

#### SUBSAMPLING.PY ####

# Python script to subsample randomly without replacement the number of next
# generation sequencing reads obtained after a PCR amplification. Nreads is the
# number of total reads to subsample (treatments: 6,10,20,40,80,100,20O).
# Replicates are the number of repetitions for each treatment (100). Output
# files are presented in files with prefix 'sub-NreadsX-RepX-'. Afterwards,
# this script also calls to other additional scripts to genotype and analyze
# the resulting subsampled fastq files.

## Importing python modules and packages ##

import sys
import os
import random
import subprocess
from random import seed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Setting the paths where all input files and output files will be stored. ##

folder_path_in = sys.argv[1] # 'Your path to directory /Data with all fastq files '
folder_path_out = sys.argv[2] # 'Your path to directory /Sub_output to save new subsampled fastq files'

## Subsampling fastq files: ##

# Setting number of total reads to subsample:
Nreads = [6,10,20,40,80,100,150,200]
# Creating a folder "Sub_output" to save subsampled fastq files:
os.mkdir(folder_path_out)

for i in Nreads:
    replicates = 100
    seed(1)

    # Creating replicates of random number for each file(replicates=100):
    for j in range(1,replicates+1):
        # Opening files from the path
        for fastq in os.listdir(folder_path_in):
            if fastq.endswith(".fastq"):
                ruta = folder_path_in + '/' + fastq
                seq_list = list(SeqIO.parse(ruta, "fastq"))
                totseq = len(seq_list)
                # random.sample () subsamples sequences without replacement
                choice = random.sample(seq_list, i)
                #Writing the output list to a new subsampled fastq file
                ruta_out = folder_path_out + '/' + 'sub-' + 'Nreads' + str(i) \
                + '-' + 'Rep' + str(j) + '-' + fastq
                SeqIO.write(choice, ruta_out, "fastq")

## Organizing files and running Megasat (calling script in bash): ##
subprocess.call("2_Managing_and_Genotyping.sh", shell=True)

## Comparing of Megasat results and exporting final results (calling script in R):##
subprocess.call("Rscript 3_Comparing_genotypes.R", shell=True)
