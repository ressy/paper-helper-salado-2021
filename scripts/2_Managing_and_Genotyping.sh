#!/bin/bash

#### MANAGING_and_GENOTYPING_GENOTYPING.sh ####
# Bash script to concatenate fastq files by Sample-PCR and creating folders by
# treatment (Nreads) where concatenated fastq files are moved and where Megasat
# will be run afterwards

## Creating three variables: ##

# dir = the name of the directory where subsampled fastq files are saved;
dir="Sub_output"

# FILENAME__PCR = a list which keeps part of filenames until PCR to rename
# concatenated fastqs; Name: sub-Nreads*-Rep*-Sample-PCR
FILENAME_PCR=($(find "$dir" -type f -name '*fastq' | rev | cut -d'/' -f1 \
	| rev | cut -d'.' -f1 | cut -d'_' -f1 | uniq))

# DIRNAME_REP = a list which keeps part of filenames until Replicate to rename
# directories. Name sub-Nreads*-Rep*
DIRNAME_REP=($(find "$dir" -maxdepth 1 -type f -name '*fastq'| rev \
	| cut -d'/' -f1 | rev | cut -d'.' -f1 | cut -d'_' -f1 |cut -d'-' -f1-3 \
	| sort -u))

## Concatenating subsampled fastq files by filename (by PCR): ##
# Creating a new directory to save concatenated fastq files mkdir "$dir"/concatenated
for i in "${FILENAME_PCR[@]}"; do
	cat /dev/null "$dir"/"$i"* > "$dir"/concatenated/"${i[@]}.fastq"
done

## Creating a directory for each treatment and replicate, and moving
# concatenated fastq files (by PCR) to their corresponding directory: ##

for j in "${DIRNAME_REP[@]}"; do
	mkdir "$dir"/concatenated/"${j[@]}"
	mv $(find "$dir"/concatenated/ -name "${j[@]}-*fastq") \
		"$dir"/concatenated/"${j[@]}"
	# Renaming fastq files preparing format of Genotype.txt (output file) of
	# Megasat:
	rename 's/sub-Nreads\d+-Rep\d+-//' "$dir"/concatenated/"${j}"/*.fastq
	# Running Megasat for each treatment:
	perl MEGASAT_Genotype.pl primer_input_megasat_Mic_coverage.csv 2 5 1 \
		"$dir"/concatenated/"${j}" "$dir"/concatenated/"${j}"
done

## Creating a new directory to copy Megasat's results on it and renaming each
# file with the name of the treatment and repetition (according to its
# container folder) ##

# Creating a new directory to save Megasat results
mkdir -p Megasat_output || exit 1
inpath='Your path to directory /concatenated with concatenated and subsampled fastq files'
for l in "$inpath"/*/*/Genotype.txt; do
	cp "$l" "Megasat_output/$( basename "$( dirname "$l")" )_Genotype.txt"
done
