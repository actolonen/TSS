#!/bin/bash

# Andrew Tolonen

# Script to map RNA-seq reads to genome using Bowtie. Forward direction reads are trimmed of first 3 bp. Reverse direction reads are untrimmed.

# files to modify to run bowtie script
trimmedForward=reads_mmlvTrimmed.fastq
sam=reads_mmlvTrimmed.sam

# load modules
module load bowtie2/2.1.0
echo loaded $BOWTIE2_EXEDIR/bowtie2

# bowtie database
bowtieDB=PATH # Path to bowtie database for reference genome

# run bowtie 
$BOWTIE2_EXEDIR/bowtie2 \
-x $bowtieDB \
-U $trimmedForward \
-S $sam
 
