#!/bin/bash

# Andrew Tolonen

# script to make bam of Capp-Switch reads using sam file as input  

###

# 4 files to customize to run script 

# sam file made by bowtie
sam=reads_mmlvTrimmed.sam

# bam file output
bam=reads_mmlvTrimmed.bam

# sorted version of bam file
sortedbam=reads_mmlvTrimmed_sorted
sortedbamin=reads_mmlvTrimmed_sorted.bam

###

# print run time
date

# load samtools
module load samtools/1.2
echo loaded $SAMTOOLS_EXEDIR/samtools

# convert sam to bam

$SAMTOOLS_EXEDIR/samtools \
view -bS $sam \
-o $bam

# sort the bam file

# read 1
$SAMTOOLS_EXEDIR/samtools \
sort $bam \
$sortedbam

# make a bam index file

$SAMTOOLS_EXEDIR/samtools \
index $sortedbamin

# cleanup sam and unsorted bam files
#/bin/rm $sam
/bin/rm $bam
