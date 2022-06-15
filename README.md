# Transcription Start Site (TSS) identification from Capp-Switch reads
This pipeline enables the identification, classification, and annotation of TSSs from Illumina sequencing reads, as published in [Boutard et al, 2016](https://www.nature.com/articles/ncomms13783).  

## Step 1: remove the MMLV tail from FASTQ reads

As [shown here](https://www.nature.com/articles/ncomms13783/figures/1), the MMLV reverse transcriptase appends a 5' tag to the sequencing reads, which needs to removed in order for the reads to map correctly to the reference genome file.

01_remove_MMLV_tail.plx takes the FASTQ file as input and outputs a modified FASTQ with the 5' MMLV tag removed.

## Step 2: map reads to reference genome

A trimmed FASTQ file with the MMLV tags removed is mapped to the reference genome.

02_bowtie2.bash takes the trimmed FASTQ file and maps the reads to the reference genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). The output is a SAM file of mapped reads.

## Step 3: make BAM file of mapped reads

A SAM file of mapped reads is converted to a BAM file. BAMs are binary files that retain the same information as the SAM, but save disk space.

03_makeBam.bash converts the SAM to BAM using [SAM tools](http://www.htslib.org/doc/samtools-view.html).

## Step 4:  

