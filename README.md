# Transcription Start Site (TSS) identification from Capp-Switch reads
This pipeline enables the identification, classification, and annotation of TSSs from Illumina sequencing reads, as published in [Boutard et al, 2016](https://www.nature.com/articles/ncomms13783).  

## Step 1: remove the MMLV tail from FASTQ reads

As [shown here](https://www.nature.com/articles/ncomms13783/figures/1), the MMLV reverse transcriptase appends a 5' tag to the sequencing reads, which needs to removed in order for the reads to map correctly to the reference genome file.

01_remove_MMLV_tail.plx takes the FASTQ files as input and outputs a modified FASTQ with the 5' MMLV tag removed.
