# Transcription Start Site (TSS) identification from Capp-Switch reads
This pipeline enables the identification, classification, and annotation of TSSs from Illumina sequencing reads sequencin using Capp-Switch sequencing, as published in [Boutard et al, 2016](https://www.nature.com/articles/ncomms13783).  

## Step 1: Remove the MMLV tail from FASTQ reads

As [shown here](https://www.nature.com/articles/ncomms13783/figures/1), the MMLV reverse transcriptase appends a 5' tag to the sequencing reads, which needs to removed in order for the reads to map correctly to the reference genome file.

01_remove_MMLV_tail.plx takes the FASTQ file as input and outputs a modified FASTQ with the 5' MMLV tag removed.

## Step 2: Map reads to reference genome

A trimmed FASTQ file with the MMLV tags removed is mapped to the reference genome.

02_bowtie2.bash takes the trimmed FASTQ file and maps the reads to the reference genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). The output is a SAM file of mapped reads.

## Step 3: Make BAM file of mapped reads

A SAM file of mapped reads is converted to a BAM file. BAMs are binary files that retain the same information as the SAM, but save disk space.

03_makeBam.bash converts the SAM to BAM using [SAM tools](http://www.htslib.org/doc/samtools-view.html).

## Step 4: Count 5' read starts

The number of potential TSS (5' read start positions) are counted in the SAM file. 

04_countStarts.plx counts the number of TSS in the SAM file and outputs a table: Genome position, Strand (+/-), number of 5' read starts.

## Step 5: Cluster read starts to identify TSS

5' read starts that are highly clustered likely correspond to a single TSS. Thus, we identified highly clustered read starts and call the TSS as the position with the greatest number of reads starts.

05_clusterTSS.plx identified all TSS with >$mincounts reads that are <= $clustersize bp within range of each other. Only the genome position with the highest # reads is retained.

## Step 6: Compare replicates

Sequencing from replicate cultures is needed to reliably identify TSS. This step compares the TSS from replicate cultures.

06_compareReplicates.plx compares TSS from two replicate cultures. Input is the file of clustered TSS. Output is file of all positions with >10 reads/million reads in both samples.

## Step 7: Classify TSS

As shown in [panel B](https://journals.asm.org/doi/10.1128/spectrum.02288-21?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#fig3), TSS are classified into 1 of 4 categories:
* InterS (intergenic TSS with downstream gene in same orientation)
* InterA (intergenic TSS with downstream gene opposite orientation)
* IntraS (intragenic TSS in gene with same orientation)
* IntraA (intragenic TSS in gene with opposite orientation).

07_classifyTSS.plx takes alist of TSS, and uses the .gff file to one of the 4 categories shown above.

## Step 8: Annotate TSS

To aid in TSS interepretation, the annotation of the gene associated with the TSS is added.

08_addAnnotation.plx takes the classified TSS file as input and uses the genome .ppt file to add gene annotations.

## Step 9: Visualize TSS

The BAM file is revised to better visualize TSS using a tool such as IGV.

09_visualizeTSS.plx revised the BAM file to be compatible with IGV.


09_visualizeTSS.plx

