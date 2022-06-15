#!/usr/bin/perl

# Andrew Tolonen

# Script to classify TSS into 1 of 4 categories: intragenic sense, intragenic antisense, intergenic sense (downstream gene is same orientation), or intergenic antisense.

###

# input and output files

# input of clustered 5' read positions: genome position, strand, reads
open(INA, "<tss_replicateComparison.csv") or die "cant open input reads\n";

# input of gene boundaries
open(GENES, "<genome.gff") or die "cant open gene boundary input\n";

# output file of classified TSS: genome position, strand, reads (avg reads per million), class, associated gene, distance to gene feature
open(OUT, ">tss_classified.csv") or die "cant open out\n";
###

# subroutines
GetInfo();
CompareReads();

###

sub GetInfo
{

# make %reads hash: genomePosition_strand -> [sample 1 reads, sample 2 reads, sample 1 reads per million, sample 2 reads per million, average reads per million]

while ($line = <INA>)
{
 chomp($line);
 if ($line =~ /^\d/)
 {
  @fields = split("\t", $line);
  $key = "$fields[0]"."_"."$fields[1]"; # hash key is genome position + strand orientation.
  $sample1 = $fields[2]; # reads in sample 1 
  $sample2 = $fields[3]; # reads in sample 2
  $reads{$key} = [$sample1, $sample2]; 
 }
}
close INA;


# count the total number of reads in each sample
$totalsample1=0;
$totalsample2=0;
$numpos = 0;
foreach $x (keys(%reads))
{
 $totalsample1 = $totalsample1 + $reads{$x}[0];
 $totalsample2 = $totalsample2 + $reads{$x}[1];
 $numpos++;
}
print "$numpos total positions in both samples\n";
print "$totalsample1 total reads in sample 1; $totalsample2 reads in sample 2\n";

# calc reads per million foreach TSS and append to hash 
foreach $x (keys(%reads))
{
 $sample1_norm = 1e6*($reads{$x}[0]/$totalsample1);
 $sample2_norm = 1e6*($reads{$x}[1]/$totalsample2);
 $avg_norm = ($sample1_norm + $sample2_norm) / 2;
 push(@{$reads{$x}}, $sample1_norm);
 push(@{$reads{$x}}, $sample2_norm);
 push(@{$reads{$x}}, $avg_norm);
}

# make %genes hash: gene->[strand, start, stop]

$genecounter=0;
LOOP:while ($line = <GENES>)
{
 chomp($line);
 if ($line =~ /^NC_010001.1/)
 {
  @fields = split("\t", $line);
  if ($fields[2] eq "gene")
  {
   if ($fields[8] =~ /locus_tag=(Cphy_[\d]{4})/)
   {
    $genecounter++;
    $genename = $1;
    $strand = $fields[6];
    if ($strand eq "+")
    {
     $start = $fields[3];
     $stop = $fields[4];
     $genes{$genename} = [$strand, $start, $stop];
     next LOOP;
    }
    if ($strand eq "-")
    {
     $start = $fields[4];
     $stop = $fields[3];
     $genes{$genename} = [$strand, $start, $stop];
     next LOOP;
    }
    print "strand error on line:\n $line\n";exit;
   }
  }
 }
}
close GENES;
print "Genome contains $genecounter genes\n";
}
###

sub CompareReads
{

# compare each read to gene hash and append to reads hash: 1. gene that contains (intragenic) or is downstream of read (intergenic). 2. classification: intragenic_sense, intragenic_antisense, intergenic_sense, intergenic_antisense

LOOP:foreach $x (keys(%reads))
{
# grab info about read
 $x =~ /(\d+)_([+-])/;
 $genomepos = $1;
 $strand = $2;
 $numreads = $reads{$x}[0]; 
 $readcounter = $readcounter + $numreads;

# test if read is intragenic by comparing read to each gene
 foreach $y (keys(%genes))
 {
  if ($genes{$y}[0] eq "+") # if gene on pos strand
  {
   if (($genomepos >= $genes{$y}[1]) and ($genomepos <= $genes{$y}[2])) # if read is >start position and <stop position
   {
    if ($strand eq "+") # if gene and reads in same orientation
    { 
     $class = "intragenic_sense";
     push(@{$reads{$x}}, $class); # append read class to reads hash
     push(@{$reads{$x}}, $y); # append gene name 
     push(@{$reads{$x}}, $genes{$y}[1]); # append gene start
     $dist = $genomepos - $genes{$y}[1]; # distance to ATG
     push(@{$reads{$x}}, $dist); # append distance to start
     next LOOP;
    }
    if ($strand eq "-")  # read is negative, gene is positive
    {
     $class = "intragenic_antisense";
     push(@{$reads{$x}}, $class); # append reads class
     push(@{$reads{$x}}, $y); # append gene name
     push(@{$reads{$x}}, $genes{$y}[1]); # append gene start
     $dist = $genomepos - $genes{$y}[1]; # distance to ATG
     push(@{$reads{$x}}, $dist); # append distance to start
     next LOOP;
    }
   }
  }
  if ($genes{$y}[0] eq "-") # if gene on neg strand
  {
   if (($genomepos >= $genes{$y}[2]) and ($genomepos <= $genes{$y}[1]))
   {
    if ($strand eq "-")
    { 
     $intraS = $intraS +$numreads;
     $class = "intragenic_sense";
     push(@{$reads{$x}}, $class); # append read class to reads hash
     push(@{$reads{$x}}, $y); # append gene to reads hash
     push(@{$reads{$x}}, $genes{$y}[1]); # append gene start
     $dist = $genes{$y}[1] - $genomepos; # distance to ATG
     push(@{$reads{$x}}, $dist); # append distance to reads hash
     next LOOP;
    }
    if ($strand eq "+") # read is positive, gene is negative 
    {
     $intraA = $intraA + $numreads;
     $class = "intragenic_antisense";
     push(@{$reads{$x}}, $class); # append class to reads hash
     push(@{$reads{$x}}, $y); # append gene to reads hash
     push(@{$reads{$x}}, $genes{$y}[1]); # append gene start
     $dist = $genes{$y}[1] - $genomepos; # distance to ATG
     push(@{$reads{$x}}, $dist); # append distance to reads hash
     next LOOP;
    }
   }
  }
 }
# for intergenic reads, get the name and distance to downstream gene
 $dist = 1000000; # initialize distance to nearest read to a big number
 $neargene = "none";
 $class = "none";
 
 if ($strand eq "+") # if reads on pos strand
 {
  foreach $y (keys(%genes)) # compare read to all genes
  {
   if (($genes{$y}[0] eq "+") and ($genes{$y}[1] > $genomepos)) # if gene on pos strand and downstream of TSS
   {
    $tempdist = $genes{$y}[1] - $genomepos; # distance from start codon to read
    $tempclass = "intergenic_sense";
    if ($tempdist < $dist)
    {
     $dist = $tempdist; # update distance to nearest gene
     $neargene = $y; # update gene name
     $class = $tempclass; # update class
     $genestart = $genes{$y}[1];
    }
   }
   if (($genes{$y}[0] eq "-") and ($genes{$y}[2] > $genomepos)) # of gene on neg strand and downstream of TSS
   {
    $tempdist = $genes{$y}[2] - $genomepos;
    $tempclass = "intergenic_antisense";
    if ($tempdist < $dist) # if gene is closest yet tested
    {
     $dist = $tempdist; # update distance to nearest gene
     $neargene = $y; # update gene name
     $class = $tempclass; # update class
     $genestart = $genes{$y}[1];
    }
   }
  }
  push(@{$reads{$x}}, $class); # append gene class to reads hash
  push(@{$reads{$x}}, $neargene); # append gene name to reads hash
  push(@{$reads{$x}}, $genestart); # append gene start 
  push(@{$reads{$x}}, $dist);
  next LOOP;
 }

 if ($strand eq "-") # if TSS on neg strand
 {
  foreach $y (keys(%genes))
  {
   if (($genes{$y}[0] eq "-") and ($genes{$y}[1] < $genomepos)) # if gene on pos strand and downstream of TSS
   {
    $tempdist = $genomepos - $genes{$y}[1];
    $tempclass = "intergenic_sense";
    if ($tempdist < $dist) # if gene is closest yet tested
    {
     $dist = $tempdist; # update distance to nearest gene
     $neargene = $y; # update gene name
     $class = $tempclass; # update class
     $genestart = $genes{$y}[1];
    }
   }
   if (($genes{$y}[0] eq "+") and ($genes{$y}[2] < $genomepos)) # of gene on neg strand and downstream of TSS
   {
    $tempdist = $genomepos - $genes{$y}[2];
    $tempclass = "intergenic_antisense";
    if ($tempdist < $dist) # if gene is closest yet tested
    {
     $dist = $tempdist; # update distance to nearest gene
     $neargene = $y; # update gene name
     $class = $tempclass; # update class
     $genestart = $genes{$y}[1];
    }
   }
  }
  push(@{$reads{$x}}, $class); # append gene name to reads hash
  push(@{$reads{$x}}, $neargene); # append gene name to reads hash
  push(@{$reads{$x}}, $genestart); # append gene start 
  push(@{$reads{$x}}, $dist);
  next LOOP;
 }
}

# print output
$readcounter=0;
$intraS=0;
$intraA=0;
$interS=0;
$interA=0;

print OUT "Genome_Position\tStrand\tReadsPerMillion\tClass\tGene\tGeneStart\tDistance\n";

foreach $x (keys(%reads))
{
 $x =~ /(\d+)_([+-])/;
 $genomepos = $1;
 $strand = $2;
 $gene = $reads{$x}[6];
 $class = $reads{$x}[5];
 $numreads = $reads{$x}[4];
 $genestart = $reads{$x}[7]; 
 $dist = $reads{$x}[8];
 $readcounter = $readcounter + $numreads;
 if ($class eq "intergenic_sense")
 {
  $interS = $interS + $numreads; # increment read count
  $interSG{$gene}++;
 }
 if ($class eq "intergenic_antisense")
 {
  $interA = $interA + $numreads;
  $interAG{$gene}++;
 }
 if ($class eq "intragenic_sense")
 {
  $intraS = $intraS + $numreads;
  $intraSG{$gene}++;
 }
 if ($class eq "intragenic_antisense")
 {
  $intraA = $intraA + $numreads;
  $intraAG{$gene}++;
 }
 print OUT "$genomepos\t$strand\t$numreads\t$class\t$gene\t$genestart\t$dist\n";
}

$interSkeys = 0;
foreach $x (keys(%interSG))
{
 $interSkeys++;
}

$interAkeys = 0;
foreach $x (keys(%interAG))
{
 $interAkeys++;
}

$intraSkeys = 0;
foreach $x (keys(%intraSG))
{
 $intraSkeys++;
}

$intraAkeys = 0;
foreach $x (keys(%intraAG))
{
 $intraAkeys++;
}

# print number of genes and reads in each category
print "intergenic_sense: $interSkeys genes and $interS total reads per million\n";
print "intergenic_antisense: $interAkeys genes and $interA total reads per million\n";
print "intragenic_sense: $intraSkeys genes and $intraS total reads per million\n";
print "intragenic_antisense: $intraAkeys genes and $intraA total reads per million\n";
}
