#!/usr/bin/perl
#
# Andrew Tolonen
#
# Script to compare 5' start positions for 2 Capp-Switch samples
# input=file of clustered start positions
# output=file of all positions with >10 reads/million reads in both samples.
# input and output files
$sample1 ="reads_allStartPositions_clustered_1.csv";
$sample2 ="reads_allStartPositions_clustered_2.csv";
$outputfile = "tss_replicateComparison.csv";

open(IN1, "<$sample1") or die "cant open sample 1 input\n";
open(IN2, "<$sample2") or die "cant open sample 2 input\n";
open(OUT, ">$outputfile") or die "cant open output file\n";

# 
# read sample 1. make hash for forward and reverse strand. genomePosition -> [sample1 reads, sample 2 reads]

while ($line = <IN1>)
{
 chomp($line);
 unless ($line =~ /^\d/)
 {
  next;
 }
 @fields = split("\t", $line);
 $pos = $fields[0]; # set hash ID as genome position 
 $strand = $fields[1]; # grab strand
 if ($strand eq "+")
 {
  $plusstarts{$pos}[0] = $fields[2]; # first hash elt is strand
  $plusstarts{$pos}[1]=0; #initialize to zero
 }
 if ($strand eq "-")
 {
  $negstarts{$pos}[0] = $fields[2]; # first hash elt is strand
  $negstarts{$pos}[1]=0; #initialize to zero
 }
}

close IN1;

# append start positions for sample 2

LOOP:while ($line = <IN2>)
{
 chomp($line);
 unless ($line =~ /^\d/)
 {
  next;
 }
 @fields = split("\t", $line);
 $pos = $fields[0]; # grab genome position 
 $strand = $fields[1]; # grab strand
 if ($strand eq "+")
 {
  foreach $x (keys(%plusstarts))
  {
   if ($pos == $x) # if genome position exists in sample 1
   {
    $plusstarts{$x}[1] = $fields[2];
    next LOOP;
   }
  }
  $plusstarts{$pos}[0] = 0;
  $plusstarts{$pos}[1] = $fields[2];
 }
 if ($strand eq "-")
 {
  foreach $x (keys(%negstarts))
  {
   if ($pos == $x) # if genome position exists in sample 1
   {
    $negstarts{$x}[1] = $fields[2];
    next LOOP;
   }
  }
  $negstarts{$pos}[0] = 0;
  $negstarts{$pos}[1] = $fields[2];
 }
}  
close IN1;

# count number of reads and positions in each sample

$readsS1=0;
$readS2=0;
$posS1=0;
$posS2=0;

foreach $x (keys(%plusstarts))
{
 if ($plusstarts{$x}[0] > 0) # if +strand reads in sample 1
 {
  $posS1++;
  $readsS1 = $readsS1 + $plusstarts{$x}[0]; 
 }
 if ($plusstarts{$x}[1] > 0) # if +strand reads in sample 2 
 {
  $posS2++;
  $readsS2 = $readsS2 + $plusstarts{$x}[1]; 
 }
}

foreach $x (keys(%negstarts))
{
 if ($negstarts{$x}[0] > 0) # if -strand reads in sample 1
 {
  $posS1++;
  $readsS1 = $readsS1 + $negstarts{$x}[0]; 
 }
 if ($negstarts{$x}[1] > 0) # if -strand reads in sample 2 
 {
  $posS2++;
  $readsS2 = $readsS2 + $negstarts{$x}[1]; 
 }
}

print "Sample 1: $posS1 positions, $readsS1 reads.\n";
print "Sample 2: $posS2 positions, $readsS2 reads.\n";

# print output genome positions found in both samples. 

$positions = 0;
$sample1reads = 0;
$sample2reads = 0;

print OUT "GenomePosition\tStrand\tSample1\tSample2\n";
foreach $x (keys(%plusstarts))
{
 if (($plusstarts{$x}[0]>0) and ($plusstarts{$x}[1]>$0)) # if position found in both samples
 { 
  print OUT "$x\t+\t$plusstarts{$x}[0]\t$plusstarts{$x}[1]\n";
  $positions++; 
  $sample1reads = $sample1reads + $plusstarts{$x}[0];
  $sample2reads = $sample2reads + $plusstarts{$x}[1];
 }
}
foreach $x (keys(%negstarts))
{
 if (($negstarts{$x}[0]>0) and ($negstarts{$x}[1]>0))
 {
  print OUT "$x\t-\t$negstarts{$x}[0]\t$negstarts{$x}[1]\n";
  $positions++; 
  $sample1reads = $sample1reads + $negstarts{$x}[0];
  $sample2reads = $sample2reads + $negstarts{$x}[1];
 }
}
 
print "Output $positions genome positions with $sample1reads sample 1 reads and $sample2reads sample 2 reads\n";

