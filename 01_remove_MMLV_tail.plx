#!/usr/bin/perl

# Andrew Tolonen

# Script to process fastq files: If the first 3 bp correspond to the MMLV tail
# then the tail is removed and the read is kept. Otherwise, the read is rejected.


$inputfile = "reads.fastq";
$outputfile = "reads_mmlvTrimmed.fastq";

open(IN, "<$inputfile") or die "cant open fastq input file";
open(OUT, ">$outputfile") or die "cant open fastq output file";

# remove the 12 most frequent 5' 3-mers that account for 95% of sequence
@mmlv = ("GGG", "CGG", "GGC", "GCG", "CCG", "GGT", "GCC", "CGC", "CGT", "GTG", "CCC", "GCT");

$reads = 0;
$printedreads = 0;
while ($line = <IN>)
{ 
 if ($line =~ /^\@/) # read up to seq ID line
 {
  $reads++; # increment total reads found
  chomp($line);
  $seqID = $line; # save seqID line
  $line = <IN>; # grab sequence line
  chomp($line);
  $sequence = $line; # save sequence line
  $three = substr($sequence, 0, 3);  # grab first 3 bp of sequence line
  $line = <IN>; # grab plus line
  chomp($line); 
  $plus = $line; # save "+" line
  $line = <IN>; # grab quality line
  chomp($line);
  $qual = $line; # save quality line
  foreach $x (@mmlv)
  {
   if ($three eq $x) # if first 3 bases match an mmlv sequence
   {
    $sequence = substr($sequence, 3); # remove first 3 bases
    $qual = substr($qual, 3); # remove first 3 quality scores
    print OUT "$seqID\n$sequence\n$plus\n$qual\n";
    $printedreads++; # increment reads that start with mmlv sequence
   }
  }
 }
}
$percent = 100*($printedreads/$reads);
print "input file: $inputfile\n";
print "file contain $reads total reads\n";
print "$printedreads ($percent%) reads started with an MMLV sequence. These reads were trimmed and printed to output file.\n";
