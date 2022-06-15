#!/usr/bin/perl

# Andrew Tolonen

# Script to read Bowtie output (sam file) to count the number of 5' read ends at each position in the C.phy genome. Output is a text file with 3 columns: genome position, strand, number of 5' read ends.

###

# input files to modify for each sample

#sam file input
$inputfile = "reads_mmlvTrimmed.sam";

# start position output
$outputfile = "reads_allStartPositions.csv";

$readID = "M2:"; # sequence that identifies read names
#$threshold = 10; # min number of reads per position to print to output
###

open(IN, "<$inputfile") or die "cant open bowtie-aligned sam file\n";
open(OUT, ">$outputfile") or die "cant open start position output file\n";

# read bowtie input. build hash: genome position -> [strand, #reads]

$totalreads = 0;
$fiveprimeends = 0;

LOOP:while ($line = <IN>)
{
 chomp($line);
 if ($line =~ /^$readID/)
 {
  @fields = split("\t", $line);
  $strand = $fields[1]; # forward or reverse strand
  if ($strand == 4) # read has no alignments
  {
   next LOOP;
  }
  $position = $fields[3]; # position of leftmost character on forward strand
  $seq = $fields[9];
  $seqlength = length($seq); # length of aligned sequence
  unless (($strand == 0) or ($strand == 16)) 
  {
   print "read $fields[0] failed because it maps to neither forward nor reverse strand\n";
   exit;
  }
  if ($strand == 16) # sequence maps to reverse strand
  {
   $position = $position + $seqlength; # calculate 5' end of reverse alignments
  } 
  $id = "$position"."_"."$strand";
  $fiveprime{$id}++;
  $totalreads++;
 }
}
close IN;  

#print "Found $totalreads reads with 5' ends in input file $inputfile.\n";

print OUT "Genome_Position\tStrand\tReads\n";

foreach $x (keys(%fiveprime))
{
 $x =~ /(\d+)_(\d+)/;
 $end = $1;
 $strand = $2;
 if ($strand == 0)
 {
  $strand = "+";
 }
 if ($strand == 16)
 {
  $strand = "-";
 }
# if ($fiveprime{$x}>$threshold)
# {
  print OUT "$end\t$strand\t$fiveprime{$x}\n";
  $fiveprimeends++;
# }
}
close OUT;

#print "Input file contains $fiveprimeends genome positions with more than $threshold mapped 5' ends.\n";
print "$inputfile contains $fiveprimeends genome positions with a total of $totalreads mapped 5' ends.\n";

