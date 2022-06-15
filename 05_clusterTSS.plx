#!/usr/bin/perl

# Andrew Tolonen

# Script to cluster TSS. All TSS with more than $mincounts reads that are within $clustersize bp are clustered. Only the genome position with the highest # reads is retained. 
# Input is text file genome position, strand, number of 5' read ends.  
###

# input and output files to modify for each sample.
$inputfile = "reads_allStartPositions.csv";
$outputfile = "reads_allStartPositions_clustered.csv";

###

$clustersize = 5; # cluster TSS within 5 bp

getStarts();
clusterStarts();
sortStarts();

#######

sub getStarts
{

# count the total number of reads to define the min reads threshold

open(IN, "<$inputfile") or die "cant open in\n";

$allreads = 0;
while ($line = <IN>)
{
 chomp($line);
 if ($line =~ /^\d/)
 {
  @fields = split("\t", $line);
  $reads = $fields[2];
  $allreads = $allreads + $reads;
 }
}
close IN;

$threshold = 10*($allreads/1000000);
print "Sample contains $allreads total reads, set min reads threshold to $threshold reads.\n";

###

# read positions with >threshold 5' starts into two hashes, one for forward stand and one for reverse strand: position -> #reads

open(IN, "<$inputfile") or die "cant open in\n";


$readcounter=0;
while ($line = <IN>)
{
 chomp($line);
 if ($line =~ /^\d/)
 {
  @fields = split("\t", $line);
  $position = $fields[0];
  $strand = $fields[1];
  $reads = $fields[2];
  if (($strand eq "+") and ($reads > $threshold))
  {
   $fiveprimeplus{$position} = $reads;
   $readcounter = $readcounter + $reads;
  }
  if (($strand eq "-") and ($reads > $threshold))
  {
   $fiveprimeneg{$position} = $reads; 
   $readcounter = $readcounter + $reads;
  }
 }
}
close IN;

# count number of genome positions with read starts before clustering. 

$pluspos = 0; # count keys in %fiveprimeplus
foreach $x (keys(%fiveprimeplus))
{
 $pluspos++;
}

$negpos = 0; # count keys in %fiveprimeneg
foreach $x (keys(%fiveprimeneg))
{
 $negpos++;
}

print "Genome positions with <$threshold reads were removed:\nplus-strand=$pluspos positions.\nminus-strand=$negpos positions.\ntotal reads=$readcounter.\n";

}

#######

sub clusterStarts
{

# cluster TSS on forward strand

$maxpos = 5000000;

for ($testpos = 1; $testpos<=$maxpos; $testpos++) # test all potential tss positions in genome
{
 if (exists($fiveprimeplus{$testpos})) # if there are reads defined at that position
 {
  $testreads = $fiveprimeplus{$testpos};
  $minustwo=0;
  $minusone=0;
  $plusone=0;
  $plustwo=0;
  if (exists($fiveprimeplus{$testpos-2}))
  {
   $minustwo = $fiveprimeplus{$testpos-2};
  }
  if (exists($fiveprimeplus{$testpos-1}))
  {
   $minusone = $fiveprimeplus{$testpos-1};
  }
  if (exists($fiveprimeplus{$testpos+1}))
  {
   $plusone = $fiveprimeplus{$testpos+1};
  }
  if (exists($fiveprimeplus{$testpos+2}))
  {
   $plustwo = $fiveprimeplus{$testpos+2};
  }
  @counts = ($testreads, $minustwo, $minusone, $plusone, $plustwo); 
  @sortedcounts = sort { $a <=> $b } @counts;
  $maxcounts = pop(@sortedcounts); 
  unless ($testreads == $maxcounts)
  {
   delete($fiveprimeplus{$testpos});
  }
 } 
}

# cluster TSS on reverse strand

for ($testpos = 1; $testpos<=$maxpos; $testpos++) # test all potential tss positions in genome
{
 if (exists($fiveprimeneg{$testpos})) # if there are reads defined at that position
 {
  $testreads = $fiveprimeneg{$testpos};
  $minustwo=0;
  $minusone=0;
  $plusone=0;
  $plustwo=0;
  if (exists($fiveprimeneg{$testpos-2}))
  {
   $minustwo = $fiveprimeneg{$testpos-2};
  }
  if (exists($fiveprimeneg{$testpos-1}))
  {
   $minusone = $fiveprimeneg{$testpos-1};
  }
  if (exists($fiveprimeneg{$testpos+1}))
  {
   $plusone = $fiveprimeneg{$testpos+1};
  }
  if (exists($fiveprimeneg{$testpos+2}))
  {
   $plustwo = $fiveprimeneg{$testpos+2};
  }
  @counts = ($testreads, $minustwo, $minusone, $plusone, $plustwo); 
  @sortedcounts = sort { $a <=> $b } @counts;
  $maxcounts = pop(@sortedcounts); 
  unless ($testreads == $maxcounts)
  {
   delete($fiveprimeneg{$testpos});
  }
 } 
}

}

#######

sub sortStarts
{
open(OUT, ">$outputfile") or die "cant open out\n";

# sort and print positions after clustering

$pluscounter=0;
$negcounter=0;
$readcounter = 0;

print OUT "GenomePosition\tStrand\tReads\n";

LOOP:for ($i = 1; $i<=$maxpos; $i++) # range of genome positions with TSS
{ 
 if (exists($fiveprimeplus{$i}))
 { 
  $pluscounter++;
  $readcounter = $readcounter+$fiveprimeplus{$i};
  print OUT "$i\t+\t$fiveprimeplus{$i}\n"; 
 } 
 if (exists($fiveprimeneg{$i}))
 { 
  $negcounter++;
  $readcounter = $readcounter+$fiveprimeneg{$i};
  print OUT "$i\t-\t$fiveprimeneg{$i}\n"; 
 } 
} 
$allcounter = $pluscounter + $negcounter;
print "After clustering:\n$pluscounter read start positions on the plus-strand.\n$negcounter positions on the minus-strand.\ntotal positions=$allcounter.\ntotal reads=$readcounter.\n";

}
