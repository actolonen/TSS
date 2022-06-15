#!/usr/bin/perl

# Andrew Tolonen

# Script to add annotations to classified TSS.

# input and output files.

open(TSS, "<tss_classified.csv") or die "Can't open classified TSS\n";
open(ANN, "<genome.ptt") or die "can't open annotation file\n";
open(OUT, ">tss_classified_annotation.csv") or die "Can't open out\n";

###
# grab annotation

while ($line = <ANN>)
{
 chomp($line);
 if ($line =~ /Cphy_[\d]{4}/)
 {
  @fields = split("\t", $line);
  $gene = $fields[5];
  $annotation = $fields[8];
  $ann{$gene} = $annotation;
 }
}
close ANN;

$counter = 0;
foreach $x (keys(%ann))
{
 $counter++;
}
print "Grabbed annotation for $counter genes\n";

###

# add annotation to classified TSS

print OUT "Genome_position\tStrand\tReads\tClass\tGene\tGeneStart\tDistance\tAnnotation\n";
LOOP:while ($line = <TSS>)
{
 chomp($line);
 if ($line =~ /(Cphy_[\d]{4})/)
 {
  $gene = $1;
  foreach $y (keys(%ann))
  {
   if ($gene eq $y)
   {
    print OUT "$line\t$ann{$y}\n";
    next LOOP; 
   }
  }
  print OUT "$line\tNo annotation found: putative pseudogene\n";
 }
}

