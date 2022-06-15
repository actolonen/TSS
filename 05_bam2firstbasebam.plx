#!/usr/bin/perl

# Laurence Ettwiller (modified by Andrew Tolonen)

use strict;
use warnings;
#use Getopt::Long qw(GetOptions);

# Laurence Ettwiller (modified by Andrew Tolonen)

#this program is primary used for visualizing the TSS on IGV or other genome browser with the read being narrow down to 1 bp resolution. 
#input=bam file X.bam (preferably generated through local alignment such as Bowtie --local)
#output= X_start.bam and X_start.bam.bai

###

# input variables: 1. bam file 2. genome file (fasta). 3. library-type (F)
my $bamfile = "reads_mmlvTrimmed_sorted.bam";
my $genome = "genome.fna.fai"; #  file showing the size of the chromosomes for the organism for which your BED files are based
my $lib_type = "F";
my $samdir = "samtools-1.2/el6-x86_64-generic/bin";
my $beddir = "bedtools-2.24.0/el6-x86_64-generic/bin";

###
 
# load modules
my $loadsamtools = "module load samtools/1.2";
`$loadsamtools`;
print "loaded $samdir/samtools\n";

my $loadbedtools = "module load bedtools/2.24.0";
`$loadbedtools`;
print "loaded $beddir/bedtools\n";

###


# grab relevant reads according to library type.

my $resulting_bam;
if ($lib_type eq "F") #single read library with R1 being the most 5' end of the transcripts. 
{
    $resulting_bam = $bamfile;
}
elsif ($lib_type eq "FR")
{
    my $tmp = "R1";
    #extract R1 from bam
    my $library_command = "$samdir/samtools view -f64 -b $bamfile | $samdir/samtools sort - $tmp";
    $resulting_bam = $tmp.".bam";
    system($library_command);
}
elsif ($lib_type eq "RF")
{
    my $tmp = "R2";
    #extract R1 from bam                                                                                                                                                      
    my $library_command = "$samdir/samtools view -f128 -b $bamfile | $samdir/samtools sort - $tmp";
    $resulting_bam = $tmp.".bam";
    system($library_command);
}

###

my $generic = $bamfile;
$generic =~ s/\.bam//; # remove '.bam' from file name
$generic =~ s/.*\///g; # delete anything before a slash (get file name). 

print STDERR "Generating the bam files - be patient it may take a while - \n";
my $file_tmp = "bamtmp";
my $bed = "bedtmp";
my $newbam = $generic."_start";

#first convert your bam to bed. 
my $command = "$beddir/bamToBed -cigar  -i $resulting_bam > $file_tmp"; 
print STDERR "$command\n";
system($command);
parse_bed($file_tmp, $bed); # $file_tmp is bedtools output, $bed is temp file of parsed bedtools output

#then go back to bam again - sort it at the same time. 

my $command2 = "$beddir/bedToBam -i $bed -g $genome | $samdir/samtools sort - $newbam";
my $newbam_withbam = $newbam.".bam";
my $command3 = "$samdir/samtools index $newbam_withbam";

#excecute the commands

system($command2);
system($command3);

#removing tmp files that are not necessary anymore. 
unlink($file_tmp);
unlink($bed);

###
 
sub parse_bed {
    my ($file, $bed)=@_; # read input files: bedfile (text of bamfile) and bedtools output file
    open (FILE, $file) or die; # open bam file
    #write the first mapped postion into a bed new bed file 
    open (OUT, ">$bed") or die; 
    
    foreach my $line (<FILE>)
    {
	chomp $line;
	my $start;
	my @tmp = split /\t/, $line;
	my $chr = $tmp[0];
	my $orientation = $tmp[5];
	#if the orientation is + take the start of the read
	if ($orientation eq "+")
	{
	    $start = $tmp[1];
	    $tmp[2] = $start+1;
	}
	#if the orientation is - , take the end of the read
	else{
	    $start = $tmp[2];
	    $tmp[1] = $start-1;
	}
	#fixe up some other issues. 
	$tmp[4] =1;
	$tmp[6] ="1M";
	
	my $new_line = join("\t", @tmp);
	print OUT "$new_line\n";
    }
    close FILE;
    close OUT;
}

