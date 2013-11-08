#TITLE: snp_coords_list.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/06/2012
#PURPOSE: converts vcf to list coordinates of SNPs with contigs as chromosome number and coordinate relative to start of contig as coordinate
#!/usr/bin/perl -w
use strict;

my $vcf = "all.mergedSNPs.vcf";
my $output = "all.mergedSNPs_coords.vcf";
my $offset = 0;
my %contigs;
my $prev_end = 0;


open OUT, ">$output";


open VCF, "<$vcf";
while (my $line = <VCF>) {
    chomp $line;
    my $output_line;
    if ($line =~ /^\d/) {                        #if the line begins with a 1, which is the case for everything but the header
        my @line = split(/\t/, $line);          #split the line into an array at any value that is not a digit
        my $contig = $line[0];
        my $coordinate = $line[1];              #the coordinate is the second value in the array
        print OUT "$contig,$coordinate\n";
}