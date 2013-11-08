#TITLE: snp_coords_list_indu.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/19/2012
#PURPOSE: obtains list of Locus, contig, and coordinate relative to the start of the contig for all SNPs in a vfc file.


#!/usr/bin/perl -w
use strict;

my $vcf = "all.mergedSNPs_SGB.vcf";
my $output = "india_snps.csv";
my %contigs;


open OUT, ">$output";


open VCF, "<$vcf";
while (my $line = <VCF>) {
    chomp $line;
    my $output_line;
    if ($line =~ /^\d/) {                        #only read if the line begins with a digit, which is the case for all lines but the header
        my @line = split(/\t/, $line);          #split the line into an array dividing by tabs
        my $contig = $line[0];
        my $coordinate = $line[1];         
        my $id = $line[2];
        my $sample = $line[12];
        ##my $ = $line[12];
        my @sample = split /\D+/, $sample;
        next if ((($sample[0] == 0) and ($sample[1] == 0)) or ($sample =~ /\./));
        print OUT "$id,$contig,$coordinate\n";
    }
}