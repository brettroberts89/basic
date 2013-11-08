#TITLE: 6-2_polymorphisms.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2012
#PURPOSE: finds polymorphisms in EG6-2 in the all.mergedSNPs.vcf file, and outputs the location of these SNPs

#!/usr/bin/perl -w
use strict;

my $input = "jatropha_ApeKI_ref.txt";
my $vcf = "all.mergedSNPs.vcf";
my $output = "6-2_polymorphisms.csv";
my %contigs;
my %snps;
my %polymorphs;

open VCF, "<$vcf";
open IN, "<$input";
open OUT, ">$output";

while (my $line = <IN>) {                                                       #from index file, create hash with keys contigs and values of the start and stop positions with the offset
    chomp $line;
    my ($contig, $start_offset, $end_offset) = split(/\D+/, $line);
    my $length = $end_offset - $start_offset;
    $contigs{$contig} = "$start_offset\t$end_offset\n";
    #print OUT "$contig,$start_no_offset,$end_no_offset\n";
}

print OUT "Coordinate,Contig,Contig_Coordinate,Reference_Base,EG6_2_Genotype\n";            #interpret vcf notation to list alleles present in EG6-2
while (my $line = <VCF>) {
    chomp $line;
    if ($line =~ /^1/) {
        my @line = split(/\t/, $line);
        my $coord = $line[1];
        my $ref_base = $line[3];
        my $alt_base = $line[4];
        my @alt_base = split /,/,$alt_base;
        my $six_two = $line[10];
        my @six_two = split /\D+/, $six_two;
        my $base_a;
        my $base_b;
        if (($six_two[0] != "0") or ($six_two[1] != "0")) {
            $base_a = $ref_base if ($six_two[0] == "0");
            $base_a = $alt_base[0] if ($six_two[0] == "1");
            $base_a = $alt_base[1] if ($six_two[0] == "2");
            $base_a = $alt_base[2] if ($six_two[0] == "3");
            $base_b = $ref_base if ($six_two[1] == "0");
            $base_b = $alt_base[0] if ($six_two[1] == "1");
            $base_b = $alt_base[1] if ($six_two[1] == "2");
            $base_b = $alt_base[2] if ($six_two[1] == "3");
            foreach my $contig (keys %contigs) {
                my ($contig_start, $contig_end) = split /\t/, $contigs{$contig};
                if (($coord >= $contig_start) and ($coord <= $contig_end)) {
                    my $contig_coord = $coord - $contig_start + 1;
                    print OUT "$coord,$contig,$contig_coord,$ref_base,$base_a/$base_b\n";
                    last;
                }
            }
        }
    }
}