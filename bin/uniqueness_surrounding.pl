#!/usr/bin/perl -w

#TITLE: vcf_uniqueness.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/31/2013
#PURPOSE: searches region 100 bp on either side of variants specified in vcf file to calculate uniqueness of region

use strict;
use warnings;

my $fasta = $ARGV[0];
my $vcf = $ARGV[1];

open FASTA, "<$fasta";
my $contig;
my $sequence;
my %sequences;

while (my $line = <FASTA>) {
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig = $line[0];
        $contig =~ s/>//;
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequence =~ s/\W//g;
        if (exists $sequences{$contig}) {
            $sequences{$contig} .= $sequence;
        }
        else {
            $sequences{$contig} = $sequence;
        }
    }
}

open VCF, "<$vcf";

while (<VCF>) {
    next if ($_ =~ /^\#/);
    chomp;
    my ($chrom, $pos, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/, $_);
    my $left_seq = substr($sequences{$chrom}, ($pos-101), 100);
    my $right_seq = substr($sequences{$chrom}, $pos, 100);
    
    
}





