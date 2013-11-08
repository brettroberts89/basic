#!/usr/bin/perl -w

#TITLE: extract_cladek_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/27/2013
#PURPOSE: takes in list of individuals and fastq and outputs only fastq entries for reads originating from the individuals in the list.

use strict;
use warnings;

my $USAGE = "USAGE: perl extract_cladek_reads.pl INDIVS FASTQ OUTPUT\n";

if (!$ARGV[0] or !$ARGV[1]) {
    die $USAGE;
}

my $indivs = $ARGV[0];

open INDIVS, "<$indivs";

my @indivs;
my %handles;

while (my $line = <INDIVS>) {
    chomp $line;
    push @indivs, $line;
}

open FASTQ, "<$ARGV[1]";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @seq_id = split /:|\s/, $line_1;
    my $indiv = $seq_id[-1];
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    foreach my $a (@indivs) {
        if ($indiv eq $a) {
            my $outfile = $indiv . "_CLADEK.fastq";
            open OUT, ">>$outfile";
            print OUT "$line_1\n$seq\n$plus\n$quality\n";
            last;
        }
    }
}