#!/usr/bin/perl -w

#TITLE: blast_hits_fastq.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/23/2013
#PURPOSE: gets all reads that hit from blast tabular file and makes fastq with just those reads

use strict;
use warnings;

my %hits;

open BLAST, "$ARGV[0]";
my $fastq = $ARGV[1];
my $output = $ARGV[2];

while (<BLAST>) {
    chomp;
    my @line = split(/\t/, $_);
    my $read_id = $line[0];
    $hits{$read_id} = 1;
}

open FASTQ, "<$fastq";
open OUT, ">$output";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @line1 = split(/\s/, $line_1);
    my $seq_id = $line1[0];
    $seq_id =~ s/^\@//;
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
    if (exists $hits{$seq_id}) {
        print OUT "$line_1\n$seq\n$plus\n$quality\n";
    }
}


