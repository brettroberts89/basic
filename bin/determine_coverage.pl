#TITLE: determine_coverage.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 02/21/13
#PURPOSE: takes in genome length and fastq file with reads, calculates coverage of genome

#!/usr/bin/perl -w
use strict;

my $genome_length = $ARGV[0];
my $fastq = $ARGV[1];
my $total_bases = 0;

open FASTQ, "<$fastq";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
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
    $total_bases += length($seq);
}

my $coverage = $total_bases/$genome_length;

print "Coverage is $coverage X.";