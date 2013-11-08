#TITLE: calculate_gap_lengths.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/15/13
#PURPOSE: get avg and stdev of gap lengths in scaffolds

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $output = $input;
$output =~ s/\.final\.scaffolds\.fasta/_gaplengths.txt/;

my $contig;
my $sequence;
my %sequences;

open IN, "<$input";
while (my $line = <IN>) {                 
    chomp $line;
    if ($line =~ /^>/) {
        $contig = $line;
        $contig =~ s/>//;
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = $sequences{$contig} . $sequence;
    }
}

open OUT, ">$output";

foreach my $scaffold (keys %sequences) {
    my @sequence = split(//,$sequences{$scaffold});
    my $prev_base = 0;
    my $gap_length = 0;
    foreach my $base (@sequence) {
        if ($prev_base =~ /N/i) {
            if ($base =~ /N/i) {
                $gap_length += 1;
            }
            else {
                print OUT "$gap_length\n";
                print "$scaffold,$gap_length\n";
                $gap_length = 0;
            }
        }
        else {
            if ($base =~ /N/i) {
                $gap_length += 1;
            }
            else {
                next;
            } 
        }
        $prev_base = $base;
    }
}