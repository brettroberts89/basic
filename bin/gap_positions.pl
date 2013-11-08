#TITLE: gap_positions.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/15/13
#PURPOSE: get avg and stdev of gap lengths in scaffolds

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $sam = $ARGV[1];
my $output = $input;
my $sanger_primers = 'Sanger_Forward_primers.txt';


my $contig;
my $sequence;
my %sequences;
my %sanger;

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

my @gaps;

foreach my $scaffold (keys %sequences) {
    my @sequence = split(//,$sequences{$scaffold});
    my $prev_base = 0;
    my $gap_length = 0;
    my $gap_start;
    foreach my $pos (0..(scalar(@sequence)-1)) {
        my $base = $sequence[$pos];
        if ($prev_base =~ /N/i) {
            if ($base =~ /N/i) {
                $gap_length += 1;
            }
            else {
                push(@gaps, "$gap_start:$pos");
                $gap_length = 0;
            }
        }
        else {
            if ($base =~ /N/i) {
                $gap_start = $pos;
            }
            else {
                next;
            } 
        }
        $prev_base = $base;
    }
}

open SAM, "<$sam";

my @sequence = split(//,$sequences{'SGBv1_scaffold0_SGBv9'});

open SANGER, "<$sanger_primers";

while (my $line = <SANGER>) {
    chomp $line;
    $sanger{$line} = 1;
}

while (my $line = <SAM>) {
    chomp $line;
    my @line = split(/\t/, $line);
    foreach my $gap (@gaps) {
        my ($gap_start, $gap_stop) = split (/\:/, $gap);
        if (($line[3] < $gap_start) and ($line[7] > $gap_stop)) {
            if (exists $sanger{$line[0]}) {
                print ">$line[0]\n";
                my @amplicon = @sequence[$line[3]..$line[7]];
                my $amplicon = join('',@amplicon);
                print"$amplicon\n";
            }
        }
    }
}
