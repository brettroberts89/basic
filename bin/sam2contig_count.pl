#!/usr/bin/perl -w

#TITLE: sam2contig_count.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/27/2013
#PURPOSE: takes in sam file and outputs tab delimited list of all contigs with reads aligning and number of reads aligned

use strict;
use warnings;
use List::Util qw(sum);

my $USAGE = "Usage: perl sam2contig_count.pl SAM";

if (!$ARGV[0]) {
    die $USAGE;
}

my %scaffold_counts;

foreach my $sam (@ARGV) {
    open SAM, "<$sam";
    while (my $line = <SAM>) {
        chomp $line;
        next if ($line =~ /^\@SQ/);
        next if ($line =~ /^\@PG/);
        my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
        $scaffold_counts{$scaffold} += 1;
    }
}

foreach my $scaffold (keys %scaffold_counts) {
    print "$scaffold\t$scaffold_counts{$scaffold}\n"
}