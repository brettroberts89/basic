#!/usr/bin/perl -w

#TITLE: first_100bp.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/03/2013
#PURPOSE: gets first 100 bp of sequence, outputs in primer3 input format

use strict;
use warnings;

open FASTA, "<$ARGV[0]";
my $contig;
my $sequence;
my %sequences;

while (my $line = <FASTA>) {                  
    chomp $line;
    if ($line =~ /^>/) {                  
        $contig = $line;
        $contig =~ s/>//;
        $contig=~ s/\|.*//;
    }
    else{             
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} .= "$line";
    }
}

foreach my $seq (keys %sequences) {
    my $front = substr($sequences{$seq}, 0, 100);
    print "SEQUENCE_ID=$seq\nSEQUENCE_TEMPLATE=$front\n=\n";
}