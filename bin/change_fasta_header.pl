#!/usr/bin/perl -w

#TITLE: change_fasta_header.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/16/2013
#PURPOSE: changes fasta headers

use strict;
use warnings;

open FASTA, "<$ARGV[0]";
my $contig;
my $sequence;
my %sequences;
my $count = 1;

while (my $line = <FASTA>) {                  
    chomp $line;
    if ($line =~ /^>/) {
        $count++;
        $contig = 'chr' . $count;
    }
    else{             
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} .= "$line";
    }
}

foreach my $a (sort keys %sequences) {
    print ">$a\n$sequences{$a}\n";
}