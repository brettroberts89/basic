#!/usr/bin/perl -w

#TITLE: extract_fasta_seqs.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/29/2013
#PURPOSE: makes new file containing only specified sequences from fasta file.

use strict;
use warnings;

my $USAGE = "USAGE: perl extract_fasta_seqs.pl FASTA CSV OUTPUT\n";

unless ($ARGV[0] and $ARGV[1] and $ARGV[2]) {
    die $USAGE;
}

my $fasta = $ARGV[0];
my $csv = $ARGV[1];
my $output = $ARGV[2];

open FASTA, "<$fasta";

my $contig;
my $sequence;
my %sequences;

while (my $line = <FASTA>) {
    chomp $line;
    if ($line =~ /^>/) {
        $contig = $line;
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

open CSV, "<$csv";
open OUT, ">$output";

my %printed;

while (my $line = <CSV>) {
    chomp $line;
    my ($seq_id1,$linkage_group,$pos,$scaffold1,$scaffold2,$scaffold1_length,$scaffold2_length,$pos1,$pos2,$mult_mapping1,$mult_mapping2) = split(/,/, $line);
    next if ($scaffold1 !~ /scaffold/);
    foreach my $contig (keys %sequences) {
        if ($contig =~ /$scaffold1\|/) {
            next if (exists $printed{$contig});
            print "$scaffold1\n";
            print OUT ">$contig\n$sequences{$contig}\n";
            $printed{$contig} = 1;
            last;
        }
    }
}