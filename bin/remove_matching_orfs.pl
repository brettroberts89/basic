#TITLE: split_mapped_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 01/16/2012
#PURPOSE: splits reads into categories based on whether or not they mapped to the genome. Outputs file with list of sequence id's

#!/usr/bin/perl -w
use strict;

my $blast_out = $ARGV[0];
my $fasta = $ARGV[1];
my %match;
my $output = $fasta;
$output =~ s/\.fa/_filtered\.fa/;

open BLAST, "<$blast_out";
while (my $line =<BLAST>) {
    chomp $line;
    my @line = split(/\t/, $line);
    if ($line[2] == 100) {
        if (($line[7] - $line[6] + 1) == $line[3]) {
            $match{$line[0]} = 1;
        }
    }
}

my $contig;
my $sequence;
my %sequences;

open FASTA, "<$fasta";

while (my $line = <FASTA>) {                 
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig = $line[0];
        $contig =~ s/>//;
        #print "$contig\n";
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        if (exists $sequences{$contig}) {
            $sequences{$contig} = $sequences{$contig} . $sequence;
        }
        else {
            $sequences{$contig} = $sequence;
        }
    }
}

open OUT, ">$output";
foreach my $orf (keys %sequences) {
    if (exists $match{$orf}) {
        next;
    }
    else {
        print "$orf\n";
        print OUT ">$orf\n$sequences{$orf}\n";
        #print ">$orf\n$sequences{$orf}\n";
    }
}