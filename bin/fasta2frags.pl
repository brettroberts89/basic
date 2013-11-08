#TITLE: fasta2frags.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 7/9/13
#PURPOSE: splits a fasta sequence into 100 bp fragments

#!/usr/bin/perl -w
use strict;

my $fasta = $ARGV[0];
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
close FASTA;

open OUT, ">$ARGV[1]";

foreach my $seq (keys %sequences) {
    my $sequence = $sequences{$seq};
    my @index = (0..(length($sequence)-100));
    foreach my $pos (@index) {
        my $read = substr($sequence, $pos, 100);
        print OUT ">$pos\n$read\n\n";
    } 
}









