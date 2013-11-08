#TITLE: fasta_gc_content.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 06/21/13
#PURPOSE: calculates the gc % for each sequence in a fasta

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

foreach my $transcript (keys %sequences) {
    my $gc = 0;
    my $length = 0;
    my @sequence = split(//, $sequences{$transcript});
    foreach my $base (@sequence) {
        $gc += 1 if ($base =~ /g|c/i);
        $length += 1;
    }
    my $gc_pct = $gc/$length;
    print OUT "$transcript\t$gc_pct\n";
}