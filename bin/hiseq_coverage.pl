#TITLE: hiseq_coverage.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 05/11/2013
#PURPOSE: determines the perentage of the genome covered by illumina reads

#!/usr/bin/perl -w
use strict;

my $genome = $ARGV[0];
open GENOME, "<$genome";

my $contig;
my $sequence;
my %sequences;
my $genome_length = 0;

while (my $line = <GENOME>) {                  
    chomp $line;
    if ($line =~ /^>/) {                  
        $contig = $line;
        $contig =~ s/>//;
    }
    else{             
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
        my $length = length($sequence);
        $genome_length += $length;
    }
}

my $uncovered = 0;
my $bedgraph = $ARGV[1];

open BEDGRAPH, "<$bedgraph";
while (<BEDGRAPH>) {
    if ($_ =~/0$/) {
        my @line = split /\t/, $_;
        my $scaffold = $line[0];
        my $min = $line[1];
        my $max = $line[2];
        for (my $i = $min; $i <= $max; $i++) {
            next if (substr($sequences{$scaffold}, $i, 1) eq 'N');
            $uncovered += 1;
        }
    }
}

my $uncovered_pct = $uncovered/$genome_length;
print "$uncovered reads uncovered out of $genome_length for uncovered proportion of $uncovered_pct\n";