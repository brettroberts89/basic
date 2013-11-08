#TITLE: onetenth_fa.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 11/28/12
#PURPOSE: creates file 1/10th the size of the original fasta file for testing.

#!/usr/bin/perl -w
use strict;
my $contig;
my $sequence;
my %sequences;
my $count = 0;

my $input = $ARGV[0];
my $output = $input;
$output =~ s/\.fasta//;
$output = $output . 'small.fasta';

open IN, "<$input";
open OUT, ">$output";


while (my $line = <IN>) {                   #go through Sato assembly, creating a hash entry for each contig with the sequence as the value
    chomp $line;
    my $prev_contig;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $prev_contig = $contig;
        $contig = $line;
        $count += 1;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
    if ($count == 1) {
        print OUT "$prev_contig\n$sequences{$prev_contig}\n";
    }
    elsif ($count == 10) {
        $count = 0;
    }
    else{
        next;
    }
}