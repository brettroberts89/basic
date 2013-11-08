#TITLE: split_long_seqs.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 5/7/13
#PURPOSE: splits contigs over length specified

#!/usr/bin/perl -w
use strict;

my $file = $ARGV[1];
open FILE, "<$file";
my $output = $ARGV[2];
open OUT, ">$output";
my $cutoff = $ARGV[0];

my %sequences;
my $contig;
my $sequence;

while (my $line = <FILE>) {                  
    chomp $line;
    if ($line =~ /^>/) {                  
        $contig = $line;
        $contig =~ s/>//;
        $contig=~ s/\|.*//;
    }
    else{             
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
}

foreach my $cont (sort(keys %sequences)) {
    if (length($sequences{$cont}) > $cutoff) {
        my $seq_length = length($sequences{$cont});
        my @abc = ('a'..'z','aa'..'zz');
        my $num_subseqs = (length($sequences{$cont}))/($cutoff-50);
        print "$num_subseqs\n";
        if (($num_subseqs - int($num_subseqs)) > 0) {
            $num_subseqs = int($num_subseqs) + 1;
        }
        my $subseq_length = length($sequences{$cont})/$num_subseqs;
        if (($subseq_length - int($subseq_length)) > 0) {
            $subseq_length = int($subseq_length) + 1;
        }
        my @sequence = split(//, $sequences{$cont});
        my $start = 0;
        print "$seq_length\t$subseq_length\t$num_subseqs\n";
        for (my $k = 0; $k <= ($num_subseqs - 1); $k++) {
            my $lexi = shift(@abc);
            my $substring_id = $cont . $lexi;
            my $end = $start + $subseq_length + 99;
            print "$k\t$start\t$end\n";
            my @subseq = @sequence[$start..$end];
            my $subseq = join('',@subseq);
            print OUT ">$substring_id\n$subseq\n";
            $start = $end - 99;
        }
        
    }
}


