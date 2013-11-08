#TITLE: sort_fasta.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 5/6/13
#PURPOSE: sorts fasta based on contig name

#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];
open FILE, "<$file";
my $output = $ARGV[1];
open OUT, ">$output";

my %sequences;
my $contig;
my $sequence;

while (my $line = <FILE>) {                  
    chomp $line;
    if ($line =~ /^>/) {                  
        $contig = $line;
        $contig =~ s/>//;
    }
    else{             
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
}

foreach my $cont (sort(keys %sequences)) {
    print OUT ">$cont\n$sequences{$cont}\n";
}
