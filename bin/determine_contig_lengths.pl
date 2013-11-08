#TITLE: determine_contig_lengths.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 11/07/2012
#PURPOSE: creates .csv file with contig names and lengths

#!/usr/bin/perl -w
use strict;

my @files = @ARGV;


foreach my $file (@files) {
    my %contig_lengths;
    my $contig;
    my $sequence;

    open FASTA, "<$file";
    my @file = split /\//, $file;
    my $file_mod = $file[-1];
    $file_mod =~ s/\.fasta//;
    my $assembly = $file_mod;
    my $output = $assembly . "_lengths.csv";
    open OUT, ">$output";
    print OUT "contig,length\n";
    while (my $line = <FASTA>) {
        chomp $line;
        if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
            $contig = $line;
            $contig =~ s/>//;
        }
        else{           #adds on to length of contig while the the line begins with a nucleotide    
            $sequence = $line;
            $sequence =~ s/\s+//;
            $contig_lengths{$contig} += length($sequence);
        }
    }
    foreach my $contig_a (keys %contig_lengths) {
        print OUT "$contig_a,$contig_lengths{$contig_a}\n";
    }
}



