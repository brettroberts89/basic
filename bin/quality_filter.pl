#TITLE: quality_filter.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 12/17/2012
#PURPOSE: filters 454 reads with length <400 bp and avg quality <20.

#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
my $fasta = $ARGV[0];
my $qual = $ARGV[1];

open QUAL, "<$qual";
open FASTA, "<$fasta";

my $contig;
my $sequence;
my %sequences;
my %qualities;
my %bad_reads;

while (my $line = <QUAL>) {
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/^>//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $line;
        if (exists $qualities{$contig}) {
            $qualities{$contig} = "$qualities{$contig} " . "$line";
        }
        else {
            $qualities{$contig} = $line;
        }
    }
}

while (my $line = <FASTA>) {
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/^>//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
}

foreach my $read (keys %qualities) {
    print "okay\n";
    my $qualities = $qualities{$read};
    my @qualities=split /\s+/, $qualities;
    my $length = scalar(@qualities);
    print "$read\n$length\n$qualities\n";
    if (scalar(@qualities) < 400) {
        $bad_reads{$read} = scalar(@qualities);
        next;
    }
    my $average_quality = sum(@qualities)/scalar(@qualities);
    if ($average_quality < 20) {
        $bad_reads{$read} = $average_quality;
    } 
}

my $output_fasta = $fasta;
my $output_qual = $qual;

$output_fasta =~ s/\.fasta/\_qualfilt\.fasta/;
$output_qual =~  s/\.fasta\.qual/\_qualfilt\.fasta\.qual/;
open OUT_FASTA, ">$output_fasta";
open OUT_QUAL, ">$output_qual";
foreach my $seq (keys %sequences) {
    if (!exists $bad_reads{$seq}) {
        print OUT_FASTA ">$seq\n$sequences{$seq}\n";
        print OUT_QUAL ">$seq\n$qualities{$seq}\n";
    }
}
    