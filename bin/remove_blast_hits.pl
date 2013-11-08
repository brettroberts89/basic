#TITLE: remove_blast_hits.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 06/26/13
#PURPOSE: removes sequences that had a good blast hit from a fasta file

#!/usr/bin/perl -w
use strict;

my $blast = $ARGV[0];
my $fasta = $ARGV[1];
my $output = $ARGV[2];

my %dump;
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

open BLAST, "<$blast";
while (my $line = <BLAST>) {
    chomp $line;
    my @line = split(/\t/, $line);
    my $seq_id = $line[0];
    $dump{$seq_id} = 1;
}

open OUT, ">$output";

foreach my $sample (keys %sequences) {
    next if (exists $dump{$sample});
    print OUT ">$sample\n$sequences{$sample}\n";
}
