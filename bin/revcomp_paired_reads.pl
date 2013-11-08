#TITLE: revcomp_paired_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 10/03/2012
#PURPOSE: takes in csv file with list of read pairs, and if

#!/usr/bin/perl -w
use strict;

my $csv = $ARGV[0];
my $fastq = $ARGV[1];
my $output = $csv;
$output =~ s/.csv/.fastq/;
my %seq_ids;

open CSV, "<$csv";
open OUT, ">$output";


#Create list of all seqIDs for the reads that are to be included in the separated fastq file
while (my $line = <CSV>) {
    next if ($line =~ /^seqID/);
    my @line = split /,/, $line;
    my $seq_id = $line[0];
    $seq_ids{$seq_id} = 1;
}

open FASTQ, "<$fastq";

#extract entry from fastq file for seqIDs that were read from the csv. Reverse complement if necessary.
while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    $line_1 =~ s/\@//;
    my $seq_id = $line_1;
    $seq_id =~ s/\/(1|2)$//;
    my $seq;
    my $plus;
    my $quality;
    my %repeats;
    while (my $line_2 = <FASTQ>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if (!exists $seq_ids{$seq_id});
    $seq = &reverse_complement($seq);
    $quality = reverse($quality);
    print OUT "$line_1\n$seq\n$plus\n$quality\n";
}
    
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}