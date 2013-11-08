#TITLE: reverse_complement_fastq.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 02/27/13
#PURPOSE: reverse complements the sequences in a fastq and reverses the quality scores.

#!/usr/bin/perl -w
use strict;

my $fastq = $ARGV[0];
my $output = $fastq;
$output =~ s/\.fastq/_revcomp\.fastq/;
open FASTQ, "<$fastq";
open OUT, ">$output";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my $seq_id = $line_1;
    my $seq;
    my $plus;
    my $quality;
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