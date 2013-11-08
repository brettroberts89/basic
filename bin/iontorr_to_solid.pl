#TITLE: iontorr_to_solid.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 01/30/2013
#PURPOSE: converts fastq in ion torrent format to read orientation and naming scheme of solid reads

#!/usr/bin/perl -w
use strict;

my $fastq = $ARGV[0];
my $fastq_out = $fastq;
$fastq_out =~ s/\.fastq/_mipready\.fastq/;
open OUT, ">$fastq_out";

open FASTQ, "<$fastq";
while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @seq_id = split /\s+/, $line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    print "$seq_id\n";
    $line_1 =~ s/\/1/\_R3/;
    $line_1 =~ s/\/2/\_F3/;
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
    next if ($seq !~ /\D/);
    if ($line_1 =~ /\/2/) {
        $quality = reverse $quality;
        $seq = &reverse_complement($seq);
    }
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