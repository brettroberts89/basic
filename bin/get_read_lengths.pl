#TITLE: get_read_length.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/12/13
#PURPOSE: create csv file with contig names and lengths from fastq file

#!/usr/bin/perl -w
use strict;

my $fastq = $ARGV[0];
my $output = $fastq;
$output =~ s/\.fastq/_hist\.csv/;
open FASTQ, "<$fastq";
open OUT, ">$output";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my $seq;
    my $plus;
    my $quality;
    my $length;
    while (my $line_2 = <FASTQ>) {
        chomp $line_2;
        $seq = $line_2;
        $length = length($seq);
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
    print OUT "$line_1,$length\n";
}