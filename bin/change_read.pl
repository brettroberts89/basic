#TITLE: change_read.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 02/28/13
#PURPOSE: change /1 to /2 and /2 to /1 in seq_ID

#!/usr/bin/perl -w
use strict;

my $fastq = $ARGV[0];
my $output = $fastq;
$output =~ s/\.fastq/_switched\.fastq/;
open FASTQ, "<$fastq";
open OUT, ">$output";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my $seq_id = $line_1;
    if ($line_1 =~ /\/1/) {
        $line_1 =~ s/\/1/\/2/;
    }
    else {
        $line_1 =~ s/\/2/\/1/;
    }
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
    print OUT "$line_1\n$seq\n$plus\n$quality\n";
}