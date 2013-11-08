#TITLE: length_filter.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/26/2013
#PURPOSE: removes reads below a certain length

#!/usr/bin/perl -w
use strict;
my @input = @ARGV;


foreach my $file (@ARGV) {
    open FASTQ, "<$file";
    my $output = $file;
    $output =~ s/\.fastq/\_lenfiltered\.fastq/;
    open OUT, ">$output";
    while (my $line_1 = <FASTQ>) {
        chomp $line_1;
        my $seq;
        my $plus;
        my $quality;
        my %repeats;
        my $quality_slice;
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
        if ($seq =~ m/(N+$)/ig) {
            my @seq = split //,$seq;
            my $n_start = pos($seq) - length($1);
            my $n_stop = pos($seq);
            my @quality = split//, $quality;
            my @quality_slice = splice(@quality,$n_start,length($1));
            $quality_slice = join "", @quality_slice;
            $quality_slice = "\Q$quality_slice";
            $quality =~ s/$quality_slice//;
            $seq =~ s/N+$//i;
        }
        if (length($seq) >= 150) {
            print OUT "$line_1\n$seq\n$plus\n$quality\n";
        }
    }
}   