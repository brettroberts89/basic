#!/usr/bin/perl -w
use strict;

my @input = @ARGV;

foreach my $file (@input) {
    print "$file\n";
    my @file = split /\//, $file;
    my $file_mod = $file[-1];
    my $read_type = $file_mod;
    my $output = "/home/dilaraally/bigdata/mtgenome/intermediate/extracted_fastq/final/$read_type";
    open OUT, ">$output";
    open FASTQ, "<$file";
    while (my $line_1 = <FASTQ>) {
        chomp $line_1;
        $line_1 = '@' . '$line_1';
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
        print OUT "$line_1\n$seq\n$plus\n$quality\n";
    }
}