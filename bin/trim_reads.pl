#TITLE: trim_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 12/17/2012
#PURPOSE: removes adaptor sequence from reads for mt_genome, also removes N's @ end of sequence

#!/usr/bin/perl -w
use strict;
my @input = @ARGV;


foreach my $file (@ARGV) {
    open FASTQ, "<$file";
    my $output = $file;
    $output =~ s/\.fastq/\_clipped\.fastq/;
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
        if ($seq =~ m/(GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC)/ig) {
            my @seq = split //,$seq;
            my $adapter_start = pos($seq) - length($1);
            my $adapter_stop = pos($seq);
            my @quality = split//, $quality;
            my @quality_slice = splice(@quality,$adapter_start,length($1));
            $quality_slice = join "", @quality_slice;
            $quality_slice = "\Q$quality_slice";
            $quality =~ s/$quality_slice//;
            $seq =~ s/GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC//i;
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
        print OUT "$line_1\n$seq\n$plus\n$quality\n";
    }
}   