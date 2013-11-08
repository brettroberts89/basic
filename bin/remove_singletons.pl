#!/usr/bin/perl -w

#TITLE: remove_singletons.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/8/2013
#PURPOSE: removes singletons from paired fastq files


use strict;
use warnings;


my $fastq1_in = $ARGV[0];
my $fastq2_in = $ARGV[1];
my $fastq1_out = $ARGV[2];
my $fastq2_out = $ARGV[3];
my $singleton_out = $ARGV[4];
my %read1;

open FASTQ1IN, "<$fastq1_in";
while (my $line_1 = <FASTQ1IN>) {
    chomp $line_1;
    my @seq_id = split /\s+/, $line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    $seq_id =~ s/\/1//;
    $seq_id =~ s/f\.ab1//;
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ1IN>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ1IN>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ1IN>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if ($seq !~ /\D/);
    $read1{$seq_id} = "$line_1\n$seq\n$plus\n$quality\n";
}

close FASTQ1IN;
open FASTQ2IN, "<$fastq2_in";
open FASTQ1OUT, ">$fastq1_out";
open FASTQ2OUT, ">$fastq2_out";
open SINGLETONOUT, ">$singleton_out";

while (my $line_1 = <FASTQ2IN>) {
    chomp $line_1;
    my @seq_id = split /\s/,$line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    $seq_id =~ s/\/2//;
    $seq_id =~ s/r\.ab1//;
    my $seq;
    my $plus;
    my $quality;
    my %repeats;
    while (my $line_2 = <FASTQ2IN>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ2IN>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ2IN>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if ($seq !~ /\D/);
    if (exists $read1{$seq_id}) {
        print FASTQ1OUT "$read1{$seq_id}";
        print FASTQ2OUT "$line_1\n$seq\n$plus\n$quality\n";
        delete $read1{$seq_id};
    }
    else {
        print SINGLETONOUT "$line_1\n$seq\n$plus\n$quality\n";
    }
}

foreach my $key (keys %read1) {
    print SINGLETONOUT "$read1{$key}";
}