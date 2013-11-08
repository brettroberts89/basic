#!/usr/bin/perl -w

#TITLE: split_fastq.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 9/12/2013
#PURPOSE: splits a single fastq with paired reads into 2 files: read1, read2. singletons are appended to the end of read2 fastq

use strict;
use warnings;

my $fastq_in = $ARGV[0];
die if ($fastq_in !~ /\.fastq/);
my $fastq1_out = $fastq_in;
my $fastq2_out = $fastq_in;
$fastq1_out =~ s/\.fastq/\_A.fastq/;
$fastq2_out =~ s/\.fastq/\_B.fastq/;
my %read1;

open FASTQ1OUT, ">$fastq1_out";
open FASTQ2OUT, ">$fastq2_out";

open FASTQIN, "<$fastq_in";
while (my $line_1 = <FASTQIN>) {
    chomp $line_1;
    my $seq_id = $line_1;
    $seq_id =~ s/\@//;
    $seq_id =~ s/\/1//;
    $seq_id =~ s/f\.ab1//;
    $seq_id =~ s/\/2//;
    $seq_id =~ s/r\.ab1//;
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQIN>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQIN>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQIN>) {
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
        $read1{$seq_id} = "$line_1\n$seq\n$plus\n$quality\n";
    }
}

foreach my $key (keys %read1) {
    print FASTQ2OUT "$read1{$key}";
}