#TITLE: split_pgm.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 11/19/2012
#PURPOSE: splits file with both forward and reverse reads into two files: one with forward reads and one with reverse reads, and makes sure that the reads are paired properly

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my @input = split /\//, $input;
my $file_end = $input[-1];
$file_end =~ s/\..+$//;
my $read_type = $file_end;

my %read1;
my %read2;
my $read1 = $file_end . "_A.fastq";
my $read2 = $file_end . "_B.fastq";
my ($seq_idprev,$line_1prev,$line_2prev,$line_3prev,$line_4prev);


open FASTQ, "<$input";
open READ1, ">$read1";
open READ2, ">$read2";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my $seq_id = $line_1;
    $seq_id =~ s/\@//;
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
    if ($line_1 =~ /\/1/) {
        $seq_idprev = $seq_id;
        $line_1prev = $line_1;
        $line_2prev = $seq;
        $line_3prev = $plus;
        $line_4prev = $quality;
    }
    if ($line_1 =~ /\/2/) {
        if ($seq_id eq $seq_idprev) {
            print READ1 "$line_1prev\n$line_2prev\n$line_3prev\n$line_4prev\n";
            print READ2 "$line_1\n$seq\n$plus\n$quality\n";
        }
        else {
            print "Error: $line_1 was not paired\n";
        }
    }
}