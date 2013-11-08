#TITLE: split_solid.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 11/19/2012
#PURPOSE: stakes reads from two files and organizes them into two files: one with forward reads and one with reverse reads, and makes sure that the reads are paired properly

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $input2 = $ARGV[1];

my %read1;
my $read1 = $input;
$read1 =~ s/\.fastq/_final\.fastq/;
my $read2 = $input2;
$read2 =~ s/\.fastq/_final\.fastq/;


open FASTQ1, "<$input";
print "okay\n";
while (my $line_1 = <FASTQ1>) {
    print "okay\n";
    chomp $line_1;
    my $seq_id = $line_1;
    $seq_id =~ s/\s1\://;
    print "$seq_id\n";
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ1>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ1>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ1>) {
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

close FASTQ1;
open READ1, ">$read1";
open READ2, ">$read2";
open FASTQ2, "<$input2";

while (my $line_1 = <FASTQ2>) {
    chomp $line_1;
    my $seq_id = $line_1;
    $seq_id =~ s/\s2\://;
    my $seq;
    my $plus;
    my $quality;
    my %repeats;
    while (my $line_2 = <FASTQ2>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ2>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ2>) {
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
        print READ1 "$read1{$seq_id}";
        print READ2 "$line_1\n$seq\n$plus\n$quality\n";
    }
}
