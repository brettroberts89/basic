#TITLE: split_mapped_reads_alt.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 01/16/2012
#PURPOSE: splits reads into categories based on whether or not they mapped to the genome. Outputs file with list of sequence id's

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];

my %mapped;
my %unmapped;

my ($seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate);
my ($mate_mapped,$mate_strand, $alt_mate);
my %mates;
my %unpaired;
open IN, "<$input";
my @file = split /\//, $input;
my $file_mod = $file[-1];
$file_mod =~ s/\.sam//;
my $read_type = $file_mod;
my $prev_seq;
while (my $line = <IN>) {
    next if ($line =~ /^\@SQ/);
    next if ($line =~ /^\@PG/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $read = $seq_id;
    #print "Reading read: $seq_id in file: $input\n";
    $seq_id =~ s/\s1\:|\s2\://;
    $seq_id =~ s/_F3|_R3//;
    $seq_id =~ s/\/1|\/2//;
    $insert =~ s/\-//;
    my $read_mapped;
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    $read_mapped = 1 if ($split_flag[-3] == 0);
    $read_mapped = 0 if ($split_flag[-3] == 1);
    if ($read_mapped == 1) {
        $read = $read . '/1' if ($read ne $prev_seq);
        $read = $read . '/2' if ($read eq $prev_seq);
        $mapped{$read} = 1;
    }
    else {
        $unmapped{$read} = 1;
    }
    $prev_seq = $read;
}

#here begins the section on going through the fastq and extracting entries and writing them somewhere else

my $fastq1 = $ARGV[1];
my $fastq2 = $ARGV[2];

my %read1;
my $read1 = $fastq1;
my $read_mapped = $read1;
$read_mapped =~ s/\.fastq/_mapped\.fastq/;
my $read_unmapped = $read1;
$read_unmapped =~ s/\.fastq/_unmapped\.fastq/;
open READMAPPED, ">$read_mapped";

open FASTQ1, "<$fastq1";
while (my $line_1 = <FASTQ1>) {
    chomp $line_1;
    my @seq_id = split /\s+/, $line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    #print "Processing $seq_id in fastq\n";
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
    if (exists $mapped{$seq_id}) {
        print READMAPPED "$line_1\n$seq\n$plus\n$quality\n";
    }
}

close FASTQ1;
open FASTQ2, "<$fastq2";

while (my $line_1 = <FASTQ2>) {
    chomp $line_1;
    my @seq_id = split /\s/,$line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
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
    if (exists $mapped{$seq_id}) {
        print READMAPPED "$line_1\n$seq\n$plus\n$quality\n";
    }
}


sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}



