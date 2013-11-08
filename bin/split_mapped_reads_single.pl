#TITLE: split_mapped_reads_single.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 05/08/2012
#PURPOSE: splits reads into categories based on whether or not they mapped to the genome. Outputs file with list of sequence id's. NOT FOR PAIRED END, JUST SINGLE END READS.

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];

my %mapped;
my %unmapped;
open IN, "<$input";
while (my $line = <IN>) {
    next if ($line =~ /^\@SQ/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $read = $seq_id;
    print "Reading read: $seq_id in file: $input\n";
    $seq_id =~ s/\s1\:|\s2\://;
    $seq_id =~ s/_F3|_R3//;
    $seq_id =~ s/\/1|\/2//;
    $insert =~ s/\-//;
    my $read_mapped;
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    $mapped{$seq_id} = 1 if ($split_flag[-3] == 0);
    $unmapped{$seq_id} = 1 if ($split_flag[-3] == 1);
}

#here begins the section on going through the fastq and extracting entries and writing them somewhere else

my $fastq = $ARGV[1];

my %read;
my $read = $fastq;
my $read_mapped = $read;
$read_mapped =~ s/\.fastq/_mapped\.fastq/;
my $read_unmapped = $read;
$read_unmapped =~ s/\.fastq/_unmapped\.fastq/;
open READMAPPED, ">$read_mapped";
open READUNMAPPED, ">$read_unmapped";

open FASTQ, "<$fastq";
while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @seq_id = split /\s+/, $line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    $seq_id =~ s/\/1//;
    print "Processing $seq_id in fastq\n";
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
    next if ($seq !~ /\D/);
        if (exists $mapped{$seq_id}) {
        print READMAPPED "$line_1\n$seq\n$plus\n$quality\n";
    }
    if (exists $unmapped{$seq_id}) {
        print READUNMAPPED "$line_1\n$seq\n$plus\n$quality\n";
    }
    else {
        print "ERROR: $seq_id not found in sam";
    }
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}