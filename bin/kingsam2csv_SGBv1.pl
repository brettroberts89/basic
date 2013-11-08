#!/usr/bin/perl -w

#TITLE: kingsam2csv.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/31/2013
#PURPOSE: creates spreadsheet for linkage analysis after mapping king primers to assembly

use strict;
use warnings;

my $sam = $ARGV[0];
open SAM, "<$sam";
my $output = $sam;
$output =~ s/\.sam/\.csv/;

my %group;
my %length;
my %long_scaffold;
my %long_anchored;

open COORDS, "<$ARGV[1]";

while (my $line = <COORDS>) {
    chomp $line;
    $line =~ s/\c//g;
    $line =~ s/\s//g;
    my ($seq_id, $pos, $lg) = split /,/, $line;
    $group{$seq_id} = "$lg,$pos";
}

open LENGTHS, "<$ARGV[2]";

while (my $line = <LENGTHS>) {
    chomp $line;
    $line =~ s/\s//g;
    my ($scaffold, $length) = split /,/, $line;
    $length{$scaffold} = $length;
    if ($length > 100000) {
        $long_scaffold{$scaffold} = 1;
    }
}

open OUT, ">$output";
print OUT "Marker_ID,Linkage_Group,Position_cM,Scaffold1,Scaffold2,Scaffold1_length,Scaffold2_length,f_primer_coord,r_primer_coord,f_primer_alt_mapping,r_primer_alt_mapping\n";

my $unmapped_count;
my $mult_mapping;
my $diff_scaffolds;
my $total_count;
my $unique;
my $single;

while (my $line_1 = <SAM>) {
    chomp $line_1;
    next if ($line_1 =~ /^\@SQ/);
    next if ($line_1 =~ /^\@PG/);
    my $mult_mapping1 = 'NA';
    my $mult_mapping2 = 'NA';
    my ($seq_id1,$flag1,$scaffold1,$pos1,$mapq1,$matches1,$mate_scaffold1,$mate_pos1,$insert1,$seq_on_ref1,$qual1,$opt1) = split /\t/, $line_1;
    while (my $line_2 = <SAM>) {
        my $linkage_group;
        my $pos;
        if (exists $group{$seq_id1}) {
            ($linkage_group, $pos) = split(/,/, $group{$seq_id1});
        }
        else {
            print "$seq_id1\n";
            last;
        }
        $total_count += 1;
        chomp $line_2;
        my ($seq_id2,$flag2,$scaffold2,$pos2,$mapq2,$matches2,$mate_scaffold2,$mate_pos2,$insert2,$seq_on_ref2,$qual2,$opt2) = split /\t/, $line_2;
        
        my $bin_flag1 = &dec2bin($flag1);
        my @split_flag1 = split //, $bin_flag1;
        my $read1_mapped;
        $read1_mapped = 1 if ($split_flag1[-3] == 0);
        $read1_mapped = 0 if ($split_flag1[-3] == 1);
        
        my $bin_flag2 = &dec2bin($flag2);
        my @split_flag2 = split //, $bin_flag2;
        my $read2_mapped;
        $read2_mapped = 1 if ($split_flag2[-3] == 0);
        $read2_mapped = 0 if ($split_flag2[-3] == 1);
        
        my $mult_switch = 0;
        
        if ($line_1 =~ /XA:Z:(.*)$/) {
            $mult_mapping1 = $1;
            $mult_mapping1 =~ s/\,/\:/g;
            $mult_switch += 1;
        }    
        if ($line_2 =~ /XA:Z:(.*)$/) {
            $mult_mapping2 = $1;
            $mult_mapping2 =~ s/\,/\:/g;
            $mult_switch += 1;
        }
        if (($mult_switch == 0) and ($scaffold1 eq $scaffold2) and ($read1_mapped == 1 and $read2_mapped == 1)) {
            $unique += 1;
            if (exists $long_scaffold{$scaffold1}) {
                $long_anchored{$scaffold1} = 1;
            }
        }
        $mult_mapping += 1 if ($mult_switch == 1);

        $diff_scaffolds += 1 if ($scaffold1 ne $scaffold2);
        if (($read1_mapped == 1 and $read2_mapped == 0) or ($read1_mapped == 0 and $read2_mapped == 1)) {
            $single += 1;
        }
        
        my $scaffold1_length = 'NA';
        my $scaffold2_length = 'NA';
        if (exists $length{$scaffold1}) {
            $scaffold1_length = $length{$scaffold1};
        }
        if (exists $length{$scaffold2}) {
            $scaffold2_length = $length{$scaffold2};
        }
        $scaffold1="NA" if ($read1_mapped == 0);
        $scaffold2="NA" if ($read2_mapped == 0);
        $pos1 = "NA" if ($read1_mapped == 0);
        $pos2 = "NA" if ($read2_mapped == 0);
        print OUT "$seq_id1,$linkage_group,$pos,$scaffold1,$scaffold2,$scaffold1_length,$scaffold2_length,$pos1,$pos2,$mult_mapping1,$mult_mapping2\n";
        last;
    }
}

my $long_scaffold = keys %long_scaffold;
my $long_anchored = keys %long_anchored;

print "Total markers: $total_count\nUnique: $unique\nSingle: $single\n$long_anchored long scaffolds anchored out of $long_scaffold\n";

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}
