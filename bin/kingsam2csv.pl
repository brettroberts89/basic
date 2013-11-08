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

open COORDS, "<$ARGV[1]";

while (my $line = <COORDS>) {
    chomp $line;
    $line =~ s/\c//g;
    $line =~ s/\s//g;
    my ($seq_id, $pos, $lg) = split /,/, $line;
    $group{$seq_id} = "$lg,$pos";
}

open OUT, ">$output";
print OUT "Marker_ID,Linkage_Group,Position_cM,Scaffold1,Scaffold2,Scaffold1_length,Scaffold2_length,f_primer_coord,r_primer_coord,f_primer_alt_mapping,r_primer_alt_mapping\n";

while (my $line_1 = <SAM>) {
    chomp $line_1;
    next if ($line_1 =~ /^\@SQ/);
    next if ($line_1 =~ /^\@PG/);
    my $mult_mapping1 = 'NA';
    my $mult_mapping2 = 'NA';
    my ($seq_id1,$flag1,$scaffold1,$pos1,$mapq1,$matches1,$mate_scaffold1,$mate_pos1,$insert1,$seq_on_ref1,$qual1,$opt1) = split /\t/, $line_1;
    while (my $line_2 = <SAM>) {
        chomp $line_2;
        my ($seq_id2,$flag2,$scaffold2,$pos2,$mapq2,$matches2,$mate_scaffold2,$mate_pos2,$insert2,$seq_on_ref2,$qual2,$opt2) = split /\t/, $line_2;
        if ($line_1 =~ /XA:Z:(.*)$/) {
            $mult_mapping1 = $1;
            $mult_mapping1 =~ s/\,/\:/g;
        }    
        if ($line_2 =~ /XA:Z:(.*)$/) {
            $mult_mapping2 = $1;
            $mult_mapping2 =~ s/\,/\:/g;
        }
        my $linkage_group;
        my $pos;
        if (exists $group{$seq_id1}) {
            ($linkage_group, $pos) = split(/,/, $group{$seq_id1});
        }
        else {
            print "$seq_id1\n";
        }
        my $scaffold1_length = 'NA';
        my $scaffold2_length = 'NA';
        $scaffold1 =~ /(.+)\|size(.+)/;
        $scaffold1 = $1;
        $scaffold1_length = $2;
        $scaffold2 =~ /(.+)\|size(.+)/;
        $scaffold2 = $1;
        $scaffold2_length = $2;
        if (!$scaffold1_length) {
            $scaffold1_length = 'NA';
        }
        if (!$scaffold2_length) {
            $scaffold2_length = 'NA';
        }
        print OUT "$seq_id1,$linkage_group,$pos,$scaffold1,$scaffold2,$scaffold1_length,$scaffold2_length,$pos1,$pos2,$mult_mapping1,$mult_mapping2\n";
        last;
    }
}
        