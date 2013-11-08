#!/usr/bin/perl -w

#TITLE: kingcsv2lddecay.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2013
#PURPOSE: reformats csv into table to show ld decay

use strict;
use warnings;

my $csv = $ARGV[0];

my ($prev_seq_id1,$prev_linkage_group,$prev_pos,$prev_scaffold1,$prev_scaffold2,$prev_scaffold1_length,$prev_scaffold2_length,$prev_pos1,$prev_pos2,$prev_mult_mapping1,$prev_mult_mapping2) = 0;
my $total_joins;
my $misjoins;
my $markers;
my $first = 1;
my $linkage_groups;
my $positions;

open CSV, "<$csv";
open OUT, ">$ARGV[1]";

while (<CSV>) {
    chomp;
    my ($seq_id1,$linkage_group,$pos,$scaffold1,$scaffold2,$scaffold1_length,$scaffold2_length,$pos1,$pos2,$mult_mapping1,$mult_mapping2) = split(/,/, $_);
    if (($scaffold1 ne $scaffold2) or ($prev_scaffold1 ne $prev_scaffold2)) {
        print "$_\n";
        next;
    }
    elsif ($first == 1) {
        $markers = $seq_id1;
        $linkage_groups = $linkage_group;
        $positions = $pos;
        $total_joins = 0;
        $misjoins = 0;
        $first = 0;
    }
    elsif (($scaffold1 eq $scaffold2) and ($prev_scaffold1 eq $prev_scaffold2) and ($scaffold1 eq $prev_scaffold2)) {
        $total_joins += 1;
        $markers .= ";$seq_id1";
        if ($prev_linkage_group ne $linkage_group) {
            $misjoins += 1;
        }
        if ($linkage_groups !~ /$linkage_group/) {
            $linkage_groups .= "/$linkage_group";
            $positions .= "/$pos";
        }
        else {
            $positions .= ";$pos"
        }
        
    }
    elsif (($scaffold1 eq $scaffold2) and ($prev_scaffold1 eq $prev_scaffold2)) {
        print OUT "$linkage_groups,$positions,$prev_scaffold1,$prev_scaffold1_length,$misjoins,$total_joins,$markers\n";
        $markers = $seq_id1;
        $total_joins = 0;
        $misjoins = 0;
        $linkage_groups = $linkage_group;
        $positions = $pos;
    }
    ($prev_seq_id1,$prev_linkage_group,$prev_pos,$prev_scaffold1,$prev_scaffold2,$prev_scaffold1_length,$prev_scaffold2_length,$prev_pos1,$prev_pos2,$prev_mult_mapping1,$prev_mult_mapping2) = split(/,/, $_);
}