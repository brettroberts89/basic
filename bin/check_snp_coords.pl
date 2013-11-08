#!/usr/bin/perl -w

#TITLE: check_snp_coords.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2013
#PURPOSE: uses two sam files from the alignment of reads to two genomes to validate the coordinates of snps transferred from one genome to the other.

use strict;
use warnings;

my $sam1 = $ARGV[0];
my $sam2 = $ARGV[1];

my $coords1 = $ARGV[2];
my $coords2 =  $ARGV[3];

my %sgbv1_snps;
my %sgbv11_snps;
my %snp_reads;

open COORDS1, "<$coords1";

while (my $line = <COORDS1>) {
    next if ($line =~ /^snpID/);
    chomp $line;
    my ($snp_id, $scaffold, $position, $ref, $alt) = split(/,/, $line);
    $sgbv1_snps{$scaffold,$position} = $snp_id;
}

open COORDS2, "<$coords2";

while (my $line = <COORDS2>) {
    next if ($line =~ /^snpID/);
    chomp $line;
    my ($snp_id, $scaffold, $position, $ref, $alt) = split(/,/, $line);
    $sgbv11_snps{$scaffold,$position} = $snp_id;
}

open SAM1, "$<$sam1";

while (my $line = <SAM1>) {
    chomp $line;
    next if ($line =~ /^\@SQ/);
    next if ($line =~ /^\@PG/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$length,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $left_coord = $pos;
    my $right_coord = $pos + $length - 1;
    foreach my $key (keys %sgbv1_snps) {
        my ($snp_scaffold, $snp_pos) = split(/,/, $key);
        if ($scaffold eq $snp_scaffold) {
            if (($snp_pos <= $right_coord) and ($snp_pos >= $left_coord)) {
                if (exists $snp_reads{$sgbv1_snps{$key}}) {
                    $snp_reads{$sgbv1_snps{$key}} .= ",$seq_id";
                }
                else {
                    $snp_reads{$sgbv1_snps{$key}} .= $seq_id;
                }
            }
        }
    }
}

open SAM2, "<$sam2";

while (my $line = <SAM2>) {
    chomp $line;
    next if ($line =~ /^\@SQ/);
    next if ($line =~ /^\@PG/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$length,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $left_coord = $pos;
    my $right_coord = $pos + $length - 1;
    if (exists $sgbv1_snps{$scaffold,$pos})  {
        if (exists $snp_reads{$sgbv1_snps{$scaffold,$pos}}) {
            $snp_reads{$sgbv1_snps{$scaffold,$pos}} .= ",$seq_id";
        }
        else {
            $snp_reads{$sgbv1_snps{$scaffold,$pos}} .= $seq_id;
        }
    }
}
