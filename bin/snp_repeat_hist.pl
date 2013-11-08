#TITLE: snp_repeat_hist.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2012
#PURPOSE: creates table of data to load into R to create a histogram indicating repeats in region surrounding SNPs, and a second table showing the percentage of the surrounding sequences that are unique and parts of reptitive regions.

#!/usr/bin/perl -w
use strict;

my $input = "snp_surrounding_100seq.csv";
my $output = "snp_hist.csv";
my $output2 = "snp_repeat_pct_hist.csv";
my (%monos, %dis, %tris);
my $uniques;
my %repeats;
my $multis;

open IN, "<$input";
open OUT, ">$output";
open HIST2, ">$output2";

print OUT "Category,Number,Percentage\n";
print HIST2 "Contig,Contig_Coord,Left_repeat_pct,Right_repeat_pct\n";

while (my $line = <IN>) {                       #create hashes for each type of repeats with one entry for each snp in order to delete redundancies
    if ($line =~ /unique/) {                        
    $uniques += 1;
    }
    elsif ($line !~ /unique/) {
        #print "$line\n";
        chomp $line;
        my @line = split /,/, $line;
        my $contig = $line[0];
        my $coord = $line[1];
        my $left_repeat_pct = $line[2];
        my $right_repeat_pct = $line[3];
        $repeats{"$contig,$coord,$left_repeat_pct,$right_repeat_pct"} = 1;
        if ($line =~ /mono/) {
            $monos{"$contig,$coord,$left_repeat_pct,$right_repeat_pct"} = 1;
        }
        if ($line =~ /di/) {
            $dis{"$contig,$coord,$left_repeat_pct,$right_repeat_pct"} = 1;
        }
        if ($line =~ /tri/) {
            $tris{"$contig,$coord,$left_repeat_pct,$right_repeat_pct"} = 1;
        }
    }
}
foreach my $key (keys %repeats) {                                                   #create hash entries when there is more than one type of repeat present for a SNP
    if ((exists $monos{$key}) and (exists $dis{$key}) and (exists $tris{$key})) {
        $multis += 1;
        #print "$key\n";
        delete $monos{$key};
        delete $dis{$key};
        delete $tris{$key};
    }
    elsif ((exists $monos{$key}) and (exists $dis{$key})) {
        $multis += 1;
        delete $monos{$key};
        delete $dis{$key};
    }
    elsif ((exists $tris{$key}) and (exists $dis{$key})) {
        $multis += 1;
        delete $tris{$key};
        delete $dis{$key};
    }
    elsif ((exists $monos{$key}) and (exists $tris{$key})) {
        $multis += 1;
        delete $monos{$key};
        delete $tris{$key};
    }
    
}                                                   #create variables with value of the number of hash elements for each repeat type
my $repeats = keys %repeats;
my $monos = keys %monos;
my $dis = keys %dis;
my $tris = keys %tris;
my $total = $repeats + $uniques;

my $uniques_pct = $uniques/$total;                  #calculate percentages from numbers just calculated
my $repeats_pct = $repeats/$total;
my $monos_pct = $monos/$total;
my $dis_pct = $dis/$total;
my $tris_pct = $tris/$total;
my $multis_pct = $multis/$total;

#create 2 files corresponding to the 2 different plots to be made
print OUT "uniques,$uniques,$uniques_pct\nmultis,$multis,$multis_pct\nmonos,$monos,$monos_pct\ndis,$dis,$dis_pct\ntris,$tris,$tris_pct\n";
foreach my $snp (keys %repeats) {
    print HIST2 "$snp\n";
}



