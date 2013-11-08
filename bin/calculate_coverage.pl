#TITLE: calculate_coverage.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 2/4/2013
#PURPOSE: calculates the coverage of each scaffold/contig from sam files.

#!/usr/bin/perl -w
use strict;

my @files = @ARGV;
my $length_file = shift @files;
my @sam_files = @files;
my $output = "/bigdata/broberts/mip_scaffolder/SGBv1_coverage.txt";
my %total_coverage;
my %coverage;
my %length;

foreach my $sam (@sam_files) {
    open SAM, "<$sam";
    while (my $line = <SAM>) {
        next if ($line =~ /^\@SQ/);
        my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
        next if ($scaffold eq "*");
        $total_coverage{$scaffold} += length($seq_on_ref);
    }
    close SAM;
}

open LENGTHFILE, "<$length_file";
open OUT, ">$output";
while (my $line = <LENGTHFILE>) {
    my ($scaffold,$length) = split /,/, $line;
    if (exists $total_coverage{$scaffold}) {
        $coverage{$scaffold} = $total_coverage{$scaffold}/$length;
    }
    else {
        $coverage{$scaffold} = 0;
    }
    $length{$scaffold} = $length;
}

my $count;

for ($count=1; $count<11906; $count++) {
    my $scaffold = "scaffold_$count";
    print OUT "$count\t$scaffold\t$length{$scaffold}\t$coverage{$scaffold}"
}