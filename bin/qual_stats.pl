#TITLE: qual_stats.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 1/29/13
#PURPOSE: determines number of entries in .qual file, avg length, and avg quality score

#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
my $qual = $ARGV[0];

open QUAL, "<$qual";

my $contig;
my $sequence;
my %sequences;

while (my $line = <QUAL>) {
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/^>//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
}

my $number_reads;
my $length_sum;
my $quality_sum;

foreach my $read (keys %sequences) {
    $number_reads += 1;
    my $qualities = $sequences{$read};
    my @qualities=split /\s+/, $qualities;
    my $read_length = @qualities;
    my $average_quality = sum(@qualities)/scalar(@qualities);
    $length_sum += $read_length;
    $quality_sum += $average_quality;
}

my $avg_quality = $quality_sum/$number_reads;
my $avg_length = $length_sum/$number_reads;

print "$qual\nNumber of Reads: $number_reads\nAvg Length: $avg_length\nAvg Quality: $avg_quality\n\n";

