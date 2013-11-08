#TITLE: remove_offset_hmp.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/16/2012
#PURPOSE: removes offset of 100 bases between each contig and assigns new coordinates, then uses information in mergedtbt.c1.hmp to determine which contig the coordinates correspond to and the contig coordinate.

#!/usr/bin/perl -w
use strict;

my $input = "jatropha_ApeKI_ref.txt";
my $hmp = "mergedtbt.c1.hmp";
my $output = "mergedtbt.c1_contig_coords.csv";
my $offset = 0;
my %contigs;
my $prev_end = 0;


open OUT, ">$output";
open IN, "<$input";

#print OUT "Contig,Contig_Length,Start_Coord_with_Offset,End_Coord_with_Offset,Start_Coord,End_Coord\n";
print OUT "vcf_pos,contig,contig_pos,contig_length\n";
while (my $line = <IN>) {                                                       #from index file, create hash with keys contigs and values of the start and stop positions with the offset
    chomp $line;
    my ($contig, $start_offset, $end_offset) = split(/\D+/, $line);
    my $length = $end_offset - $start_offset;
    $contigs{$contig} = "$start_offset\t$end_offset\n";
    #print OUT "$contig,$start_no_offset,$end_no_offset\n";
}

open HMP, "<$hmp";
while (my $line = <HMP>) {
    chomp $line;
    if ($line =~ /^S1/) {                        #if the line begins with a 1, which is the case for everything but the header
        my @line = split(/\t/, $line);          #split the line into an array at any value that is not a digit
        my $coordinate = $line[3];              #the coordinate is the second value in the array
        my $contig = "unknown";                 #sets values of contig, contig_coord, and contig_length as unknown if the coordinate does fall within any of the contig regions
        my $contig_coord = "unknown";
        my $contig_length = 'unknown';
        foreach my $a (keys %contigs) {                     
            my ($start_coord, $end_coord) = split (/\t/, $contigs{$a});
            if ($coordinate >= $start_coord and $coordinate <= $end_coord) {    #goes through all contigs and checks if the coordinate is within the region of the contig    
                $contig = $a;                                                   #variable $contig is set as the key
                $contig_coord = $coordinate - $start_coord + 1;                 #the coordinate relative to the beginning of the contig is calculated by subtracting the start coordinate of the contig from the coordinate of the snp, then adding 1 to make the first base in the contig have a coordinate of 1   
                $contig_length = $end_coord - $start_coord;                     #the length of the contig is calculated by subtracting the start coordinate of the contig from the end coordinate of the contig
                last;                                                           #stops searching since the contig has been determined
            }
            else { 
                next;
            }
        }
        print OUT "$coordinate,$contig,$contig_coord,$contig_length\n";
    }
}