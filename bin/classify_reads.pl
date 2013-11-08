#TITLE: classify_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/16/2012
#PURPOSE: takes in two files: a .fa file of the Sato Assembly and a file with mate reads that mapped to different contigs. Determines the possible lengths of the insert,
### and based on a cutoff of 250 bases, defines the pairs as mate pairs or paired end reads.

#!/usr/bin/perl -w
use strict;
use List::Util qw(max min);

my $sato = "Sato.fa";
my $reads = "mytemp.csv";
my $output = "reads_classified.csv";
my $test = "contig_lengths.txt";
my $contig;
my %contig_lengths;
my %pairs;
open(OUT, ">$output");
open(TEST, ">$test");
#my $problem1;
#my $problem2;
#open(PROBLEM, ">problem.txt");

#print PROBLEM "Contig_length,Coordinate\n";

print TEST "Contig\tLength\n";
print OUT "Read_ID,Insert1,Insert1_Mate,Insert2,Insert2_Mate,Insert3,Insert3_Mate,Insert4,Insert4_Mate\n";

open (SATO, "<$sato");

while (my $line = <SATO>) {                 ###create hash with keys contig name and values contig length
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/>//;
        $contig =~ s/\D+//;
    }
    else {           #adds on to length of contig while the the line begins with a nucleotide    
        $line =~ s/\s+//;
        my $contig_length = length($line);
        $contig_lengths{$contig} += $contig_length;
    }
}

foreach $a (keys %contig_lengths) {
    print TEST "$a\t$contig_lengths{$a}\n";
}

open (READS, "<$reads");

while (my $line = <READS>) {                #compile information from both lines of mates into one entry
    chomp $line;
    my ($id, $insert, $contig1, $read_length1, $TF1, $coordinate1, $contig2, $TF2, $coordinate2)  = split(/,/, $line);
    $id =~ s/_subrec1//;
    next if $insert != 0;
    $contig1 =~ s/\D+//;
    $contig2 =~ s/\D+//;
    if (exists $pairs{$id}) {               #if no entry exists for the pair, creates a new one
        my $read_length2 = $read_length1;
        ($insert, $contig1, $read_length1, $TF1, $coordinate1, $contig2, $TF2, $coordinate2) = split(/,/, $pairs{$id});
        $pairs{$id} = "$insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2,$read_length2";
    }
    else {                                  #if an entry already exists for the pair, inserts the length of the second read into the hash value
        $pairs{$id} = "$insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2";
    }
}

foreach my $id (keys %pairs) {              #calculate possible lengths of insert for each pair
    my ($insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2,$read_length2) = split(/,/, $pairs{$id});
    my $length_contig1 = $contig_lengths{$contig1};
    my $length_contig2 = $contig_lengths{$contig2};
#    if ($length_contig1 < $coordinate1) {
 #       print PROBLEM "$length_contig1,$coordinate1\n";
  #  }
   # if ($length_contig2 < $coordinate2) {
    #    print PROBLEM "$length_contig2,$coordinate2\n";
    #}
    my $a = $coordinate1;
    my $b = $length_contig1 - $read_length1 - $a;
    my $c = $coordinate2;
    my $d = $length_contig2 - $read_length2 - $c;
    my $insert1 = $a + $c;
    my $insert2 = $a + $d;
    my $insert3 = $b + $c;
    my $insert4 = $b + $d;
    #my $insert_max = max $insert1, $insert2, $insert3, $insert4;
    #my $insert_min = min $insert1, $insert2, $insert3, $insert4;
    
    my $mate1;                          #defines a mate pair as anything with an insert greater than 500 bp and a paired end as anything with an insert less than or equal to 500 bp
    if ($insert1 > 250) {
        $mate1 = "yes";
    }
    if ($insert1 <= 250) {
        $mate1 = "no";
    }
    my $mate2;
        if ($insert2 > 250) {
        $mate2 = "yes";
    }
    if ($insert2 <= 250) {
        $mate2 = "no";
    }
        my $mate3;                          
    if ($insert3 > 250) {
        $mate3 = "yes";
    }
    if ($insert3 <= 250) {
        $mate3 = "no";
    }
    my $mate4;
        if ($insert4 > 250) {
        $mate4 = "yes";
    }
    if ($insert4 <= 250) {
        $mate4 = "no";
    }
    print OUT "$id,$insert1,$mate1,$insert2,$mate2,$insert3,$mate3,$insert4,$mate4\n"
}



