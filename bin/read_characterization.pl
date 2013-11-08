#TITLE: read_characterization.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/23/2012
#PURPOSE: 

#!/usr/bin/perl -w
use strict;
use List::Util qw(sum min);

my $sato = "Sato.fa";
my $pairs = "alnpairRun001_unmapped_F.csv";
my $mates = "alnmateRun001_unmapped_F.csv";
my $output = "read_characteristics.csv";
my $output2 = "read_characteristics_avgs.csv";
my $contig;
my %contig_lengths;
my %pairs;
my %mates;
my $diff_contigs_pair;
my $diff_contigs_mate;
my @pair_min_inserts;
my @mate_min_inserts;

open PAIRS, "<$pairs";
open MATES, "<$mates";
open OUT, ">$output";
open (SATO, "<$sato");

while (my $line = <SATO>) {                 ###create hash with keys contig name and values contig length
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/>//;
        $contig =~ s/\D+//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide
        $line =~ s/\s+//;
        my $contig_length = length($line);
        $contig_lengths{$contig} += $contig_length;
    }
}

while (my $line = <PAIRS>) {                #compile information from both lines of mates into one entry
    chomp $line;
    my ($id, $b, $insert, $contig1, $read_length1, $TF1, $coordinate1, $contig2, $TF2, $coordinate2)  = split(/,/, $line);
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

my $total_pairs = keys %pairs;
my @pair_lengths;


foreach my $id (keys %pairs) {
    my ($insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2,$read_length2) = split /,/, $pairs{$id};
    if (defined $read_length2) {
        if ($insert == 0) {
            $diff_contigs_pair += 1;
            my $length_contig1 = $contig_lengths{$contig1};
            my $length_contig2 = $contig_lengths{$contig2};
            my $a = $coordinate1;
            my $b = $length_contig1 - $read_length1 - $a;
            my $c = $coordinate2;
            my $d = $length_contig2 - $read_length2 - $c;
            my $insert1 = $a + $c;
            my $insert2 = $a + $d;
            my $insert3 = $b + $c;
            my $insert4 = $b + $d;
            my $min_insert = min $insert1, $insert2, $insert3, $insert4;
            my @contig_lengths = ($length_contig1, $length_contig2);
            push @pair_lengths, @contig_lengths;
            push @pair_min_inserts, $min_insert;
            print OUT "$id,paired_end,$length_contig1,$coordinate1,$length_contig2,$coordinate2,$min_insert\n"
        }
    }
}

while (my $line = <MATES>) {                #compile information from both lines of mates into one entry
    chomp $line;
    my ($id, $b, $insert, $contig1, $read_length1, $TF1, $coordinate1, $contig2, $TF2, $coordinate2)  = split(/,/, $line);
    $id =~ s/_subrec1//;
    next if $insert != 0;
    $contig1 =~ s/\D+//;
    $contig2 =~ s/\D+//;
    if (exists $mates{$id}) {               #if no entry exists for the pair, creates a new one
        my $read_length2 = $read_length1;
        ($insert, $contig1, $read_length1, $TF1, $coordinate1, $contig2, $TF2, $coordinate2) = split(/,/, $mates{$id});
        $mates{$id} = "$insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2,$read_length2";
    }
    else {                                  #if an entry already exists for the pair, inserts the length of the second read into the hash value
        $mates{$id} = "$insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2";
    }
}

my $total_mates = keys %mates;
my @mate_lengths;


foreach my $id (keys %mates) {
    my ($insert,$contig1,$read_length1,$TF1,$coordinate1,$contig2,$TF2,$coordinate2,$read_length2) = split /,/, $pairs{$id};
    if (defined $read_length2) {
        if ($insert == 0) {
            $diff_contigs_mate += 1;
            my $length_contig1 = $contig_lengths{$contig1};
            my $length_contig2 = $contig_lengths{$contig2};
            my $a = $coordinate1;
            my $b = $length_contig1 - $read_length1 - $a;
            my $c = $coordinate2;
            my $d = $length_contig2 - $read_length2 - $c;
            my $insert1 = $a + $c;
            my $insert2 = $a + $d;
            my $insert3 = $b + $c;
            my $insert4 = $b + $d;
            my $min_insert = min $insert1, $insert2, $insert3, $insert4;
            my @contig_lengths = ($length_contig1, $length_contig2);
            push @mate_lengths, @contig_lengths;
            push @mate_min_inserts, $min_insert;
            print OUT "$id,mate_pair,$length_contig1,$coordinate1,$length_contig2,$coordinate2,$min_insert\n"
        }
    }
}

my $pair_freq = $diff_contigs_pair/$total_pairs;
my $mate_freq = $diff_contigs_mate/$total_mates;
my $pair_avg_length = sum(@pair_lengths)/(length @pair_lengths);
my $mate_avg_length = sum(@mate_lengths)/(length @mate_lengths);
my $pair_avg_insert = sum(@pair_min_inserts)/(length @pair_min_inserts);
my $mate_avg_insert = sum(@mate_min_inserts)/(length @mate_min_inserts);

print OUT2 "Type,Diff_contig_freq,Avg_contig_length,Avg_min_insert";
print OUT2 "pair,$pair_freq,$pair_avg_length,$pair_avg_insert\nmate,$mate_freq,$mate_avg_length,$mate_avg_insert";

