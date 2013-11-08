#TITLE: parse_sam_primers.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 10/17/2012
#PURPOSE:

#!/usr/bin/perl -w
use strict;

my $sam = $ARGV[0];
my @file = split /\//, $sam;
my $file_mod = $file[-1];
$file_mod =~ s/\.sam//;
my $read_type = $file_mod;
my $output = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/linkage_map/$read_type" . '_sam.csv';
my $scaffolds = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/linkage_map/SGBv1_scaff.csv";
my $supplement = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/linkage_map/SupMat_Wangetal_2011.csv";
open OUT, ">$output";
open SAM, "<$sam";
open SUPP, "<$supplement";

my ($seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate);
my ($mate_matches,$mate_strand, $alt_mate);
my %mates;
my %unpaired;
my %scaffold_lengths;
my %supp;

open SCAFFOLDS, "<$scaffolds";

while (my $line = <SCAFFOLDS>) {
    chomp $line;
    my ($scaffold, $length) = split /,/, $line;
    $length =~ s/\D//;
    $scaffold_lengths{$scaffold} = $length;
}

while (my $line = <SUPP>) {
    next if ($line =~ /^Original_Order/);
    chomp $line;
    my ($num,$linkage_group,$marker,$rel_pos,$motif,$f_primer,$r_primer,$product_size,$genbank_id,$f_length,$r_length) = split /,/, $line;
    $supp{$marker} = "$linkage_group,$rel_pos,$motif,$product_size";
}

print OUT "Scaffold,Scaffold_length,F_primer_pos,R_primer_pos,Orientation,Linkage_group,Marker_name,Relative_pos_cM,Marker_type,Insert_size,Amplicon_size,F_primer_CIGAR,R_primer_CIGAR,Multimapping,Type\n";
while (my $line = <SAM>) {
    next if ($line =~ /^\@/);
    chomp $line;
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $read = "$seq_id,$flag,$scaffold,$mate_scaffold,$pos,$mate_pos,$insert,$seq_on_ref";
    print "Reading read: $seq_id\n";
    $insert =~ s/\*/NA/;
    $matches =~ s/\*/NA/;
    $pos =~ s/\*/NA/;
    $insert =~ s/\-//;
    my $alternatives = "NA";
    if ($line =~ /\s+XA:Z:(.+)/) {
        $alternatives = $1;
        $alternatives =~ s/\,/\:/g;
    }
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    my $read_strand = "NA";
    $read_strand = "+" if ($split_flag[-5] == 0);
    $read_strand = "-" if ($split_flag[-5] == 1);
    $unpaired{$read} = "$read_strand,$matches,$alternatives";
    if ($seq_id eq $seq_id_mate) {
        $mates{$seq_id_mate} = "$scaffold_mate,$scaffold,$pos_mate,$pos,$insert,$mate_strand,$read_strand,$mate_matches,$matches,$alt_mate,$alternatives";
        my $mate_read = "$seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate";
        delete $unpaired{$mate_read};
        delete $unpaired{$read};
    }
    ($seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate) = split /,/, $read;
    ($mate_strand, $mate_matches, $alt_mate) = split /,/, $unpaired{$read};
}

foreach my $unpaired_seq (keys %unpaired) {
    print "$unpaired_seq is unpaired\n";
}

foreach my $paired_seq (keys %mates) {
    my ($scaffold_mate,$scaffold,$pos_mate,$pos,$insert,$mate_strand,$read_strand,$mate_matches,$matches,$alt_mate,$alternatives) = split /,/, $mates{$paired_seq};
    if (($mate_matches ne "NA") and ($matches ne "NA")) {
        my $multimapping = "n";
        my $type = "NA";
        if ($scaffold_mate eq $scaffold) {
            my $orientation;
            if (($mate_strand eq "+") and ($read_strand eq "-")) {
                $orientation = "+";
            }
            elsif (($mate_strand eq "-") and ($read_strand eq "+")) {
                $orientation = "-";
            }
            elsif (($mate_strand eq "+") and ($read_strand eq "+")) {
                $orientation = "-";
            }
            elsif (($mate_strand eq "-") and ($read_strand eq "-")) {
                $orientation = "-";
            }
            else {
                $orientation = "NA";
            }
            if (($alt_mate ne "NA") and ($alternatives ne "NA")) {
                $multimapping = "y";
                my @mate_alt_scaffold;
                my @alt_scaffold;
                my @alts_mate = split /;/, $alt_mate;
                my @alts = split /;/, $alternatives;
                foreach my $alt (@alts_mate) {
                    my @array = split /:/, $alt;
                    push @mate_alt_scaffold, $array[0];
                }
                my $same_scaffold = 0;
                my @matches;
                foreach my $alt (@alts) {
                    my @array = split /:/, $alt;
                    my $scaffold_match = $array[0];
                    foreach my $mate_scaffold (@mate_alt_scaffold) {
                        if ($mate_scaffold eq $scaffold_match) {
                            push @matches, $scaffold_match;
                            $type = "samecontigforFandR=$scaffold_match";
                            $same_scaffold += 1;
                            last;
                        }
                    }
                }
                my $matches = join ";", @matches;
                if ($same_scaffold == 0) {
                    $type = "diffcontigsforFandR";
                }
                else {
                    $type = "samecontigforFandR=$matches";
                }
            }
            elsif ($alt_mate ne "NA") {
                $multimapping = "y";
                $type = "F";
            }
            elsif ($alternatives ne "NA") {
                $multimapping = "y";
                $type = "R";
            }
            my $scaffold_length = "$scaffold_lengths{$scaffold_mate}";
            my ($linkage_group,$rel_pos,$motif,$product_size) = split /,/, $supp{$paired_seq};
            print OUT "$scaffold_mate,$scaffold_length,$pos_mate,$pos,$orientation,$linkage_group,$paired_seq,$rel_pos,$motif,$insert,$product_size,$mate_matches,$matches,$multimapping,$type\n";
            print "$scaffold_mate,$scaffold_length,$pos_mate,$pos,$orientation,$linkage_group,$paired_seq,$rel_pos,$motif,$insert,$product_size,$mate_matches,$matches,$multimapping,$type\n";
        }
    }
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}