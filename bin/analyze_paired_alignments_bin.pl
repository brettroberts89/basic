#TITLE: analyze_paired_alignments_bin.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/26/2012
#PURPOSE: looks at paired reads and determines if they are in mate pair or paired end orientation, as well as insert size and what percentage of the reads fell into the categories, as well as what percentage did not map.
#USAGE: perl analyze_paired_alignments_bin.pl alignment.sam /output_directory/

#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];

my ($seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate);
my ($mate_mapped,$mate_strand, $alt_mate);
my %mates;
my %unpaired;
open IN, "<$file";
my @file = split /\//, $file;
my $file_mod = $file[-1];
$file_mod =~ s/\.sam//;
my $read_type = $file_mod;
my $directory = $ARGV[1];
my $error=$directory . $read_type . '_error.csv';
open ERROR, ">$error";
my $output = $directory . $read_type . '_classified.csv';
open OUT, ">$output";
my $unmapped = $directory . $read_type . '_unmapped.csv';
open UNMAPPED, ">$unmapped";
my $diff_scaffolds = $directory . $read_type . '_diffscaffolds.csv';
open TWOSCAFFOLDS, ">$diff_scaffolds";
print OUT "seq_ID,Mate_pair,Paired_end,Same_strand,insert_greater_than_400bp,insert_size,read1_alt_mappings,read2_alt_mappings\n";
print UNMAPPED "seq_ID,read1_scaffold,read2__scaffold,read1_pos,read2_pos,insert,read1_mapped,read2_mapped,read1_strand,read2_strand,read1_alt_mappings,read2_alt_mappings\n";
print TWOSCAFFOLDS "seq_ID,read1_scaffold,read2__scaffold,read1_pos,read2_pos,insert,read1_mapped,read2_mapped,read1_strand,read2_strand,read1_alt_mappings,read2_alt_mappings\n";
while (my $line = <IN>) {
    next if ($line =~ /^\@SQ/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    my $read = "$seq_id,$flag,$scaffold,$mate_scaffold,$pos,$mate_pos,$insert,$seq_on_ref";
    print "Reading read: $seq_id in file: $file\n";
    my $num_alt = 0;
    $seq_id =~ s/_F3|_R3//;
    $seq_id =~ s/\/1|\/2//;
    $insert =~ s/\-//;
    if ($line =~ /\s+XA:Z:(.+)\s+/) {
        my @alternatives = split /\;/, $1;
        $num_alt = @alternatives;
    }
    my $read_mapped;
    my $read_strand;
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    $read_mapped = 1 if ($split_flag[-3] == 0);
    $read_mapped = 0 if ($split_flag[-3] == 1);
    $read_strand = "+" if ($split_flag[-5] == 0);
    $read_strand = "-" if ($split_flag[-5] == 1);
    $unpaired{$read} = "$read_mapped,$read_strand,$num_alt";
    if ($seq_id eq $seq_id_mate) {
        print "$seq_id\t$seq_id_mate\n";
        $mates{$seq_id_mate} = "$scaffold_mate,$scaffold,$pos_mate,$pos,$insert_mate,$mate_mapped,$read_mapped,$mate_strand,$read_strand,$alt_mate,$num_alt";
        my $mate_read = "$seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate";
        delete $unpaired{$mate_read};
        delete $unpaired{$read};
    }
    ($seq_id_mate,$flag_mate,$scaffold_mate,$mate_scaffold_mate,$pos_mate,$mate_pos_mate,$insert_mate,$seq_on_ref_mate) = split /,/, $read;
    $seq_id_mate =~ s/_F3|_R3//;
    $seq_id_mate =~ s/\/1|\/2//;
    $insert_mate =~ s/\-//;
    ($mate_mapped,$mate_strand, $alt_mate) = split /,/, $unpaired{$read};
}
my $same_scaffold;
my $different_scaffold;
my $total_mates_mapped;
my $one_mapped;
my $neither_mapped;
my $matepair_count;
my $pairedend_count;
my $samestrand_count;
my $biginsert_count;
my $smallinsert_count;
my $ff_count;
my $rr_count;
my $pct = $directory . $read_type . '_pct.txt';
open PCT, ">$pct";

foreach my $mate (keys %mates) {
    $total_mates_mapped += 1;
    my ($scaffold,$mate_scaffold,$pos,$mate_pos,$insert,$read_mapped,$mate_mapped,$read_strand,$mate_strand,$read_alt,$mate_alt) = split /,/, $mates{$mate};
    my $number_mapped = $read_mapped + $mate_mapped;
    if ($number_mapped == 2) {
        if ($scaffold eq $mate_scaffold) {
            $same_scaffold += 1;
            my $big_insert = "y" if ($insert > 400);
            $big_insert = "n" if ($insert <= 400);
            $biginsert_count += 1 if ($insert > 400);
            $smallinsert_count += 1 if ($insert <= 400);
            my $mate_pair;
            my $paired_end;
            my $same_strand;
            my $ff;
            my $rr;
            if ($pos < $mate_pos) {
                if (($read_strand eq "+") and ($mate_strand eq "-")) {
                    $paired_end = "y";
                    $mate_pair = "n";
                    $pairedend_count += 1;
                    $ff = "n";
                    $rr = "n";
                }
                elsif (($read_strand eq "-") and ($mate_strand eq "+")) {
                    $mate_pair = "y";
                    $paired_end = "n";
                    $matepair_count += 1;
                    $ff = "n";
                    $rr = "n";
                }
                elsif (($read_strand eq "+") and ($mate_strand eq "+")) {
                    $ff = "y";
                    $rr = "n";
                    $paired_end = "n";
                    $mate_pair = "n";
                    $ff_count += 1;
                }
                elsif (($read_strand eq "-") and ($mate_strand eq "-")) {
                    $ff = "n";
                    $rr = "y";
                    $paired_end = "n";
                    $mate_pair = "n";
                    $rr_count += 1;
                }
                print OUT "$mate,$mate_pair,$paired_end,$same_strand,$big_insert,$insert,$read_alt,$mate_alt\n";
            }
            elsif ($pos > $mate_pos) {
                if (($read_strand eq "+") and ($mate_strand eq "-")) {
                    $mate_pair = "y";
                    $paired_end = "n";
                    $matepair_count += 1;
                    $ff = "n";
                    $rr = "n";
                }
                elsif (($read_strand eq "-") and ($mate_strand eq "+")) {
                    $paired_end = "y";
                    $mate_pair = "n";
                    $pairedend_count += 1;
                    $ff = "n";
                    $rr = "n";
                }
                elsif (($read_strand eq "+") and ($mate_strand eq "+")) {
                    $paired_end = "n";
                    $mate_pair = "n";
                    $samestrand_count += 1;
                    $ff = "n";
                    $rr = "y";
                    $rr_count += 1;
                }
                elsif (($read_strand eq "-") and ($mate_strand eq "-")) {
                    $paired_end = "n";
                    $mate_pair = "n";
                    $samestrand_count += 1;
                    $ff = "y";
                    $rr = "n";
                    $ff_count += 1;
                }
                print OUT "$mate,$mate_pair,$paired_end,$same_strand,$big_insert,$insert,$read_alt,$mate_alt\n";
            }
            elsif ($pos = $mate_pos) {
                if ((($read_strand eq "+") and ($mate_strand eq "+")) or (($read_strand eq "-") and ($mate_strand eq "-"))) {
                    $same_strand = "y";
                    $paired_end = "n";
                    $mate_pair = "n";
                    $samestrand_count += 1;
                }
                else {
                    $paired_end = "y";
                    $mate_pair = "n";
                    $same_strand = "n";
                    $pairedend_count += 1;
                }
                print OUT "$mate,$mate_pair,$paired_end,$same_strand,$big_insert,$insert,$read_alt,$mate_alt\n";
            }
            else {
                print ERROR "$mate has unexpected alignment: pos:$pos mate_pos:$mate_pos\n";
            }   
        }
        else {
            $different_scaffold += 1;
            print TWOSCAFFOLDS "$mate,$mates{$mate}\n";
        }
    }
    elsif ($number_mapped == 1) {
        $one_mapped += 1;
        print UNMAPPED "$mate,$mates{$mate}\n";
    }
    elsif ($number_mapped == 0) {
        $neither_mapped += 1;
        print UNMAPPED "$mate,$mates{$mate}\n";
    }  
}
foreach my $unpaired (keys %unpaired) {
    my ($seq_id,$flag,$scaffold,$mate_scaffold,$pos,$mate_pos,$insert,$seq_on_ref) = split /,/, $unpaired;
    print ERROR "$seq_id has no mate\n";
}
if (($same_scaffold + $different_scaffold + $one_mapped + $neither_mapped) != $total_mates_mapped) {
    print ERROR "Mappings don't add up: $same_scaffold + $different_scaffold + $one_mapped + $neither_mapped <-> $total_mates_mapped\n";
}
else {
    my $pct_same_scaffold = $same_scaffold/$total_mates_mapped;
    my $pct_diff_scaffolds = $different_scaffold/$total_mates_mapped;
    my $pct_one_mapped = $one_mapped/$total_mates_mapped;
    my $pct_neither_mapped = $neither_mapped/$total_mates_mapped;
    my $pct_matepair = $matepair_count/$same_scaffold;
    my $pct_pairedend = $pairedend_count/$same_scaffold;
    my $pct_ff = $ff_count/$same_scaffold;
    my $pct_rr = $rr_count/$same_scaffold;
    my $pct_biginsert = $biginsert_count/$same_scaffold;
    my $pct_smallinsert = $smallinsert_count/$same_scaffold;
    print PCT "Percent that mapped to same scaffold: $pct_same_scaffold\nPercent that mapped to different scaffolds: $pct_diff_scaffolds\nPercent with only one mapped read: $pct_one_mapped\nPercent with no mapped reads: $pct_neither_mapped\n\nPercent in matepair orientation (RF): $pct_matepair\nPercent in paired end orientation (FR): $pct_pairedend\nPercent on same strand (FF): $pct_ff\nPercent on same strand (RR): $pct_rr\nPercent with insert > 400: $pct_biginsert\nPercent with insert size <= 400: $pct_smallinsert\n";
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}



