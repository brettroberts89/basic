#!/usr/bin/perl -w

#TITLE: gbs_mapping.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/25/2013
#PURPOSE: maps SNPs called in a genome to a position in another version of the genome

use strict;

my $old_genome = $ARGV[0];
my $new_genome = $ARGV[1];
my $snp_csv = $ARGV[2];

my $old_genome_file = "/shared/sgbiofuels/SGBassemblies/$old_genome/$old_genome.fa";
my $new_genome_file = "/shared/sgbiofuels/SGBassemblies/$new_genome/$new_genome.fa";
my $index = $new_genome_file;
$index =~ s/\.fa|\.fasta//;

my $before_seqs = '/bigdata/broberts/gbs/snp_mapping/before_seqs.fa';
my $after_seqs = '/bigdata/broberts/gbs/snp_mapping/after_seqs.fa';

open OLD, "<$old_genome_file";
my $contig;
my $sequence;
my %sequences;
my %snps;

while (my $line = <OLD>) {
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig = $line[0];
        $contig =~ s/>//;
        $contig =~ s/\D//g;
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequence =~ s/\W//g;
        if (exists $sequences{$contig}) {
            $sequences{$contig} .= $sequence;
        }
        else {
            $sequences{$contig} = $sequence;
        }
    }
}

close OLD;

open BEFORE, ">$before_seqs";
open AFTER, ">$after_seqs";
open SNPS, "<$snp_csv";

while (my $line = <SNPS>) {
    next if ($line =~ /^snpID/);
    chomp $line;
    my ($snp_id, $scaffold, $pos, $ref_nuc, $alt_nucs) = split(/,/, $line);
    $scaffold -= 1;
    my $before = substr($sequences{$scaffold}, ($pos - 101), 100);
    my $after = substr($sequences{$scaffold}, $pos, 100);
    print BEFORE ">$snp_id\n$before\n";
    print AFTER ">$snp_id\n$after\n";
    $snps{$snp_id} = "$ref_nuc,$alt_nucs";
}
close SNPS;
close BEFORE;
close AFTER;

my $sam = "/bigdata/broberts/gbs/snp_mapping/bwa_snps_" . $old_genome . "_" . $new_genome . ".sam";
system("bash -c 'bwa sampe $index <(bwa aln $index $before_seqs) <(bwa aln $index $after_seqs) $before_seqs $after_seqs > $sam'");

open SAM, "<$sam";
my $output = "/bigdata/broberts/gbs/snp_mapping/summary.snps.$new_genome.csv";
open OUT, ">$output";
print OUT "snpID,scaffold,bp_position,ref,alt\n";


while (my $line_1 = <SAM>) {
    my $snp_pos;
    next if ($line_1 =~ /^\@SQ/);
    next if ($line_1 =~ /^\@PG/);
    my ($seq_id1,$flag1,$scaffold1,$pos1,$mapq1,$matches1,$mate_scaffold1,$mate_pos1,$insert1,$seq_on_ref1,$qual1,$opt1) = split /\t/, $line_1;
    LINE2: while (my $line_2 = <SAM>) {
        my ($seq_id2,$flag2,$scaffold2,$pos2,$mapq2,$matches2,$mate_scaffold2,$mate_pos2,$insert2,$seq_on_ref2,$qual2,$opt2) = split /\t/, $line_2;
        if (($pos1 + 100) == ($pos2 - 1)) {
            $snp_pos = $pos2-1;
            print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
        }
        elsif (($pos1 - 1) == ($pos2 + 100)) {
            $snp_pos = $pos1-1;
            print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
        }
        else {
            print "$opt1\n$opt2\n\n";
            my $mult_mapping1;
            my $mult_mapping2;
            if ($opt1 =~ /XA:Z:(.*)$/) {
                $mult_mapping1 = $1;
            }
            if ($opt2 =~ /XA:Z:(.*)$/) {
                $mult_mapping1 = $2;
            }
            if ($mult_mapping1) {
                my @mult_mapping1 = split(/;/, $mult_mapping1);
                foreach my $mult_mapping (@mult_mapping1) {
                    my @mult_mapping = split(/,/, $mult_mapping);
                    my $scaffold = $mult_mapping[0];
                    my $pos = $mult_mapping[1];
                    $pos =~ s/[\+\-]//;
                    print "$pos2\t$pos\n";
                    if (($pos + 100) == ($pos2 - 1)) {
                        $snp_pos = $pos2-1;
                        print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
                        last LINE2;
                    }
                    elsif (($pos - 1) == ($pos2 + 100)) {
                        $snp_pos = $pos-1;
                        print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
                        last LINE2;
                    }
                }
            }
            if ($mult_mapping2) {
                my @mult_mapping2 = split(/;/, $mult_mapping2);
                foreach my $mult_mapping (@mult_mapping2) {
                    my @mult_mapping = split(/,/, $mult_mapping);
                    my $scaffold = $mult_mapping[0];
                    my $pos = $mult_mapping[1];
                    $pos =~ s/[\+\-]//;
                    if (($pos1 + 100) == ($pos - 1)) {
                        $snp_pos = $pos-1;
                        print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
                        last LINE2;
                    }
                    elsif (($pos1 - 1) == ($pos2 + 100)) {
                        $snp_pos = $pos1-1;
                        print OUT "$seq_id1,$scaffold1,$snp_pos,$snps{$seq_id1}\n";
                        last LINE2;
                    }
                }
            }
            print "ERROR: Could not identify $seq_id1 in $new_genome.\n";
        }
        last;
    }
}

