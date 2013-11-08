#TITLE: fastq_base_distribution.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 6/19/13
#PURPOSE: creates matrix of frequency of each base pair at first 20 positions in the the read and 20 positions prior to read.

#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);

my $fastq = $ARGV[0];
my $sam = $ARGV[1];
my $fasta = $ARGV[2];
my $output = $ARGV[3];

my $contig;
my $sequence;
my %sequences;

open FASTA, "<$fasta";

while (my $line = <FASTA>) {
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig = $line[0];
        $contig =~ s/>//;
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
close FASTA;

print "Finished reading fasta.\n";

my @bases = ('G','C','A','T');
my %basecount;
my %mapped;

open SAM, "<$sam";

while (my $line = <SAM>) {
    next if ($line =~ /^\@SQ/);
    next if ($line =~ /^\@PG/);
    my @line = split(/\t/, $line);
    my $flag = $line[1];
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    my $read_mapped;
    $read_mapped = 1 if ($split_flag[-3] == 0);
    $read_mapped = 0 if ($split_flag[-3] == 1);
    my $strand;
    $strand = '+' if ($split_flag[-5] == 0);
    $strand = '-' if ($split_flag[-5] == 1);
    
    if ($read_mapped == 1) {
        my $seq_id = $line[0];
        my $contig = $line[2];
        my $cigar = $line[5];
        $mapped{$seq_id} = 1;
        if ($strand eq '+') {
            my $start_pos = $line[3] - 1;
            my $contig_seq = substr($sequences{$contig}, ($start_pos - 20), 20);
            my @contig_seq = split(//,$contig_seq);
            my $pos = -20;
            foreach my $contig_base (@contig_seq) {
                if ($pos < 0) {
                $basecount{$pos,$contig_base} += 1;
                $pos += 1;
                }
            }
        }
        else {
            my @cigar_nums = split(/\D/, $cigar);
            my $cigar_sum = sum(@cigar_nums);
            my $start_pos = $line[3] + $cigar_sum - 2;
            my $contig_seq = substr($sequences{$contig}, ($start_pos + 1), 20);
            $contig_seq = &reverse_complement($contig_seq);
            my @contig_seq = split(//,$contig_seq);
            my $pos = -20;
            foreach my $contig_base (@contig_seq) {
                if ($pos < 0) {
                $basecount{$pos,$contig_base} += 1;
                $pos += 1;
                }
            }
        }
    }
}

print "Finished reading sam\n";

close SAM;
open FASTQ, "<$fastq";

while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @seq_id = split /\s+/, $line_1;
    my $seq_id = $seq_id[0];
    $seq_id =~ s/\@//;
    #print "Processing $seq_id in fastq\n";
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if (!exists $mapped{$seq_id});
    my @seq = split(//, $seq);
    my $pos = 0;
    
    foreach my $read_base (@seq) {
        if ($pos < 20) {
            $basecount{$pos,$read_base} += 1;
            $pos += 1;
        }
    }
}
close FASTQ;

print "Finished reading fastq\n";

open OUT, ">$output";

my @count = (-20..19);
print OUT "Position\tG\tC\tA\tT\n";
foreach my $pos (@count) {
    print OUT "$pos\t$basecount{$pos,'G'}\t$basecount{$pos,'C'}\t$basecount{$pos,'A'}\t$basecount{$pos,'T'}\n"
}


sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}
