#!/usr/bin/perl -w

#TITLE: gbs_read_filter.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/23/2013
#PURPOSE: identifies reads with barcode, removes barcode, filters for quality, identifies and removes common adapter, and separates reads based on adapter sequence

use strict;

my %read;
my %barcodes;
my %gender;
my %indiv;
my $count;


open GENDER, "<$ARGV[2]";

while (<GENDER>) {
    chomp;
    my @line = split /,/, $_;
    $gender{$line[1]} = $line[2];
}

open KEY, "<$ARGV[1]";

while (<KEY>) {
    chomp;
    my @line = split /\t/, $_;
    if (exists $barcodes{"$line[0],$line[1]"}) {
        $barcodes{"$line[0],$line[1]"} .= ",$line[2]";
    }
    else {
        $barcodes{"$line[0],$line[1]"} = $line[2];
    }
    $indiv{$line[2]} = $line[3];
}

my $fastq_mid = $ARGV[0];
$fastq_mid =~ s/\.fastq/_qual_length_filt\.fastq/;

system("fastq_quality_trimmer -Q33 -t 20 -l 36 -i $ARGV[0] | fastq_quality_filter -Q33 -q 20 -p 70 > $fastq_mid");

my $female_reads = $fastq_mid;
my $male_reads = $fastq_mid;
$female_reads =~ s/\.fastq/_F\.fastq/;
$male_reads =~ s/\.fastq/_M\.fastq/;

open FASTQ, "<$fastq_mid";
open MALE, ">$male_reads";
open FEMALE, ">$female_reads";

LINE: while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my @seq_id = split /:|\s/, $line_1;
    my $flowcell = $seq_id[2];
    my $lane = $seq_id[3];
    my $filter = $seq_id[8];
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
    next if ($seq !~ /\D/);
    next if ($filter eq 'Y');

    my $barcodes = $barcodes{"$flowcell,$lane"};
    my @barcodes = split(/,/, $barcodes);
    @barcodes = reverse @barcodes;
    
    my $original_seq = $seq;
    #####remove common adapter and any following sequence if present
    $seq =~ s/^N+//;
    $seq =~ s/C[AT]GAGATC.*$//;
    $seq =~ s/AGATCGGA.*$//;
    
    my $barcode;
    my @seqs;
    
    
    foreach my $code (@barcodes) {
        if ($seq =~ /^\Q$code\E(C[AT]GC.*)$/) {
            $seq = $1;
            if ($seq =~ /GC[AT]GC/) {
                my $find = 'GC([AT])GC';
                my $replace = '"G,C" . $1 . "GC"';
                $seq =~ s/$find/$replace/gee;
                @seqs = split(/,/, $seq);
            }
            $barcode = $code;
            last;
        }
    }
    
    if (!$barcode) {
        foreach my $code (@barcodes) {
            my $temp_code = $code;
            $code =~ s/^.//;
            if ($seq =~ /^\Q$code\E(C[AT]GC.*)$/) {
                $seq = $1;
                if ($seq =~ /GC[AT]GC/) {
                    my $find = 'GC([AT])GC';
                    my $replace = '"G,C" . $1 . "GC"';
                    $seq =~ s/$find/$replace/gee;
                    @seqs = split(/,/, $seq);
                }
                $barcode = $temp_code;
                last;
            }
        }
    }
    
    if (!$barcode) {
        foreach my $code (@barcodes) {
            my $revcomp_code = &reverse_complement($code);
            if ($seq =~ /^(.*G)C[AT]G\Q$revcomp_code\E/) {
                $seq = $1;
                if ($seq =~ /GC[AT]GC/) {
                    my $find = 'GC([AT])GC';
                    my $replace = '"G,C" . $1 . "GC"';
                    $seq =~ s/$find/$replace/gee;
                    @seqs = split(/,/, $seq);
                }
                $barcode = $code;
                last;
            }
        }
    }
    
    my $indiv;
    
    if (defined $barcode) {
        $indiv = $indiv{$barcode};
        $line_1 .= $indiv;
    }
    else {
        #print "No barcode:\n$line_1\n$original_seq\n$plus\n$quality\n\n";
        next;
    }
    
    if (@seqs) {
        foreach my $subseq (@seqs) {
            next if (length($subseq) < 36);
            $original_seq =~ /\Q$subseq\E/;
            my $start = $-[0];
            my $qual = substr($quality, $start, length($subseq));
            if ($gender{$indiv} eq 'F') {
                print FEMALE "$line_1\n$subseq\n$plus\n$qual\n";
            }
            elsif ($gender{$indiv} eq 'M') {
                print MALE "$line_1\n$subseq\n$plus\n$qual\n";
            }
            else {
               print "Unknown gender for $indiv\n"; 
            }
        }
    }
    else {
        next if (length($seq) < 36);
        $original_seq =~ /\Q$seq\E/;
        my $start = $-[0];
        my $qual = substr($quality, $start, length($seq));
        if ($gender{$indiv} eq 'F') {
            print FEMALE "$line_1\n$seq\n$plus\n$qual\n";
        }
        elsif ($gender{$indiv} eq 'M') {
            print MALE "$line_1\n$seq\n$plus\n$qual\n";
        }
        else {
           print "Unknown gender for $indiv\n"; 
        }
    }
}

print "$count reads\n";



sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}




