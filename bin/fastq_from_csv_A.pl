#TITLE: revcomp_paired_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 10/03/2012
#PURPOSE: takes in csv file with list of read pairs, and if

#!/usr/bin/perl -w
use strict;

my $sam = $ARGV[0];
my $fastq = $ARGV[1];
my $output = $sam;
$output =~ s/.sam/_A_.fastq/;
my %seq_ids;

open SAM, "<$sam";
open OUT, ">$output";


#Create list of all seqIDs for the reads that are to be included in the separated fastq file
while (my $line = <SAM>) {
    next if ($line =~ /^\@/);
    my ($seq_id,$flag,$scaffold,$pos,$mapq,$matches,$mate_scaffold,$mate_pos,$insert,$seq_on_ref,$qual,$opt) = split /\t/, $line;
    print "Reading read: $seq_id in file: $sam\n";
    $insert =~ s/\-//;
    my $read_mapped;
    my $bin_flag = &dec2bin($flag);
    my @split_flag = split //, $bin_flag;
    $read_mapped = 1 if ($split_flag[-3] == 0);
    $read_mapped = 0 if ($split_flag[-3] == 1);
    $seq_ids{$seq_id} = 1 if ($read_mapped == 1);
}

open FASTQ, "<$fastq";

#extract entry from fastq file for seqIDs that were read from the csv. Reverse complement if necessary.
while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    $line_1 =~ s/\@//;
    my $seq_id = $line_1;
    $seq_id =~ s/\/(1|2)$//;
    my $seq;
    my $plus;
    my $quality;
    my %repeats;
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
    next if (!exists $seq_ids{$seq_id});
    print "Reading $seq_id\n";
    print OUT "$line_1\n$seq\n$plus\n$quality\n";
}
    
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}