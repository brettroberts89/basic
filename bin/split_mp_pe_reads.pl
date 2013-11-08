#TITLE: split_mp_pe_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/05/2013
#PURPOSE: takes in csv file and fastq files. outputs fastq files of matepair reads, paired end reads, and unmapped reads. Determines type of reads based on insert size. Reverse complements reads to get correct orientation if necessary.

my $csv = $ARGV[0];
my $fastq1 = $ARGV[1];
my $fastq2 = $ARGV[2];

open CSV, "<$csv";

my (%matepair,%mp_rev,%pairedend,%pe_rev);

while (my $line = <CSV>) {
    chomp $line;
    my @line = split(/,/,$line);
    my $seq_id = $line[0];
    my $matepair = $line[1];
    my $pairedend = $line[2];
    my $big = $line[4];
    if ($matepair == "y") {
        if ($big == "y") {
            $matepair{$seq_id} = 1;
        }
        else {
            $pe_rev{$seq_id} = 1;
        }
    }
    elsif ($pairedend == "y") {
        if ($big == "n") {
            $pairedend{$seq_id} = 1;
        }
        else {
            $mp_rev{$seq_id} = 1;
        }
    }
}

while (my $line_1 = <FASTQ1>) {
    chomp $line_1;
    my $seq_id = $line_1;
    $seq_id =~ s/\@//;
    $seq_id =~ s/\/1//;
    $seq_id =~ s/\/2//;
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ1>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ1>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ1>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    if (exists $mp_rev{$seq_id}) {
        $seq = &reverse_complement($seq);
        $quality = reverse($quality);
        print MP "$line_1\n$seq\n$plus\n$quality\n";
    }
    elsif (exists $matepair{$seq_id}) {
        print MP "$line_1\n$seq\n$plus\n$quality\n";
    }
    elsif (exists $pe_rev{$seq_id}) {
        $seq = &reverse_complement($seq);
        $quality = reverse($quality);
        print PE "$line_1\n$seq\n$plus\n$quality\n";
    }
    elsif (exists $pairedend{$seq_id}) {
        print PE "$line_1\n$seq\n$plus\n$quality\n";
    }
    else {
        print UNMAPPED "$line_1\n$seq\n$plus\n$quality\n";
    }
}


sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}