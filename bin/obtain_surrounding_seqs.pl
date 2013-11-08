#TITLE: obtain_surrounding_seqs.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/04/2012
#PURPOSE: obtains sequences in Sato genome of 100 bp to left and right of subsetted SNPs

#!/usr/bin/perl -w
use strict;

my $input = "SNP_subsets_combined.csv";
my $sato = "Sato.fa";
my $output = "snp_combined_w_seqs.csv";
my $contig;
my $sequence;
my %sequences;
my @hundred = (0..99);                                        #create array of numbers 0 to 100 as counter


open OUT, ">$output";
print OUT "vcf_coord,Contig_ID,SNP_position_in_contig,100bp_left,100bp_right\n";
open SATO, "<$sato";

while (my $sato_line = <SATO>) {                   #go through Sato assembly, creating a hash entry for each contig with the sequence as the value
    chomp $sato_line;
    if ($sato_line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $sato_line;
        $contig =~ s/>//;
        $contig =~ s/\D+//;
        $contig =~ s/^0+//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $sato_line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$sato_line";
    }
}



open IN, "<$input";

while (my $line = <IN>) {
    chomp $line;
    next if ($line =~ /^\D/);
    my ($vcf_coord,$contig_id,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right,$geno_rich) = split(/,/, $line);
    my $hundred_left = $contig_coord - 101;                        #determining start position of 100 bases to the left of the SNP
    if ($hundred_left < 0) {
        $hundred_left = 0;
    }
    my $hundred_right = $contig_coord;                            #determining start position of 100 bases to the right of the SNP, works out since substring starts with coordinate 0
    my $left_seq = substr $sequences{$contig_id}, $hundred_left, 100; #create substring of the 100 bases to the left of the SNP
    my $right_seq = substr $sequences{$contig_id}, $hundred_right, 100; #create substring of the 100 bases to the right of the SNP
    print OUT "$vcf_coord,$contig_id,$contig_coord,$hundred_left,$hundred_right,$left_seq,$right_seq\n";
}