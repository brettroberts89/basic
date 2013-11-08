#TITLE: subset_SSRs.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/10/2012
#PURPOSE: obtains sequences of 150 bp to left and right of SSRs on the SGB assembly. Also removes SSRs previously tested, obtained from the Sato assembly, by using the alignment of HRM primers for these SSRs to the SGB assembly. Also selects for SSRs with GC content in teh 150 bp surrounding >.35 and dinucleotide repeates > 10 or trinucleotide repeats > 8.

#!/usr/bin/perl -w
use strict;

my $primer_aln = "aln_hrmprimers_SGB.sam";
my $sgb_ssrs_di = "SGB_results_di.csv";
my $sgb_ssrs_tri = "SGB_results_tri.csv";
my $sgb_fa = "SGB.v1.fasta";
my $same = "SSRs_in_both_assemblies.csv";                               #file containing SSRs that were already tested, and were found between two primers on the SGB assembly
my $output = "SSR_subset_with_150bp.csv";
my $mapped = "SSR_coords.csv";
my %tested_ssrs;
my %remove_ssr;

open PRIMERS, "<$primer_aln";
open MAPPED, ">$mapped";

print MAPPED "SSR_ID,Scaffold,Start_pos,End_pos\n";

while (my $line = <PRIMERS>) {                                                  #create hash of SSRs that were already tested, retaining information on the pattern and amplified region in the SGB assembly
    if ($line =~ /^S/) {
        chomp $line;
        my @line = split /\t/, $line;
        if (($line[6] == "=") and ($line [8] != 0)) {
            next if ($line[7] > $line[3]);
            $line[5] =~ s/M//;
            my $right_coord = $line[3] + $line[5] - 1;
            my ($ssr_num, $pattern) = split /\-/, $line[0];
            print MAPPED "$ssr_num,$line[2],$line[7],$right_coord\n";
            $tested_ssrs{$ssr_num} = "$pattern,$line[2],$line[7],$right_coord";
        }
    }
}


open SGB_SSRS_DI, "<$sgb_ssrs_di";
open SGB_SSRS_TRI, "<$sgb_ssrs_tri";
open SAME, ">$same";
print SAME "SSR_ID,Scaffold,Repeat_length,Start_pos,End_pos\n";

while (my $line = <SGB_SSRS_DI>) {
    next if ($line =~ /^\D/);
    chomp $line;
    my ($row,$scaffold_num,$scaffold_length,$sgb_pattern,$repeat_length,$start_pos,$end_pos,$repeat_char) = split /,/, $line;
    foreach my $tested_ssr_num (keys %tested_ssrs) {
        my ($primer_pattern, $primer_scaffold, $primer_left, $primer_right) = split /,/, $tested_ssrs{$tested_ssr_num};
        my $rev_pattern = reverse $primer_pattern;
        my $revcomp_pattern = &reverse_complement($primer_pattern);
        my $comp_pattern = reverse $revcomp_pattern;
        if (($primer_scaffold eq $scaffold_num) and ($start_pos >= $primer_left) and ($end_pos <= $primer_right) and (($sgb_pattern eq $primer_pattern) or ($sgb_pattern eq $rev_pattern) or ($sgb_pattern eq $revcomp_pattern) or ($sgb_pattern eq $comp_pattern))) {          #call two SSRs the same if it falls between the left and right coordinates of the primers and the pattern is the same;
            print SAME "$tested_ssr_num,$primer_scaffold,$start_pos,$end_pos,$sgb_pattern\n";
            $remove_ssr{"$scaffold_num,$scaffold_length,$sgb_pattern,$repeat_length,$start_pos,$end_pos"} = 2;
            last;
        }
        else {
            $remove_ssr{"$scaffold_num,$scaffold_length,$sgb_pattern,$repeat_length,$start_pos,$end_pos"} = 1;
        }
    }
}


open SGB, "<$sgb_fa";
open OUT, ">$output";

print OUT "Scaffold,Scaffold_length,Pattern,Repeat_length,Start_pos,End_pos,GC_left,GC_right,Sequence,Sequence_length\n";

my $contig;
my $sequence;
my %sequences;
my %final;

while (my $line = <SGB>) {                   #go through SGB assembly, creating a hash entry for each contig with the sequence as the value
    chomp $line;
    if ($line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $line;
        $contig =~ s/>//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$line";
    }
}

foreach my $sgb_ssr_di (keys %remove_ssr) {                                        #at this point, have removed SSRs that were already tested
    if ($remove_ssr{$sgb_ssr_di} == 1) {
        my @left_snp;
        my @right_snp;
        my $gc_right;
        my $gc_left;
        my ($scaffold_num,$scaffold_length,$pattern,$repeat_length,$start_pos,$end_pos) = split /,/, $sgb_ssr_di;
        next if ($repeat_length < 10);                              #subset for repeat length >= 10
        
        my $hundred_left = $start_pos - 150;                        #determining start position of 150 bases to the left of the SNP
        if ($hundred_left < 0) {
            $hundred_left = 0;
        }
        my $ssr_start = $start_pos;
        $start_pos += 1;
        my $hundred_right = $end_pos + 2;                            #determining start position of 100 bases to the right of the SNP, works out since substring starts with coordinate 0
        $end_pos += 2;
        my $left_seq = substr $sequences{$scaffold_num}, $hundred_left, 150; #create substring of the 100 bases to the left of the SNP
        my $right_seq = substr $sequences{$scaffold_num}, $hundred_right, 150; #create substring of the 100 bases to the right of the SNP
        my $ssr = substr $sequences{$scaffold_num}, $ssr_start, ($repeat_length*2);
        my @left_seq = split(//, $left_seq);
        my @right_seq = split(//, $right_seq);
        my $full_seq = $left_seq . $ssr . $right_seq;
        my $seq_length = length $full_seq;
        my $length_left = @left_seq;
        my $length_right = @right_seq;
        next if ($length_left < 50);
        next if ($length_right < 50);
        foreach my $base (@left_seq) {                              #calculate gc content of surrounding region
            if ($base =~ /G|C/i) {
                $gc_left += 1/$length_left;
            }
        }
        foreach my $base (@right_seq) {
            if ($base =~ /G|C/i) {
                $gc_right += 1/$length_right;
            }
        }
        if (($gc_left >= .35) and ($gc_right >= .35)) {                 #subset for GC content > .35
            my $plus_one = $start_pos + 1;
            my $minus_one = $start_pos - 1;
            if (exists $final{"$scaffold_num,$scaffold_length,$minus_one"}) {
                next;
            }
            elsif (exists $final{"$scaffold_num,$scaffold_length,$plus_one"}) {
                delete $final{"$scaffold_num,$scaffold_length,$plus_one"};
                $final{"$scaffold_num,$scaffold_length,$start_pos"} = "$pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length";
            }
            else {
                $final{"$scaffold_num,$scaffold_length,$start_pos"} = "$pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length";
            }
        }
    }    
}


while (my $line = <SGB_SSRS_TRI>) {
    next if ($line =~ /^\D/);
    chomp $line;
    my ($row,$scaffold_num,$scaffold_length,$pattern,$repeat_length,$start_pos,$end_pos,$repeat_char) = split /,/, $line;
    my @left_snp;
    my @right_snp;
    my $gc_right;
    my $gc_left;
    next if ($repeat_length < 8);                              #subset for repeat length >= 8
    
    my $hundred_left = $start_pos - 150;                        #determining start position of 150 bases to the left of the SNP
    if ($hundred_left < 0) {
        $hundred_left = 0;
    }
    my $ssr_start = $start_pos;
        $start_pos += 1;
        my $hundred_right = $end_pos + 3;                            #determining start position of 100 bases to the right of the SNP, works out since substring starts with coordinate 0
        $end_pos += 3;
    my $left_seq = substr $sequences{$scaffold_num}, $hundred_left, 150; #create substring of the 100 bases to the left of the SNP
    my $right_seq = substr $sequences{$scaffold_num}, $hundred_right, 150; #create substring of the 100 bases to the right of the SNP
    my $ssr = substr $sequences{$scaffold_num}, $ssr_start, ($repeat_length*3);
    my @left_seq = split(//, $left_seq);
    my @right_seq = split(//, $right_seq);
    my $full_seq = $left_seq . $ssr . $right_seq;
    my $seq_length = length $full_seq;
    my $length_left = @left_seq;
    my $length_right = @right_seq;
    next if ($length_left < 50);
    next if ($length_right < 50);
    foreach my $base (@left_seq) {                              #calculate gc content of surrounding region
        if ($base =~ /G|C/i) {
            $gc_left += 1/$length_left;
        }
    }
    foreach my $base (@right_seq) {
        if ($base =~ /G|C/i) {
            $gc_right += 1/$length_right;
        }
    }
    if (($gc_left >= .35) and ($gc_right >= .35)) {                 #subset for GC content > .35
        my $plus_one = $start_pos + 1;
        my $minus_one = $start_pos - 1;
        my $minus_two = $start_pos - 2;
        my $plus_two = $start_pos + 2;
        if (exists $final{"$scaffold_num,$scaffold_length,$minus_two"}) {
            next;
        }
        elsif (exists $final{"$scaffold_num,$scaffold_length,$minus_one"}) {
            next;
        }
        elsif (exists $final{"$scaffold_num,$scaffold_length,$plus_one"}) {
            delete $final{"$scaffold_num,$scaffold_length,$plus_one"};
            $final{"$scaffold_num,$scaffold_length,$start_pos"} = "$pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length";
        }
        elsif (exists $final{"$scaffold_num,$scaffold_length,$plus_two"}) {
            delete $final{"$scaffold_num,$scaffold_length,$plus_two"};
            $final{"$scaffold_num,$scaffold_length,$start_pos"} = "$pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length";
        }
        else {
            $final{"$scaffold_num,$scaffold_length,$start_pos"} = "$pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length";
        }
    }
}    


foreach my $key (keys %final) {
    my ($scaffold_num,$scaffold_length,$start_pos) = split /,/, $key;
    my ($pattern,$repeat_length,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length) = split /,/, $final{$key};
    print OUT "$scaffold_num,$scaffold_length,$pattern,$repeat_length,$start_pos,$end_pos,$gc_left,$gc_right,$full_seq,$seq_length\n";
}


sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}


