#TITLE: determine_surrounding_100seq.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2012
#PURPOSE: surveys the region 100 bp on either side of a snp for mono-, di-, and trinucleotide repeats of at least 5. Outputs the location and length of such repeat regions.

#!/usr/bin/perl -w
use strict;

my $input = "contig_coords.csv";
my $sato = "Sato.fa";
my $output = "snp_surrounding_100seq.csv";
my $contig;
my $sequence;
my %snps;
my %sequences;
my %repeats;
my %uniques;
my @hundred = (0..99);                                        #create array of numbers 0 to 100 as counter


open OUT, ">$output";
print OUT "vcf_coord,Contig,SNP_position_in_contig,Pct_repetitive_left,Pct_repetitive_right,GC_content_left_of_SNP,GC_content_right_of_SNP,type_of_repeat,repeat_length,start_position\n";
open SATO, "<$sato";

while (my $sato_line = <SATO>) {                   #go through Sato assembly, creating a hash entry for each contig with the sequence as the value
    chomp $sato_line;
    if ($sato_line =~ /^>/) {                    #starts new hash entry if line starts with >
        $contig = $sato_line;
        $contig =~ s/>//;
        $contig =~ s/\D+//;
    }
    else{           #adds on to length of contig while the the line begins with a nucleotide    
        $sequence = $sato_line;
        $sequence =~ s/\s+//;
        $sequences{$contig} = "$sequences{$contig}" . "$sato_line";
    }
}


#create list of all possible repeats                                            
my @mono_probes = ("A","C","G","T");
my @di_probes = ("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG");
my @tri_probes = ("ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATC","ATG","ATA","ATT","CAG","CAC","CAA","CAT","CGC","CGG","CGA","CGT","CTC","CTG","CTA","CTT","GAC","GAG","GAA","GAT","GCC","GCG","GCA","GCT","GTC","GTG","GTA","GTT","TAC","TAG","TAA","TAT","TCC","TCG","TCA","TCT","TGC","TGG","TGA","TGT");

open IN, "<$input";

#go through each SNP
while (my $line = <IN>) {
    chomp $line;
    next if ($line =~ /^\D/);
    my @left_snp;
    my @right_snp;
    my $gc_right;
    my $gc_left;
    my $repeat_pct_left = 0;
    my $repeat_pct_right = 0;
    my @repeat_count_left;
    my @repeat_count_right;
    my %repeat_check_left;
    my %repeat_check_right;
    my ($coord, $contig, $contig_coord, $contig_length) = split(/,/, $line);
    my $hundred_left = $contig_coord - 101;                        #determining start position of 100 bases to the left of the SNP
    if ($hundred_left < 0) {
        $hundred_left = 0;
    }
    my $hundred_right = $contig_coord;                            #determining start position of 100 bases to the right of the SNP, works out since substring starts with coordinate 0
    my $left_seq = substr $sequences{$contig}, $hundred_left, 100; #create substring of the 100 bases to the left of the SNP
    my $right_seq = substr $sequences{$contig}, $hundred_right, 100; #create substring of the 100 bases to the right of the SNP
    my @left_seq = split(//, $left_seq);
    my @right_seq = split(//, $right_seq);
    foreach my $base (@left_seq) {                              #calculate gc content of surrounding region
        if ($base =~ /G|C/i) {
            $gc_left += .01;
        }
    }
    foreach my $base (@right_seq) {
        if ($base =~ /G|C/i) {
            $gc_right += .01;
        }
    }

    
    foreach my $probe (@mono_probes) {                          #goes through each probe
        my @skip;
        foreach my $count (@hundred) {                            #goes through each position in the sequence
            if (exists $skip[$count]) {                         #if the sequence is already in a repeat, skip
                next;
            }
            else {      
                my $start_pos = $hundred_left + $count + 1;           #determine start position of reduced sequence, starting at current position
                my $region = substr $left_seq, $count;              #determine sequence starting at current position and continuing to the end of the sequence
                if ($region =~ /^(($probe){5,})/i) {                #searching at start of the region for probe sequence, with at least 5 repeats
                    my $repeat_length = length($1);                 #calculate length of the repeat
                    my @repeat = split(//,$1);
                    @repeat_count_left = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_left) {
                        if (!exists $repeat_check_left{$base}){
                            $repeat_pct_left += .01;
                        }
                    }
                    push @skip, @repeat;                            #doesn't check following start positions that are already a part of the repeat
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,mono,$repeat_length,$start_pos"} = 1;     
                }
                else {                                              #if no match, goes on to next position
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@di_probes) {                                #same with dinucleotide repeats
        my @skip;
        foreach my $count (@hundred) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $hundred_left + $count + 1;
                my $region = substr $left_seq, $count;
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/2);
                    my @repeat = split(//,$1);
                    @repeat_count_left = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_left) {
                        if (!exists $repeat_check_left{$base}){
                            $repeat_pct_left += .01;
                        }
                    }
                    pop @repeat;
                    push @skip, @repeat;
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,di,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@tri_probes) {                                   #same, with trinucleotide repeats
        my @skip;
        foreach my $count (@hundred) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $hundred_left + $count + 1;
                my $region = substr $left_seq, $count;
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/3);
                    my @repeat = split(//,$1);
                    @repeat_count_left = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_left) {
                        if (!exists $repeat_check_left{$base}){
                            $repeat_pct_left += .01;
                        }
                    }
                    pop @repeat;
                    pop @repeat;
                    push @skip, @repeat;
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,tri,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    
    
    
    foreach my $probe (@mono_probes) {                                              #on right side, with mononucleotide repeats
        my @skip;
        foreach my $count (@hundred) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $hundred_right + $count + 1;
                my $region = substr $right_seq, $count;
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = length($1);
                    my @repeat = split(//,$1);
                    @repeat_count_right = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_right) {
                        if (!exists $repeat_check_right{$base}){
                            $repeat_pct_right += .01;
                        }
                    }
                    push @skip, @repeat;
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,mono,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@di_probes) {                                            #on right side, with dinucleotide repeats
        my @skip;
        foreach my $count (@hundred) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $hundred_right + $count + 1;
                my $region = substr $right_seq, $count;
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/2);
                    my @repeat = split(//,$1);
                    @repeat_count_right = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_right) {
                        if (!exists $repeat_check_right{$base}){
                            $repeat_pct_right += .01;
                        }
                    }
                    pop @repeat;
                    push @skip, @repeat;
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,di,$repeat_length,$start_pos"} = 1;                    
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@tri_probes) {                                           #on right side, with trinucleotide repeats
        my @skip;
        foreach my $count (@hundred) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $hundred_right + $count + 1;
                my $region = substr $right_seq, $count;
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/3);
                    my @repeat = split(//,$1);
                    @repeat_count_right = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count_right) {
                        if (!exists $repeat_check_right{$base}){
                            $repeat_pct_right += .01;
                        }
                    }
                    pop @repeat;
                    pop @repeat;
                    push @skip, @repeat;
                    $repeats{"$coord,$contig,$contig_coord,$gc_left,$gc_right,tri,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    $snps{"$coord,$contig,$contig_coord,$gc_left,$gc_right"} = "$repeat_pct_left,$repeat_pct_right";
}


foreach my $snp (keys %snps) {                                                       #checking that each snp has at least one entry; if it doesn't, an entry specifying unique is made
    my ($coord,$contig, $contig_coord,$gc_left,$gc_right) = split /,/, $snp;
    my ($repeat_pct_left, $repeat_pct_right) = split /,/, $snps{$snp};
    $uniques{"$coord,$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right,unique,unique,unique"} = 1;
    foreach my $repeat (keys %repeats) {
        my @array = split /,/, $repeat;
        my $contig_repeat = $array[1];
        my $contig_coord_repeat = $array[2];
        if ($contig_repeat == $contig and $contig_coord_repeat == $contig_coord) {
            $repeats{$repeat} = "$repeat_pct_left,$repeat_pct_right";
            delete $uniques{"$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right,unique,unique,unique"};
        }
    }    
}
foreach my $unique (sort(keys %uniques)) {
    print OUT "$unique\n";
}
foreach my $repeat (sort(keys %repeats)) {
    my ($coord,$contig,$contig_coord,$gc_left,$gc_right,$type,$repeat_length,$start_pos) = split /,/, $repeat;
    my ($repeat_pct_left,$repeat_pct_right) = split /,/, $repeats{$repeat};
    print OUT"$coord,$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right,$type,$repeat_length,$start_pos\n";
}
