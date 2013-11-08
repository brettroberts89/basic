#TITLE: determine_repeats_in_reads.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/11/2012
#PURPOSE: searches a fastq file for mono-, di- and trinucleotide repeats in the sequences. Also outputs percent of sequence that is repetitive.

#!/usr/bin/perl -w
use strict;

my @input = @ARGV;
my %unmapped;
my $file = $ARGV[0];
my $unmapped = $ARGV[1];
open UNMAPPED, "<$unmapped";


#create list of all possible repeats                                            
my @mono_probes = ("A","C","G","T");
my @di_probes = ("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG");
my @tri_probes = ("ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATC","ATG","ATA","ATT","CAG","CAC","CAA","CAT","CGC","CGG","CGA","CGT","CTC","CTG","CTA","CTT","GAC","GAG","GAA","GAT","GCC","GCG","GCA","GCT","GTC","GTG","GTA","GTT","TAC","TAG","TAA","TAT","TCC","TCG","TCA","TCT","TGC","TGG","TGA","TGT");


while (my $line = <UNMAPPED>) {
    if ($line !~ /^seqID/) {
        my @line = split /,/, $line;
        my $seq_id = $line[0];
        my $first_mapped = $line[6];
        my $second_mapped = $line[7];
        if ($first_mapped == 0) {
            $unmapped{"$seq_id/1"} = 1;
        }
        if ($second_mapped == 0) {
            $unmapped{"$seq_id/2"} = 1;
        }
    }
}

my $line_count;
print "Starting to read file: $file\n";
open IN, "<$file";
my @file = split /\//, $file;
my $file_mod = $file[-1];
$file_mod =~ s/.fastq/_repeats.csv/;
my $output = "/shared/sgbiofuels/repeat_analysis/output/$file_mod";
open OUT, ">$output";
print OUT "file,seq_ID,repeat_pct,repeat_type,repeat_length,start_pos\n";
while (my $line_1 = <IN>) {
    chomp $line_1;
    $line_count += 1;
    $line_1 =~ s/\@//;
    my $seq_id = $line_1;
    my $seq;
    my $plus;
    my $quality;
    my %repeats;
    while (my $line_2 = <IN>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <IN>) {
            chomp $line_3;
            $plus = $line_3;
            print "error\n" if ($plus ne "+");
            while (my $line_4 = <IN>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if (!exists $unmapped{$seq_id});
    my $repeat_pct = 0;
    my @repeat_count;
    my %repeat_check;
    my @seq = split //, $seq;
    my @count = (0 .. (length($seq) - 1));
    my %mono_skip;
    my %di_skip;
    my %tri_skip;
    
    foreach my $probe (@mono_probes) {                          #goes through each probe
        my @skip;
        print "Searching for probe: $probe in read: $seq_id with sequence: $seq in file: $file\n";
        foreach my $count (@count) {                            #goes through each position in the sequence
            if (exists $skip[$count]) {                         #if the sequence is already in a repeat, skip
                next;
            }
            else {      
                my $start_pos =$count;           #determine start position of reduced sequence, starting at current position
                my $region = substr $seq, $count;              #determine sequence starting at current position and continuing to the end of the sequence
                if ($region =~ /^(($probe){5,})/i) {                #searching at start of the region for probe sequence, with at least 5 repeats
                    my $repeat_length = length($1);                 #calculate length of the repeat
                    my @repeat = split(//,$1);
                    @repeat_count = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count) {
                        if (!exists $repeat_check{$base}){
                            $repeat_pct += 1/(length($seq));
                        }
                    }
                    push @skip, @repeat;                            #doesn't check following start positions that are already a part of the repeat
                    $repeats{"$file_mod,$seq_id,mono,$repeat_length,$start_pos"} = 1;     
                }
                else {                                              #if no match, goes on to next position
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@di_probes) {                                #same with dinucleotide repeats
        my @skip;
        print "Searching for probe: $probe in read: $seq_id with sequence: $seq in file: $file\n";
        foreach my $count (@count) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $count + 1;
                my $region = substr $seq, $count;              #determine sequence starting at current position and continuing to the end of the sequence
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/2);
                    my @repeat = split(//,$1);
                    @repeat_count = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count) {
                        if (!exists $repeat_check{$base}){
                            $repeat_pct += 1/(length($seq));
                        }
                    }
                    pop @repeat;
                    push @skip, @repeat;
                    push @skip, @repeat;                            #doesn't check following start positions that are already a part of the repeat
                    $repeats{"$file_mod,$seq_id,di,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    foreach my $probe (@tri_probes) {                                   #same, with trinucleotide repeats
        my @skip;
        print "Searching for probe: $probe in read: $seq_id with sequence: $seq in file: $file\n";
        foreach my $count (@count) {
            if (exists $skip[$count]) {
                next;
            }
            else {
                my $start_pos = $count + 1;
                my $region = substr $seq, $count;              #determine sequence starting at current position and continuing to the end of the sequence
                if ($region =~ /^(($probe){5,})/i) {
                    my $repeat_length = (length($1)/3);
                    my @repeat = split(//,$1);
                    @repeat_count = ($count..($count+(length $1)-1));
                    foreach my $base (@repeat_count) {
                        if (!exists $repeat_check{$base}){
                            $repeat_pct += 1/(length($seq));
                        }
                    }
                    pop @repeat;
                    pop @repeat;
                    push @skip, @repeat;
                    push @skip, @repeat;                            #doesn't check following start positions that are already a part of the repeat
                    $repeats{"$file_mod,$seq_id,tri,$repeat_length,$start_pos"} = 1;
                }
                else {
                    $skip[$count] = 1;
                }
            }
        }
    }
    if ($repeat_pct == 0) {
        print OUT "$file_mod,$seq_id,0,unique,unique,unique\n";
    }
    else {
        foreach my $repeat (keys %repeats) {
            my ($file_a,$seq_id_a,$type_a,$repeat_length_a,$start_pos_a) = split /,/, $repeat;
            if ($type_a eq "di") {
                foreach my $comp_key (keys %repeats) {
                    my ($comp_file,$comp_seq_id,$comp_type,$comp_repeat_length,$comp_start_pos) = split /,/, $comp_key;
                    if (($comp_file eq $file_a) and ($comp_seq_id eq $seq_id_a) and ($comp_start_pos eq ($start_pos_a-1))) {
                        delete $repeats{$repeat};
                    }
                }
            }
            if ($type_a eq "tri") {
                foreach my $comp_key (keys %repeats) {
                    my ($comp_file,$comp_seq_id,$comp_type,$comp_repeat_length,$comp_start_pos) = split /,/, $comp_key;
                    if (($comp_file eq $file_a) and ($comp_seq_id eq $seq_id_a) and (($comp_start_pos eq ($start_pos_a-1)) or ($comp_start_pos eq ($start_pos_a-2)))) {
                        delete $repeats{$repeat};
                    }
                }
            }
        }
        foreach my $repeat (sort(keys %repeats)) {
            my ($file_a,$seq_id_a,$type_a,$repeat_length_a,$start_pos_a) = split /,/, $repeat;
            print OUT"$file_a,$seq_id_a,$repeat_pct,$type_a,$repeat_length_a,$start_pos_a\n";
        }
    }
}
        