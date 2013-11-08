#TITLE: replace_misassemblies.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 11/05/2012
#PURPOSE: takes in files with list of contigs in each linkage group and a file with misassembled scaffolds, and replaces the misassembled scaffold names with contigs and extracts the fasta entries for both scaffolds and contigs.

#!/usr/bin/perl -w
use strict;

my $final_set = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/linkage_map/final_set_problems.csv";
my $contigs = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/SGB_assembly/SGBv1contigs.fasta";
my $scaffolds = "/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/SGB_assembly/SGB.v1.fasta";
my @linkage_groups = @ARGV;
my %misassembly;
open SET, "<$final_set";

while (my $line = <SET>) {
    chomp $line;
    $line =~ s/\"//g;
    my ($scaffold, $linkage_group, $marker, $f_pos, $cm_pos, $contig, %contig_pos) = split(/,/, $line);
    next if ($misassembly{"$scaffold,$linkage_group"} =~ /$contig/);
    $misassembly{"$scaffold,$linkage_group"} = $misassembly{"$scaffold,$linkage_group"} . "$contig;";
}
foreach my $linkage_group (@linkage_groups) {
    open IN, "<$linkage_group";
    my @file_split = split(/_/, $linkage_group);
    my $lg = $file_split[0];
    my $output = "$lg.fa";
    open OUT, ">$output";
    my %reference;
    while (my $line = <IN>) {
        chomp $line;
        $line =~ s/\"//g;
        next if $line =~ "^Scaffold";
        my @line = split(/,/, $line);
        my $scaffold = $line[0];
        my $linkage_group = $line[5];
        $reference{"$scaffold,$linkage_group"} = 1;
    }
    foreach my $key (keys %reference) {
        my ($scaffold,$linkage_group) = split(/,/, $key);
        print "$scaffold\n";
        if (exists $misassembly{"$key"}) {
            my @contigs = split(/;/, $misassembly{"$scaffold,$linkage_group"});
            open CONTIGS, "<$contigs";
            foreach my $contig (@contigs) {
                while (my $fa_line = <CONTIGS>) {
                    chomp $fa_line;
                    if ($fa_line =~ /($contig)/) {
                        print OUT "\n$fa_line\n";
                        while (my $fa_line2 = <CONTIGS>) {
                            chomp $fa_line2;
                            if ($fa_line2 !~ /\>/) {
                                print OUT "$fa_line2";
                            }
                            else {
                                last;
                            }
                        }
                        last;
                    }
                }
            }
        }
        else {
            open SCAFFOLDS, "<$scaffolds";
            while (my $fa_line = <SCAFFOLDS>) {
                chomp $fa_line;
                if ($fa_line =~ /($scaffold)/) {
                    print OUT "\n$fa_line\n";
                    while (my $fa_line2 = <SCAFFOLDS>) {
                        chomp $fa_line2;
                        if ($fa_line2 !~ /\>/) {
                            print OUT "$fa_line2";
                        }
                        else {
                            last;
                        }
                    }
                    last;
                }
            }
        }
    }
}