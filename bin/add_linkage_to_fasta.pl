#!/usr/bin/perl -w

#TITLE: add_linkage_to_fasta.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/22/2013
#PURPOSE: adds linkage group and cM position of scaffold if marker aligned

use strict;
use warnings;

my $fasta = $ARGV[0];
my $table = $ARGV[1];
my %scaffolds;

open TABLE, "<$table";

while (<TABLE>) {
    my @line = split(/,/, $_);
    my $linkage = $line[0];
    my $position = $line[1];
    my $scaffold = $line[2];
    $scaffolds{$scaffold} = "_$linkage:$position";
}

open FASTA, "<$fasta";
open OUT, ">$ARGV[2]";

while (<FASTA>) {
    chomp;
    if ($_ =~ /^>(.+)\|/) {
        my $scaffold = $1;
        if (exists $scaffolds{$scaffold}) {
            $scaffold = $_ . $scaffolds{$scaffold};
            print OUT ">$scaffold\n";
        }
        else {
            print OUT "$_\n";
        }
    }
    else {
        print OUT "$_\n";
    }
}
