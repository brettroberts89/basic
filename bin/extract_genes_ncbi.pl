#TITLE: extract_genes_ncbi.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 6/19/13
#PURPOSE: extracts gene sequences from ncbi summary page

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $output = $input;
$output =~ s/\.txt/genes\.fasta/;

open IN, "<$input";

my $gene;
my $seq;

open OUT, ">$output";

while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /\/gene=\"(.+)\"/) {
        $gene = $1;
    }
    if ($line =~ /\/translation\="(.+)/) {
        $seq = $1;
        $seq =~ s/\W//g;
        #while (my $line = <IN>) {
        #    chomp $line;
        #    $line =~ s/\W//g;
        #    $seq .= $line;
        #    last if ($line =~ /\"$/);
        #}
        print ">$gene\n\n\n"
    }
}