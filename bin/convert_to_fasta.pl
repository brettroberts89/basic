#TITLE: convert_to_fasta.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 09/05/12
#PURPOSE: To create an .fa file from an excel spreadsheet with columns name and sequence
#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];
my $output= $ARGV[0];
$output =~ s/csv/fa/;
open(OUT, ">$output");

open (FILE, "<$file");
while (my $line = <FILE>) {
    chomp $line;
    my @line = split(/,/, $line);
    my $id = $line[0];
    my $seq = $line[1];
    print OUT ">$id\n$seq\n";
    }
    