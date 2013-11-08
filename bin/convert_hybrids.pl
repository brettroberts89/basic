#TITLE: convert_hybrids.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/29/2012
#PURPOSE: converts old hybrid cross codes to new hybrid codes

#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
$ARGV[0] =~ s/.csv//;
my $output = "$ARGV[0]" . "_converted.csv";
my $key = "2009-2011_Crosses_ALL.csv";
my %codes;

open IN, "<$input";
open OUT, ">$output";
open KEY, "<$key";

while (<KEY>) {                                         #make hash connecting new hybrid codes and old hybrid codes
    chomp;
    $_ =~ /(\X\d\dH\d\d\d\d).+(\w\w\w\w\sX\s\d+)/;
    my $old_code = $2;
    my $new_code = $1;
    $codes{$old_code} = $new_code;
    #print "$old_code = $new_code\n";
}

while (<IN>) {                                      #search file for occurences of old hybrid codes, convert to new codes
    chomp;
    my $new_code;
    $_ =~ /(\w\w\w\w\sX\s\d+)/;
    #print "$1\n";
    my $old_code = $1;
    if (exists $codes{$old_code}) {
        $new_code = $codes{$old_code};
        $_ =~ s/$old_code/$new_code/;
    }
    print OUT "$_\n";
}