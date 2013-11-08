#TITLE: convert_coords.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/15/2012
#PURPOSE: takes in file with columns: accession number, latitude, and longitude and outputs file with columns accession number, lat_degrees, lat_mins, lat_secs, long_degrees, long_mins, and long_secs.
        ####also converts different geographic units into standard degrees, minutes, seconds format




#!/usr/bin/perl -w
use strict;

my $input = "coordinates.txt";
my $output = "coords_formatted.txt";
open(OUT, ">$output");


open (INPUT, "<$input");

while (my $line = <INPUT>) {
    chomp $line;
    if ($line =~ /^Accession/) {
        print OUT "Accession_Family,Lat_degs,Lat_mins,Lat_secs,Lat_direction,Long_degs,Long_mins,Long_secs,Long_direction\n";
        next;
    }
    my ($accession, $lat, $long);
    my ($lat_deg, $lat_min, $lat_sec, $lat_direction);
    my ($long_deg, $long_min, $long_sec, $long_direction);
    
    ($accession, $lat, $long) = split(/\t/, $line);
    
    if ($lat eq "NA") {         ####case of having no coordinate information
        $lat_deg = "NA";
        $lat_min = "NA";
        $lat_sec = "NA";
    }
    if ($long eq "NA") {
        $long_deg = "NA";
        $long_min = "NA";
        $long_sec = "NA";
    }   
    elsif ($lat =~ /\'/ and $long =~ /\'/) {      ###case of having coordinates in deg,min,sec format
        my $lat_sec_decimal;
        $lat =~ s/^"//;                 #deleting any  " or 0 at the beginning of the degrees
        $long =~ s/^"//;
        $lat =~ s/^0//;
        $long=~ s/^0//;
        $lat =~ s/^\-//;
        $long =~ s/^\-//;
        ($lat_deg, $lat_min, $lat_sec, $lat_direction) = split(/[^0-9.NSEW]+/, $lat);           #splitting latitute info by anything that is not a digit, NSEW, or a period.
        ($long_deg, $long_min, $long_sec, $long_direction) = split(/[^0-9.NSEW]+/, $long);
        if ($lat_min =~ /\./) {
            ($lat_min, my $lat_sec_decimal) = split(/\./, $lat_min);            #splitting into mins and what will be used to calculate seconds
            if (defined $lat_sec_decimal) {                             #appending decimal to the beginning to calculate seconds
                $lat_sec_decimal = "." . "$lat_sec_decimal";
            }
            else {                                                      #setting seconds to zero if it is not defined
                $lat_sec_decimal = 0;
            }
            $lat_sec = ($lat_sec_decimal*60);                           #calculating seconds from decimal   
            $lat_sec = sprintf "%.1f", $lat_sec;                        #adjusting seconds to be in the form of one decimal place
        }
        if ($long_min =~ /\./) {
            ($long_min, my $long_sec_decimal) = split(/\./, $long_min);
            if (defined $long_sec_decimal) {
                $long_sec_decimal = "." . "$long_sec_decimal";
            }
            else {
                $long_sec_decimal = 0;
            }
            $long_sec = ($long_sec_decimal*60);
            $long_sec = sprintf "%.1f", $long_sec;
        }
        if ($lat_sec == "") {                                       #if the coordinates do not specify seconds, making it zero
            $lat_sec = 0;
        }
        if ($long_sec == "") {
            $long_sec = 0;
        }
    }
    elsif ($lat =~ /p/ and $long =~ /p/){    #case of having p coordinate system
        
        
    }
    elsif ($long =~ /\./) {                                 #case of having decimal coordinate system
        $lat =~ s/^\-//;
        $long =~ s/^\-//;
        ($lat_deg, my $lat_decimal) = split(/\./, $lat);            #take everything after the decimal to calculate mins
        if (defined $lat_decimal) {                                 #append decimal to the beginning for correct calculations
            $lat_decimal = "." . "$lat_decimal";
        }
        else {                                                      #setting it equal to zero if there is nothing after the decimal
            $lat_decimal = 0;
        }
        my $lat_min_sec = ($lat_decimal*60);                        #caluclating mins and seconds in decimal
        ($lat_min, my $lat_sec_decimal) = split(/\./, $lat_min_sec);            #splitting into mins and what will be used to calculate seconds
        if (defined $lat_sec_decimal) {                             #appending decimal to the beginning to calculate seconds
            $lat_sec_decimal = "." . "$lat_sec_decimal";
        }
        else {                                                      #setting seconds to zero if it is not defined
            $lat_sec_decimal = 0;
        }
        $lat_sec = ($lat_sec_decimal*60);                           #calculating seconds from decimal   
        $lat_sec = sprintf "%.1f", $lat_sec;                        #adjusting seconds to be in the form of one decimal place
            
        ($long_deg, my $long_decimal) = split(/\./, $long);
        $long_decimal = "." . "$long_decimal";
        my $long_min_sec = ($long_decimal*60);
        ($long_min, my $long_sec_decimal) = split(/\./, $long_min_sec);
        if (defined $long_sec_decimal) {
            $long_sec_decimal = "." . "$long_sec_decimal";
        }
        else {
            $long_sec_decimal = 0;
        }
        $long_sec = ($long_sec_decimal*60);
        $long_sec = sprintf "%.1f", $long_sec;
        if ($lat_sec == 60.0) {
            $lat_min += 1;
            $lat_sec = "0.0";
        }
        if ($long_sec == 60.0) {
            $long_min += 1;
            $long_sec = "0.0";
        }
    }
    
    else{                                           #case of having it in form of just numbers
        $lat_deg = "NA";
        $lat_min = "NA";
        $lat_sec = "NA";

        $long_deg = "NA";
        $long_min = "NA";
        $long_sec = "NA";
    }
    
    
    
    
    print OUT "$accession,$lat_deg,$lat_min,$lat_sec,N,$long_deg,$long_min,$long_sec,W\n"
}
