#TITLE: split_contigs_gi.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 04/11/13
#PURPOSE: separates contigs into files following BLAST. uses gi list for plant and gi list for microbes.

#!/usr/bin/perl -w
use strict;

my $gi_plant = "/bigdata/broberts/annotation/gi_list_plant.txt";
my $gi_microbe = "/bigdata/broberts/annotation/gi_list_microbial.txt";

my %plant;
my %microbe;
my %plant_count;
my %microbe_count;
my %other_count;
my %plant_bitscoresum;
my %plant_escoresum;
my %microbe_bitscoresum;
my %microbe_escoresum;
my @contigs;
my %plant_gi;
my %microbe_gi;

my $blast_out = $ARGV[0];
my $output = $ARGV[1];
my $plant_gi = $ARGV[2];

open PLANT, "<$gi_plant";

while (my $line = <PLANT>) {
    chomp $line;
    $plant{$line} = 1;
}

open MICROBE, "<$gi_microbe";

while (my $line = <MICROBE>) {
    chomp $line;
    $microbe{$line} = 1;
}

for (my $x = 0; $x <= 38481; $x += 1) {
    $plant_count{'contig_$x'} = 0;
    $microbe_count{'contig_$x'} = 0;
    push @contigs, "contig_$x";
}

open BLAST, "<$blast_out";

while (my $line = <BLAST>) {
    chomp $line;
    my @line = split(/\t/,$line);
    my $e_value = $line[-2];
    my $bit_score = $line[-1];
    if ($e_value <= .00001) {
        if ($line =~ /(contig_\d+)\sgi\|(\d+)\|/) {
            my $contig = $1;
            my $gi = $2;
            if (exists $plant{$gi}) {
                $plant_gi{$contig} = "$plant_gi{$contig}" . "$gi\n";
                $plant_escoresum{$contig} += $e_value;
                $plant_bitscoresum{$contig} += $bit_score;
                $plant_count{$contig} += 1;
            }
            elsif (exists $microbe{$gi}) {
                $microbe_gi{$contig} = "$microbe_gi{$contig}" . "$gi\n";
                $microbe_escoresum{$contig} += $e_value;
                $microbe_bitscoresum{$contig} += $bit_score;
                $microbe_count{$contig} += 1;
            }
            else {
                $other_count{$contig} += 1;
            }
        }
    }
}

open OUT, ">$output";
open PLANTGI, ">$plant_gi";

my $contamination_count;

print OUT "contig\tnum_plant\tnum_microbe\tnum_other\n";

foreach my $contig (@contigs) {
    my $plant_avg_escore = 0;
    my $microbe_avg_escore = 0;
    my $plant_avg_bitscore = 0;
    my $microbe_avg_bitscore = 0;
    $plant_avg_escore = $plant_escoresum{$contig}/$plant_count{$contig} if ($plant_count{$contig} > 0);
    $microbe_avg_escore = $microbe_escoresum{$contig}/$microbe_count{$contig} if ($microbe_count{$contig} > 0);
    $plant_avg_bitscore = $plant_bitscoresum{$contig}/$plant_count{$contig} if ($plant_count{$contig} > 0);
    $microbe_avg_bitscore = $microbe_bitscoresum{$contig}/$microbe_count{$contig} if ($microbe_count{$contig} > 0);
    print OUT "$contig\t$plant_count{$contig}\t$plant_avg_escore\t$plant_avg_bitscore\t$microbe_count{$contig}\t$microbe_avg_escore\t$microbe_avg_bitscore\t$other_count{$contig}\n";
    if (($plant_count{$contig} > 0) and ($microbe_count{$contig} > 0)) {
        #print "$contig\t$plant_count{$contig}\t$plant_avg_escore\t$plant_avg_bitscore\t$microbe_count{$contig}\t$microbe_avg_escore\t$microbe_avg_bitscore\t$other_count{$contig}\n";;
        print PLANTGI "$plant_gi{$contig}";
    }
    if (($plant_count{$contig} == 0) and ($microbe_count{$contig} > 0)) {
        $contamination_count += 1;
        print "$contig\n$microbe_gi{$contig}\n";
    }
}
print "Contaminated Contigs = $contamination_count\n";