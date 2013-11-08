#TITLE: split_contigs_gi.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 04/11/13
#PURPOSE: separates contigs into files following BLAST. uses gi list for plant and gi list for microbes.

#!/usr/bin/perl -w
use strict;

my $gi_plant = "/bigdata/broberts/annotation/gi_list_plant.txt";
my $gi_fungi = "/bigdata/broberts/annotation/gi_list_fungi.txt";
my $gi_archaea = "/bigdata/broberts/annotation/gi_list_archaea.txt";
my $gi_bacteria = "/bigdata/broberts/annotation/gi_list_bacteria.txt";
my $gi_viroids = "/bigdata/broberts/annotation/gi_list_viroids.txt";
my $gi_virus = "/bigdata/broberts/annotation/gi_list_virus.txt";

my %plant;
my %fungi;
my %archaea;
my %bacteria;
my %viroids;
my %virus;
my @contigs;
my %some_plant;

my $blast_out = $ARGV[0];
my $fasta = $ARGV[1];
my $output = $ARGV[2];

my $contig_fa;
my $sequence;
my %contig_length;

open FASTA, "<$fasta";
while (my $line = <FASTA>) {                 
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig_fa = $line[0];
        $contig_fa =~ s/>//;
        #print "$contig\n";
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        if (exists $contig_length{$contig_fa}) {
            $contig_length{$contig_fa} += length($sequence);
        }
        else {
            $contig_length{$contig_fa} = length($sequence);
        }
    }
}


open PLANT, "<$gi_plant";

while (my $line = <PLANT>) {
    chomp $line;
    $plant{$line} = 1;
}

open FUNGI, "<$gi_fungi";

while (my $line = <FUNGI>) {
    chomp $line;
    $fungi{$line} = 1;
}

open ARCHAEA, "<$gi_archaea";

while (my $line = <ARCHAEA>) {
    chomp $line;
    $archaea{$line} = 1;
}

open BACTERIA, "<$gi_bacteria";

while (my $line = <BACTERIA>) {
    chomp $line;
    $bacteria{$line} = 1;
}

open VIRUS, "<$gi_virus";

while (my $line = <VIRUS>) {
    chomp $line;
    $virus{$line} = 1;
}

open VIROIDS, "<$gi_viroids";

while (my $line = <VIROIDS>) {
    chomp $line;
    $viroids{$line} = 1;
}

open BLAST, "<$blast_out";

open OUT, ">$output";

print OUT "Contig\tContig_length\tStart_pos\tEnd_pos\tFrag_size\tgi\tType\n";

while (my $line = <BLAST>) {
    chomp $line;
    my @line = split(/\t/,$line);
    my $e_value = $line[-2];
    if ($e_value <= .00001) {
        if ($line =~ /(contig_\d+)\sgi\|(\d+)\|/) {
            my $contig = $1;
            my $gi = $2;
            if (exists $plant{$gi}) {
                $some_plant{$contig} = 1;
            }
        }
    }
}

seek BLAST, 0, 0;

while (my $line = <BLAST>) {
    chomp $line;
    my @line = split(/\t/,$line);
    my $e_value = $line[-2];
    my $start_pos = $line[6];
    my $end_pos = $line[7];
    my $align_length = $line[3];
    if ($e_value <= .00001) {
        if ($line =~ /(contig_\d+)\sgi\|(\d+)\|/) {
            my $contig = $1;
            next if exists $some_plant{$contig};
            my $length = $contig_length{$contig};
            my $gi = $2;
            if (exists $fungi{$gi}) {
                print OUT "$contig\t$length\t$start_pos\t$end_pos\t$align_length\t$gi\tfungi\n"
            }
            elsif (exists $archaea{$gi}) {
                print OUT "$contig\t$length\t$start_pos\t$end_pos\t$align_length\t$gi\tarchaea\n"
            }
            elsif (exists $bacteria{$gi}) {
                print OUT "$contig\t$length\t$start_pos\t$end_pos\t$align_length\t$gi\tbacteria\n"
            }
            elsif (exists $viroids{$gi}) {
                print OUT "$contig\t$length\t$start_pos\t$end_pos\t$align_length\t$gi\tviroids\n"
            }
            elsif (exists $virus{$gi}) {
                print OUT "$contig\t$length\t$start_pos\t$end_pos\t$align_length\t$gi\tvirus\n"
            }
        }
    }
}