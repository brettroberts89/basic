#TITLE: snp_random_subset.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/28/2012
#PURPOSE: creates output file with genotypic richness, GC content of surrounding regions, and coordinate information for a random subset of SNPs
#!/usr/bin/perl -w
use strict;

my $subset = "random_SNPs.csv";
my $input = "snp_surrounding_100seq.csv";
my $vcf = "all.mergedSNPs.vcf";
my $output = "snp_random_subset_final.csv";
my $test = "snp_subset_1.csv";
my $test2 = "snp_subset_2.csv";
my %subset;
my %subset2;

open SUBSET, "<$subset";
open IN, "<$input";
open OUT, ">$output";
open TEST, ">$test";
open TEST2, ">$test2";

while (my $line = <SUBSET>) {                                                   #create hash of all SNPs in random subset
    chomp $line;
    my ($vcf_pos,$contig,$contig_pos,$contig_length) = split(/,/, $line);
    $subset{$vcf_pos} = "$contig,$contig_pos,$contig_length";
}

print OUT "Possible_alleles,vcf_coord,contig_ID,SNP_position_in_contig,Repeat_pct_left_of_SNP,Repeat_pct_right_of_SNP,GC_content_left_of_SNP,GC_content_right_of_SNP,Genotypic_Richness\n";
while (my $line = <IN>) {                                                       #obtain information about surrounding sequence from snp_surrounding_100seq.csv                                       
    chomp $line;    
    my @line = split /,/, $line;
    my $vcf_coord = $line[0];
    my $contig = $line[1];
    my $contig_coord = $line[2];
    my $repeat_pct_left = $line[3];
    my $repeat_pct_right = $line[4];
    my $gc_left = $line[5];
    my $gc_right = $line[6];
    if (exists $subset{$vcf_coord}) {
        $subset2{$vcf_coord}="$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right";
        print TEST "$vcf_coord,$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right\n";
    }
}

my $count = keys %subset2;
my $count2;
print "$count\n";
open VCF, "<$vcf";
my %geno_rich;

while (my $line = <VCF>) {                                                      #obtain genotypic information from vcf file
    if ($line =~ /^1/) {    
        my %genotypes;
        my $possible_genotypes;
        my %possible_alleles;
        my @genotypes_translated;
        my @alleles;
        chomp $line;
        my @line = split(/\t/, $line);
        my $vcf_coord = $line[1];
        my $ref_base = $line[3];
        my $alt_base = $line[4];
        if (exists $subset2{$vcf_coord}) {
            my ($contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right) = split /,/, $subset2{$vcf_coord};
            my @alt_base = split /,/,$alt_base;
            my @genotypes = ($line[10],$line[28],$line[30],$line[12],$line[20],$line[47],$line[32],$line[14]);
            foreach my $sample (@genotypes) {                                       #determine genotypes of different samples from vcf file
                my @sample = split /\D+/, $sample;
                my $base_a;
                my $base_b;
                $base_a = $ref_base if ($sample[0] == "0");
                $base_a = $alt_base[0] if ($sample[0] == "1");
                $base_a = $alt_base[1] if ($sample[0] == "2");
                $base_a = $alt_base[2] if ($sample[0] == "3");
                $base_b = $ref_base if ($sample[1] == "0");
                $base_b = $alt_base[0] if ($sample[1] == "1");
                $base_b = $alt_base[1] if ($sample[1] == "2");
                $base_b = $alt_base[2] if ($sample[1] == "3");
                $genotypes{"$base_a/$base_b"} += 1;
                $possible_alleles{$base_a} = 1;
                $possible_alleles{$base_b} = 1;
                push @genotypes_translated, "$base_a/$base_b";
            }                                                                       #count number of different genotypes
            foreach my $genotype (keys %genotypes) {
                $possible_genotypes += 1;
            }
            my $geno_rich = $possible_genotypes/8;                                  #calculate genotypic richness
            foreach my $allele (keys %possible_alleles) {
                push @alleles, $allele;
            }
            my $alleles = join "/", @alleles;
            my $genotypes = join (",", @genotypes_translated);
            print OUT "$alleles,$vcf_coord,$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right,$geno_rich\n";
            $count2 += 1;
            print TEST2 "$vcf_coord,$alleles,$geno_rich,$genotypes\n"
        }
    }
}
print "$count2\n";