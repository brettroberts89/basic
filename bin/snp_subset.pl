#TITLE: snp_subset.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/28/2012
#PURPOSE: creates a list of SNPs that have a GC content > .35, and either a genomic richness > .375 or a genomic richness > .25 AND no genotype occurs more than 4 times.

#!/usr/bin/perl -w
use strict;

my $input = "snp_surrounding_100seq.csv";
my $vcf = "all.mergedSNPs.vcf";
my $output = "snp_subset_final.csv";
my $test = "snp_subset_gc_unique.csv";
my $test2 = "snp_subset_geno_rich.csv";
my %unique_gc_sub;

open IN, "<$input";
open OUT, ">$output";
open TEST, ">$test";
open TEST2, ">$test2";

print OUT "Possible_alleles,vcf_coord,contig_ID,SNP_position_in_contig,GC_content_left_of_SNP,GC_content_right_of_SNP,Genotypic_Richness\n";
while (my $line = <IN>) {                                                           
    if ($line =~ /unique/) {                                                            #subset for SNPs that are in unique regions
        chomp $line;    
        my @line = split /,/, $line;
        my $vcf_coord = $line[0];
        my $contig = $line[1];
        my $contig_coord = $line[2];
        my $repeat_pct_left = $line[3];
        my $repeat_pct_right = $line[4];
        my $gc_left = $line[5];
        my $gc_right = $line[6];
        if (($gc_left >.35) and ($gc_right >.35)) {                                     #further subset for SNPs with GC content >.35 on both sides
            $unique_gc_sub{$vcf_coord}="$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right";
            print TEST "$vcf_coord,$contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right\n";
        }
    }
}

my $count = keys %unique_gc_sub;
print "$count\n";
open VCF, "<$vcf";
my %geno_rich;

while (my $line = <VCF>) {
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
        next if $ref_base =~ /-/;
        next if $alt_base =~ /-/;
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
        if ($geno_rich > .375) {                                                #subset for SNPs with genotypic richness >.375
            foreach my $allele (keys %possible_alleles) {
                push @alleles, $allele;
            }
            my $alleles = join "/", @alleles;
            $geno_rich{$vcf_coord} = "$alleles,$geno_rich";
            my $genotypes = join (",", @genotypes_translated);
            print TEST2 "$vcf_coord,$alleles,$geno_rich,$genotypes\n"
        }
        elsif ($geno_rich > .25) {                                              #subset for SNPs with genotypic richness >.25 and no genotype occurs more than 4 times
            my $flag = 0;
            foreach my $genotype2 (keys %genotypes) {
                if ($genotypes{$genotype2} > 4) {
                    $flag = 1;
                }
            }
            if ($flag == 0) {
                foreach my $allele (keys %possible_alleles) {
                    push @alleles, $allele;
                }
                my $alleles = join "/", @alleles;
                $geno_rich{$vcf_coord} = "$alleles,$geno_rich";
                my $genotypes = join (",", @genotypes_translated);
                print TEST2 "$vcf_coord,$alleles,$geno_rich,$genotypes\n"
            }
        }
    }
}
my $count2 = keys %geno_rich;
print "$count2\n";
my $count3;

foreach my $gc_unique_coord (keys %unique_gc_sub) {                             #print out union of two subsets
    if (exists $geno_rich{$gc_unique_coord}) {
        $count3 += 1;
        my ($alleles, $geno_richness) = split (/,/, $geno_rich{$gc_unique_coord});
        my ($contig,$contig_coord,$repeat_pct_left,$repeat_pct_right,$gc_left,$gc_right) = split (/,/, $unique_gc_sub{$gc_unique_coord});
        print OUT "$alleles,$gc_unique_coord,$contig,$contig_coord,$gc_left,$gc_right,$geno_richness\n";
    }
}

print "$count3\n";



