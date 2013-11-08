#TITLE: snp_subset.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 08/28/2012
#PURPOSE: creates a list of SNPs that have a GC content > .35, and either a genomic richness > .375 or a genomic richness > .25 AND no genotype occurs more than 4 times.

#!/usr/bin/perl -w

my %snps;
my $input = "final_96SNPs_GBS_validation.csv";
my $vcf = "all.mergedSNPs.vcf";
my $output = "final_96SNPs_genotypes.csv";

open IN, "<$input";
open VCF, "<$vcf";
open OUT, ">$output";
print OUT "vcf_coord,EG6_2,TM704_2,TM338_7,MCPINDIA5C_B,VV12_10,SL851_3,TM241_22\n";

while (my $line = <IN>) {                                                       #get list of vcf coordinates of SNPs that we want info on
    chomp $line;
    my @line = split /,/, $line;
    $snps{$line[1]} = 1;
}

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
        foreach my $snp (keys %snps) {                                          #obtain genotype information
            if ($snp == $vcf_coord) {
                my $ref_base = $line[3];
                my $alt_base = $line[4];
                my @alt_base = split /,/,$alt_base;
                my @genotypes = ($line[10],$line[28],$line[30],$line[12],$line[32],$line[14],$line[19]);              #list of the samples that we want genotypes for
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
                }
                my $genotypes_final = join(",", @genotypes_translated);
                print OUT "$vcf_coord,$genotypes_final\n";
            }
        }

    }
}