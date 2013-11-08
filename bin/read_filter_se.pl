#TITLE: read_filter_se.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 06/13/2013
#PURPOSE: filters microbial, chloroplast, mitochondrial, and spike in from reads

#!/usr/bin/perl -w
use strict;

my $fastq = $ARGV[0];
my @fastq_split = split(/\//, $fastq);
my $base = $fastq_split[-1];
$base =~ s/\.fastq//;
my $split_base = $fastq;
$split_base =~ s/\.fastq//;
my $sai;
my $spike_sam = '/bigdata/broberts/contamination/aln_output/spike/bwa_' . $base . '_SPIKE.sam';
my $hg19_sam = '/bigdata/broberts/contamination/aln_output/hg19/bwa_' . $base . '_HG19.sam';
my $jatcp_sam = '/bigdata/broberts/contamination/aln_output/chloroplast/bwa_' . $base . '_Jatcp.sam';
my $jatmt_sam = '/bigdata/broberts/contamination/aln_output/mitochondria/bwa_' . $base . '_Jatmt.sam';
my $microbe_sam = '/bigdata/broberts/contamination/aln_output/microbe/bwa_' . $base . '_MICROBE.sam';
my $mikeslist_sam = '/bigdata/broberts/contamination/aln_output/mikeslist/bwa_' . $base . '_MIKESLIST.sam';


system("module load bwa");

chdir "/bigdata/broberts/contamination/aln_output/spike/";

$sai = $base . '_SPIKE.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/spike/SPIKE $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/spike/SPIKE $sai $fastq > $spike_sam");

system("echo 'Finished mapping to spike in.'");


chdir "/bigdata/broberts/contamination/aln_output/hg19/";

$sai = $base . '_HG19.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/hg19/HG19 $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/hg19/HG19 $sai $fastq > $hg19_sam");

system("echo 'Finished mapping to human genome.'");


chdir "/bigdata/broberts/contamination/aln_output/chloroplast/";

$sai = $base . '_Jatcp.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/chloroplast/Jatcp $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/chloroplast/Jatcp $sai $fastq > $jatcp_sam");

system("echo 'Finished mapping to Jatropha chloroplast.'");


chdir "/bigdata/broberts/contamination/aln_output/mitochondria/";

$sai = $base . '_Jatmt.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mitochondria/MTALL $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/mitochondria/MTALL $sai $fastq > $jatmt_sam");

system("echo 'Finished mapping to Jatropha mitochondria.'");


chdir "/bigdata/broberts/contamination/aln_output/microbe/";

$sai = $base . '_MICROBE.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/microbe/MICROBE $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/microbe/MICROBE $sai $fastq > $microbe_sam");

system("echo 'Finished mapping to microbial contaminants.'");


chdir "/bigdata/broberts/contamination/aln_output/mikeslist/";

$sai = $base . '_MIKESLIST.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mikeslist/MIKESLIST $fastq > $sai");

system("bwa samse /bigdata/broberts/contamination/aln_output/mikeslist/MIKESLIST $sai $fastq > $mikeslist_sam");

system("echo 'Finished mapping to Mikes list of common contaminants.'");

system("echo 'Beginning to split fastq.'");

chdir "/bigdata/broberts/contamination/filtered_reads/";

my $dir = "/bigdata/broberts/contamination/filtered_reads/";

system ("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $spike_sam $fastq");

my $rm_file = $split_base . '_mapped.fastq';
system("rm $rm_file");
my $mv_file = $split_base . '_unmapped.fastq';
my $mv_dest = $dir . $base . '_unmapped.fastq';
system("mv $mv_file $mv_dest");

my $unmapped_fastq = $dir . $base . '_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $hg19_sam $unmapped_fastq");

system("rm $unmapped_fastq");

$unmapped_fastq = $dir . $base . '_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $microbe_sam $unmapped_fastq");

$rm_file = $dir . $base . '_unmapped_unmapped.fastq';
system("rm $rm_file");
$rm_file = $dir . $base . '_unmapped_mapped.fastq';
system("rm $rm_file");

$unmapped_fastq = $dir . $base . '_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $mikeslist_sam $unmapped_fastq");

$rm_file = $dir . $base . '_unmapped_unmapped_unmapped.fastq';
system("rm $rm_file");
$rm_file = $dir . $base . '_unmapped_unmapped_mapped.fastq';
system("rm $rm_file");

$unmapped_fastq = $dir . $base . '_unmapped_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $jatmt_sam $unmapped_fastq");

$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped.fastq';
system("rm $rm_file");
$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_mapped.fastq';
system("rm $rm_file");

$unmapped_fastq = $base . '_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads_single.pl $jatcp_sam $unmapped_fastq");

$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';
system("rm $rm_file");
$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_mapped.fastq';
system("rm $rm_file");

$mv_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';
$mv_dest = $dir . $base . '_no_contam.fastq';
system("mv $mv_file $mv_dest");

$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped_mapped.fastq';
system("rm $rm_file");

system("echo 'Finished filtering reads.'");


