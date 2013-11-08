#TITLE: read_filter_pe.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 06/13/2013
#PURPOSE: filters microbial, chloroplast, mitochondrial, and spike in from reads

#!/usr/bin/perl -w
use strict;

my $fastq1 = $ARGV[0];
my $fastq2 = $ARGV[1];
my @fastq1_split = split(/\//, $fastq1);
my @fastq2_split = split(/\//, $fastq2);
my $base1 = $fastq1_split[-1];
my $base2 = $fastq2_split[-1];

$base1 =~ s/\.fastq//;
$base2 =~ s/\.fastq//;
my $sambase = $base1;
$sambase =~ s/_A//;
my $split_base1 = $fastq1;
my $split_base2 = $fastq2;
$split_base1 =~ s/\.fastq//;
$split_base2 =~ s/\.fastq//;
my $sai1;
my $sai2;
my $spike_sam = '/bigdata/broberts/contamination/aln_output/spike/bwa_' . $sambase . '_SPIKE.sam';
my $hg19_sam = '/bigdata/broberts/contamination/aln_output/hg19/bwa_' . $sambase . '_HG19.sam';
my $jatcp_sam = '/bigdata/broberts/contamination/aln_output/chloroplast/bwa_' . $sambase . '_Jatcp.sam';
my $jatmt_sam = '/bigdata/broberts/contamination/aln_output/mitochondria/bwa_' . $sambase . '_Jatmt.sam';
my $microbe_sam = '/bigdata/broberts/contamination/aln_output/microbe/bwa_' . $sambase . '_MICROBE.sam';
my $mikeslist_sam = '/bigdata/broberts/contamination/aln_output/mikeslist/bwa_' . $sambase . '_MIKESLIST.sam';


system("module load bwa");

chdir "/bigdata/broberts/contamination/aln_output/spike/";

$sai1 = $base1 . '_SPIKE.sai';
$sai2 = $base2 . '_SPIKE.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/spike/SPIKE $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/spike/SPIKE $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/spike/SPIKE $sai1 $sai2 $fastq1 $fastq2 > $spike_sam");

system("echo 'Finished mapping to spike in.'");


chdir "/bigdata/broberts/contamination/aln_output/hg19/";

$sai1 = $base1 . '_HG19.sai';
$sai2 = $base2 . '_HG19.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/hg19/HG19 $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/hg19/HG19 $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/hg19/HG19 $sai1 $sai2 $fastq1 $fastq2 > $hg19_sam");

system("echo 'Finished mapping to human genome.'");


chdir "/bigdata/broberts/contamination/aln_output/chloroplast/";

$sai1 = $base1 . '_Jatcp.sai';
$sai2 = $base2 . '_Jatcp.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/chloroplast/Jatcp $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/chloroplast/Jatcp $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/chloroplast/Jatcp $sai1 $sai2 $fastq1 $fastq2 > $jatcp_sam");

system("echo 'Finished mapping to Jatropha chloroplast.'");


chdir "/bigdata/broberts/contamination/aln_output/mitochondria/";

$sai1 = $base1 . '_Jatmt.sai';
$sai2 = $base2 . '_Jatmt.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mitochondria/MTALL $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mitochondria/MTALL $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/mitochondria/MTALL $sai1 $sai2 $fastq1 $fastq2 > $jatmt_sam");

system("echo 'Finished mapping to Jatropha mitochondria.'");


chdir "/bigdata/broberts/contamination/aln_output/microbe/";

$sai1 = $base1 . '_MICROBE.sai';
$sai2 = $base2 . '_MICROBE.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/microbe/MICROBE $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/microbe/MICROBE $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/microbe/MICROBE $sai1 $sai2 $fastq1 $fastq2 > $microbe_sam");

system("echo 'Finished mapping to microbial contaminants.'");


chdir "/bigdata/broberts/contamination/aln_output/mikeslist/";

$sai1 = $base1 . '_MIKESLIST.sai';
$sai2 = $base2 . '_MIKESLIST.sai';

system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mikeslist/MIKESLIST $fastq1 > $sai1");
system("bwa aln -t 8 /bigdata/broberts/contamination/aln_output/mikeslist/MIKESLIST $fastq2 > $sai2");

system("bwa sampe /bigdata/broberts/contamination/aln_output/mikeslist/MIKESLIST $sai1 $sai2 $fastq1 $fastq2 > $mikeslist_sam");

system("echo 'Finished mapping to Mikes list of common contaminants.'");

system("echo 'Beginning to split fastq.'");

chdir "/bigdata/broberts/contamination/filtered_reads/";

my $dir = "/bigdata/broberts/contamination/filtered_reads/";

system ("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $spike_sam $fastq1 $fastq2");

#my $rm_file = $split_base1 . '_mapped.fastq';
#system("rm $rm_file");
my $mv_file1 = $split_base1 . '_unmapped.fastq';
my $mv_file2 = $split_base2 . '_unmapped.fastq';
my $mv_dest1 = $dir . $base1 . '_unmapped.fastq';
my $mv_dest2 = $dir . $base2 . '_unmapped.fastq';
system("mv $mv_file1 $mv_dest1");
system("mv $mv_file2 $mv_dest2");

my $unmapped_fastq1 = $dir . $base1 . '_unmapped.fastq';
my $unmapped_fastq2 = $dir . $base2 . '_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $hg19_sam $unmapped_fastq1 $unmapped_fastq2");

#system("rm $unmapped_fastq");

$unmapped_fastq1 = $dir . $base1 . '_unmapped_unmapped.fastq';
$unmapped_fastq2 = $dir . $base2 . '_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $microbe_sam $unmapped_fastq1 $unmapped_fastq2");

#$rm_file = $dir . $base . '_unmapped_unmapped.fastq';
#system("rm $rm_file");
#$rm_file = $dir . $base . '_unmapped_mapped.fastq';
#system("rm $rm_file");

$unmapped_fastq1 = $dir . $base1 . '_unmapped_unmapped_unmapped.fastq';
$unmapped_fastq2 = $dir . $base2 . '_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $mikeslist_sam $unmapped_fastq1 $unmapped_fastq2");

#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped.fastq';
#system("rm $rm_file");
#$rm_file = $dir . $base . '_unmapped_unmapped_mapped.fastq';
#system("rm $rm_file");

$unmapped_fastq1 = $dir . $base1 . '_unmapped_unmapped_unmapped_unmapped.fastq';
$unmapped_fastq2 = $dir . $base2 . '_unmapped_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $jatmt_sam $unmapped_fastq1 $unmapped_fastq2");

#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped.fastq';
#system("rm $rm_file");
#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_mapped.fastq';
#system("rm $rm_file");

$unmapped_fastq1 = $base1 . '_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';
$unmapped_fastq2 = $base2 . '_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';

system("perl /bigdata/dilaraally/mtgenome/scripts/split_mapped_reads.pl $jatcp_sam $unmapped_fastq1 $unmapped_fastq2");

#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';
#system("rm $rm_file");
#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_mapped.fastq';
#system("rm $rm_file");

#$mv_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped_unmapped.fastq';
#$mv_dest = $dir . $base . '_no_contam.fastq';
#system("mv $mv_file $mv_dest");

#$rm_file = $dir . $base . '_unmapped_unmapped_unmapped_unmapped_unmapped_mapped.fastq';
#system("rm $rm_file");

system("echo 'Finished filtering reads.'");

