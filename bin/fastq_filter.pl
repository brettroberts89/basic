#!/usr/bin/perl -w

#TITLE: fastq_filter.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 07/22/2013
#PURPOSE: filters reads based on length, quality, duplication, and contamination


use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $help = '';
my $se = '';
my $pe = '';
my $d = '';
my $q = '';
my $l = '';
my $m = '';
my $c = '';
my $o = '';
my $s = '';
my $threads = 1;
my $min_qual = 17;
my $pct = 70;
my $min_length = 36;
my $opt_threads = '';
my $opt_min_qual = '';
my $opt_pct = '';
my $opt_min_length = '';

GetOptions('se' => \$se, 'pe' => \$pe, 'quality|q' => \$q, 'microbe|m' => \$m, 'chloroplast|c' => \$c, 'mitochondria|o' => \$o,  'pct|p' => \$opt_pct, 'min_qual|y' => \$opt_min_qual, 'min_length|l' => \$opt_min_length, 'threads|t:i' => \$threads, 'help|h' => \$help);

$threads = $opt_threads if ($opt_threads);
$min_qual = $opt_min_qual if ($opt_min_qual);
$pct = $opt_pct if ($opt_pct);
$min_length = $opt_min_length if ($opt_min_length);
$min_qual = $opt_min_qual if ($opt_min_qual);

my $usage = <<USAGE;
    
    Usage: perl fastq_filter.pl --se/--pe [OPTS] reads1.fastq reads2.fastq
    
    --se:               Single-end reads
    --pe:               Paired-end reads
    --quality/-q        quality/length filter
    --microbe/-m        microbial filter
    --chloroplast/-c    chloroplast filter
    --mitochondria/-o   mitochondrial filter
    --threads/-t:       Number of threads [1]
    --min_qual/-y:       Qualities below this value will be trimmed from the end of the reads [17]
    --pct/-p:           Minimum percent of bases in a read that must have [-q] quality [70]
    --min_length/-l:    Reads shorter than this after trimming will be discarded [36]

    
    *****Required modules: fastx_toolkit, bwa
    
USAGE

if (((!$se) and (!$pe)) or (($se) and ($pe)) or ($help) or ((!$q) and (!$l) and (!$m) and (!$c) and (!$o) and (!$s))) {
    print $usage;
}

if ($se) {
    my $fastq = $ARGV[-1];
    if ($fastq !~ /\.fastq/) {
        die "\nInvalid fastq file!\n$usage";
    }
    my $dir = getcwd;
    if ($fastq !~ /^\//) {
        $fastq = $dir . '/' . $fastq;
    }
    my $fastq_mid = $fastq;
    my @fastq_split = split(/\//, $fastq);
    my $base = $fastq_split[-1];
    $base =~ s/\.fastq//;
    $fastq_mid =~ s/\.fastq/_dup_length_qual_filt\.fastq/;
    my $fastq_dup_filter = $fastq;
    $fastq_dup_filter =~ s/\.fastq/_dup_filt\.fastq/;
    my $hg19_sam = '/shared/sgbiofuels/contamination/aln_output/hg19/bwa_' . $base . '_HG19.sam';
    my $microbe_sam = '/shared/sgbiofuels/contamination/aln_output/microbe/bwa_' . $base . '_MICROBE.sam';
    my $mikeslist_sam = '/shared/sgbiofuels/contamination/aln_output/mikeslist/bwa_' . $base . '_MIKESLIST.sam';
    my $jatcp_sam = '/shared/sgbiofuels/contamination/aln_output/jatcp/bwa_' . $base . '_JATCP.sam';
    my $jatmt_sam = '/shared/sgbiofuels/contamination/aln_output/jatmt/bwa_' . $base . '_JATMT.sam';
    my $spike_sam = '/shared/sgbiofuels/contamination/aln_output/spike/bwa_' . $base . '_SPIKE.sam';
    
    
#########   Duplicate, length, and quality filter     
    
    print "Beginning DNA filter of $fastq with $threads threads.\n";
    if ($d) {
        system("python /shared/sgbiofuels/contamination/bin/fastq_unique.py < $fastq > $fastq_dup_filter");
        print "Finished duplicate filtering.\n";
    }
    else {
        $fastq_dup_filter = $fastq;
    }
    if ($q) {
        system("fastq_quality_trimmer -t $min_qual -l $min_length -i $fastq_dup_filter | fastq_quality_filter -q $min_qual -p $pct > $fastq_mid");
        print "Finished length and quality filtering/trimming\n";
    }
    else {
        $fastq_mid = $fastq_dup_filter;
    }
    
    print "Finished duplicate, length, and quality trimming.\n";
    
#########   Use BWA to map to contaminant sequences
    
    if ($m) {
        chdir('/shared/sgbiofuels/contamination/aln_output/microbe/');
        system("bwa aln -t $threads MICROBE $fastq_mid | bwa samse MICROBE - $fastq_mid > $microbe_sam");
        print "Finished mapping to microbial contaminants.\n";
    }
    if ($c) {
        chdir('/shared/sgbiofuels/contamination/aln_output/jatcp/');
        system("bwa aln -t $threads JATCP $fastq_mid | bwa samse JATCP - $fastq_mid > $jatcp_sam");
        print "Finished mapping to Jatropha chloroplast genome.\n";
    }
    if ($o) {
        chdir('/shared/sgbiofuels/contamination/aln_output/jatmt/');
        system("bwa aln -t $threads JATMT $fastq_mid | bwa samse JATMT - $fastq_mid > $jatmt_sam");
        print "Finished mapping to Jatropha mitochondrial genome.\n";
    }
    
#########   Remove reads that aligned to any of the contaminant sequences
    
    chdir('/shared/sgbiofuels/contamination/tmp/');
    system("mv $fastq_mid /shared/sgbiofuels/contamination/tmp/");
    my $input_file = "/shared/sgbiofuels/contamination/tmp/$base" . "_dup_length_qual_filt.fastq";
    if ($m) {
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $microbe_sam $input_file");
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against microbial contaminants.\n";
    }
    if ($c) {
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $jatcp_sam $input_file");
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against Jatropha chloroplast genome.\n";
    }
    if ($o) {
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $jatmt_sam $input_file");
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against Jatropha mitochondrial genome.\n";
    }
    my $final_unmapped = $input_file;
    print "Filtering complete.\n";
    print "Num_reads\tFile\n";
    print `wc -l $fastq | awk '{print \$1/4"\t"\$2}'`;
    print `wc -l $base* | awk '{print \$1/4"\t"\$2}'`;
    my $final_file = $fastq;
    $final_file =~ s/\.fastq/_dup_length_qual_contam_filt\.fastq/;
    system("mv $final_unmapped $final_file");
}

if ($pe) {
    my $fastq1;
    my $fastq2;
    my $dir = getcwd;
    if (($ARGV[-1] =~ /\.fastq/) and ($ARGV[-2] =~ /\.fastq/)) {
        $fastq1 = $ARGV[-2];
        $fastq2 = $ARGV[-1];
    }
    else {
        system("perl /shared/sgbiofuels/contamination/bin/split_fastq.pl $ARGV[-1]");
        $fastq1 = $ARGV[-1];
        $fastq2 = $ARGV[-1];
        die if ($ARGV[-1] !~ /\.fastq/);
        $fastq1 =~ s/\.fastq/\_A.fastq/;
        $fastq2 =~ s/\.fastq/\_B.fastq/;
    }
    if ($fastq1 !~ /^\//) {
        $fastq1 = $dir . '/' . $fastq1;
    }
    if ($fastq2 !~ /^\//) {
        $fastq2 = $dir . '/' . $fastq2;
    }
    if (($fastq1 !~ /\.fastq/) or ($fastq2 !~ /\.fastq/)) {
        die "\nInvalid fastq file!\n$usage";
    }  
    my $fastq1_mida = $fastq1;
    my $fastq1_mid = $fastq1;
    my $fastq_singleton = $fastq1;
    my @fastq1_split = split(/\//, $fastq1);
    my $base1 = $fastq1_split[-1];
    $base1 =~ s/\.fastq//;
    $fastq1_mida =~ s/\.fastq/_dup_length_qual_filt\.fastq/;
    $fastq1_mid =~ s/\.fastq/_dup_length_qual_singleton_filt\.fastq/;
    my $fastq1_dup_filter = $fastq1;
    $fastq1_dup_filter =~ s/\.fastq/_dup_filt\.fastq/;
    my $fastq2_dup_filter = $fastq2;
    $fastq2_dup_filter =~ s/\.fastq/_dup_filt\.fastq/;
    my $fastq2_mida = $fastq2;
    my $fastq2_mid = $fastq2;
    my @fastq2_split = split(/\//, $fastq2);
    my $base2 = $fastq2_split[-1];
    $base2 =~ s/\.fastq//;
    $fastq2_mida =~ s/\.fastq/_dup_length_qual_filt\.fastq/;
    $fastq2_mid =~ s/\.fastq/_dup_length_qual_singleton_filt\.fastq/;
    $fastq_singleton =~ s/\.fastq/_dup_length_qual_filt_singleton\.fastq/;
    my $hg19_sam = '/shared/sgbiofuels/contamination/aln_output/hg19/bwa_' . $base1 . '_HG19.sam';
    my $microbe_sam = '/shared/sgbiofuels/contamination/aln_output/microbe/bwa_' . $base1 . '_MICROBE.sam';
    my $mikeslist_sam = '/shared/sgbiofuels/contamination/aln_output/mikeslist/bwa_' . $base1 . '_MIKESLIST.sam';
    my $jatcp_sam = '/shared/sgbiofuels/contamination/aln_output/jatcp/bwa_' . $base1 . '_JATCP.sam';
    my $jatmt_sam = '/shared/sgbiofuels/contamination/aln_output/jatmt/bwa_' . $base1 . '_JATMT.sam';
    
    my $hg19_sam_singleton = '/shared/sgbiofuels/contamination/aln_output/hg19/bwa_' . $base1 . '_singleton_HG19.sam';
    my $microbe_sam_singleton = '/shared/sgbiofuels/contamination/aln_output/microbe/bwa_' . $base1 . '_singleton_MICROBE.sam';
    my $mikeslist_sam_singleton = '/shared/sgbiofuels/contamination/aln_output/mikeslist/bwa_' . $base1 . '_singleton_MIKESLIST.sam';
    my $jatcp_sam_singleton = '/shared/sgbiofuels/contamination/aln_output/jatcp/bwa_' . $base1 . '_singleton_JATCP.sam';
    my $jatmt_sam_singleton = '/shared/sgbiofuels/contamination/aln_output/jatmt/bwa_' . $base1 . '_singleton_JATMT.sam';
    
#########   Duplicate, length, and quality filter     
    
    print "Beginning DNA filter of $fastq1 and $fastq2 with $threads threads.\n";
    if ($d) {
        system("python /shared/sgbiofuels/contamination/bin/fastq_unique.py < $fastq1 > $fastq1_dup_filter");
        system("python /shared/sgbiofuels/contamination/bin/fastq_unique.py < $fastq2 > $fastq2_dup_filter");
        print "Finished duplicate filtering.\n";
    }
    else {
        $fastq1_dup_filter = $fastq1;
        $fastq2_dup_filter = $fastq2;
    }
    if ($q) {
        system("fastq_quality_trimmer -t $min_qual -l $min_length -i $fastq1_dup_filter | fastq_quality_filter -q $min_qual -p $pct > $fastq1_mida");
        system("fastq_quality_trimmer -t $min_qual -l $min_length -i $fastq2_dup_filter | fastq_quality_filter -q $min_qual -p $pct > $fastq2_mida");
        print "Finished length and quality filtering/trimming\n";
    }
    else {
        $fastq1_mida = $fastq1_dup_filter;
        $fastq2_mida = $fastq2_dup_filter;
    }

    print "Finished duplicate, length, and quality trimming.\n";
    
    if (($d) or ($q)) {
        print "Removing singletons resulting from filtering\n";
        system("perl /shared/sgbiofuels/contamination/bin/remove_singletons.pl $fastq1_mida $fastq2_mida $fastq1_mid $fastq2_mid $fastq_singleton");
    }
    else {
        system("cp $fastq1 $fastq1_mid");
        system("cp $fastq2 $fastq2_mid");
    }
    
#########   Use BWA to map to contaminant sequences
    
    if ($m) {
        chdir('/shared/sgbiofuels/contamination/aln_output/microbe/');
        system("bash -c 'bwa sampe MICROBE <(bwa aln -t $threads MICROBE $fastq1_mid) <(bwa aln -t $threads MICROBE $fastq2_mid) $fastq1_mid $fastq2_mid > $microbe_sam'");
        system("bwa aln -t $threads MICROBE $fastq_singleton | bwa samse MICROBE - $fastq_singleton > $microbe_sam_singleton");
        print "Finished mapping to microbial contaminants.\n";
    }
    if ($c) {
        chdir('/shared/sgbiofuels/contamination/aln_output/jatcp/');
        system("bash -c 'bwa sampe JATCP <(bwa aln -t $threads JATCP $fastq1_mid) <(bwa aln -t $threads JATCP $fastq2_mid) $fastq1_mid $fastq2_mid > $jatcp_sam'");
        system("bwa aln -t $threads JATCP $fastq_singleton | bwa samse JATCP - $fastq_singleton > $jatcp_sam_singleton");
        print "Finished mapping to Jatropha chloroplast genome.\n";
    }
    if ($o) {
        chdir('/shared/sgbiofuels/contamination/aln_output/jatmt/');
        system("bash -c 'bwa sampe JATMT <(bwa aln -t $threads JATMT $fastq1_mid) <(bwa aln -t $threads JATMT $fastq2_mid) $fastq1_mid $fastq2_mid > $jatmt_sam'");
        system("bwa aln -t $threads JATMT $fastq_singleton | bwa samse JATMT - $fastq_singleton > $jatmt_sam_singleton");
        print "Finished mapping to Jatropha mitochondrial genome.\n";
    }
    
    
#########   Remove reads that aligned to any of the contaminant sequences
    
    chdir('/shared/sgbiofuels/contamination/tmp/');
    system("mv $fastq1_mida /shared/sgbiofuels/contamination/tmp/");
    system("mv $fastq2_mida /shared/sgbiofuels/contamination/tmp/");
    system("mv $fastq1_mid /shared/sgbiofuels/contamination/tmp/");
    system("mv $fastq2_mid /shared/sgbiofuels/contamination/tmp/");
    system("mv $fastq_singleton /shared/sgbiofuels/contamination/tmp/");
    my $input_file1 = "/shared/sgbiofuels/contamination/tmp/$base1" . "_dup_length_qual_singleton_filt.fastq";
    my $input_file2 = "/shared/sgbiofuels/contamination/tmp/$base2" . "_dup_length_qual_singleton_filt.fastq";
    my $input_file = "/shared/sgbiofuels/contamination/tmp/$base1" . "_dup_length_qual_filt_singleton.fastq";
    if ($m) {
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads.pl $microbe_sam $input_file1 $input_file2");
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $microbe_sam_singleton $input_file");
        $input_file1 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file2 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against microbial contaminants.\n";
    }
    if ($c) {
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads.pl $jatcp_sam $input_file1 $input_file2");
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $jatcp_sam_singleton $input_file");
        $input_file1 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file2 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against Jatropha chloroplast genome.\n";
    }
    if ($o) {
        print "perl /shared/sgbiofuels/contamination/bin/split_mapped_reads.pl $jatmt_sam $input_file1 $input_file2\n";
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads.pl $jatmt_sam $input_file1 $input_file2");
        system("perl /shared/sgbiofuels/contamination/bin/split_mapped_reads_single.pl $jatmt_sam_singleton $input_file");
        $input_file1 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file2 =~ s/\.fastq/_unmapped\.fastq/;
        $input_file =~ s/\.fastq/_unmapped\.fastq/;
        print "Finished filtering against Jatropha mitochondrial genome.\n";
    }

    my $final_unmapped1 = $input_file1;
    my $final_unmapped2 = $input_file2;
    my $final_unmapped = $input_file;
    print "Filtering complete.\n";
    print "Num_reads\tFile\n";
    print `wc -l $fastq1 | awk '{print \$1/4"\t"\$2}'`;
    print `wc -l $base1* | awk '{print \$1/4"\t"\$2}'`;
    print `wc -l $fastq2 | awk '{print \$1/4"\t"\$2}'`;
    print `wc -l $base2* | awk '{print \$1/4"\t"\$2}'`;
    my $final_file1 = $fastq1;
    $final_file1 =~ s/\.fastq/_dup_length_qual_contam_filt\.fastq/;
    system("mv $final_unmapped1 $final_file1");
    my $final_file2 = $fastq2;
    $final_file2 =~ s/\.fastq/_dup_length_qual_contam_filt\.fastq/;
    system("mv $final_unmapped2 $final_file2");
    my $final_file = $fastq1;
    $final_file =~ s/\.fastq/_dup_length_qual_contam_filt_singleton\.fastq/;
    system("mv $final_unmapped $final_file");
}
