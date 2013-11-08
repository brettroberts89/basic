#TITLE: fix_fastq.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 2/4/2013
#PURPOSE: add @ symbol to beginning of read.

my $fastq = $ARGV[0];

open FASTQ, "<$fastq";
my $alt = $fastq;
$alt =~ s/\.fastq/_fixed\.fastq/;
open OUT, ">$alt";
while (my $line_1 = <FASTQ>) {
    chomp $line_1;
    my $seq;
    my $plus;
    my $quality;
    while (my $line_2 = <FASTQ>) {
        chomp $line_2;
        $seq = $line_2;
        while (my $line_3 = <FASTQ>) {
            chomp $line_3;
            $plus = $line_3;
            while (my $line_4 = <FASTQ>) {
                chomp $line_4;
                $quality = $line_4;
                last;
            }
            last;
        }
        last;
    }
    next if ($line_1 =~ /^\@/);
    $line_1 = '@' . $line_1;
    print OUT "$line_1\n$seq\n$plus\n$quality\n";
}