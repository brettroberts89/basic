#TITLE: contig_lengths_hist.pl
#AUTHOR: BRETT ROBERTS
#DATE LAST MODIFIED: 03/02/13
#PURPOSE: create csv file with contig names and lengths from fasta file

my $fasta = $ARGV[0];
my $output = $ARGV[1];
open FASTA, "<$fasta";

while (my $line = <FASTA>) {                 
    chomp $line;
    if ($line =~ /^>/) {
        my @line = split(/\s+/, $line);
        $contig = $line[0];
        $contig =~ s/>//;
    }
    else{          
        $sequence = $line;
        $sequence =~ s/\s+//;
        $length = length($sequence);
        $sequences{$contig} += $length;
    }
}

open OUT, ">$output";
print OUT "contig,length\n";
foreach my $contig (keys %sequences) {
    print OUT "$contig,$sequences{$contig}\n";
}