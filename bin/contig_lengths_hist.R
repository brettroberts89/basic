#TITLE: contig_length_hist.R
#AUTHOR: Brett Roberts
#DATE LAST MODIFIED: 11/07/12
#PURPOSE: creates a histogram showing the distribution of contig and scaffold lengths in the two assemblies

input_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/sspace/"
output_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/sspace/"

setwd(input_dir)
gap_lengths <-read.csv("run5_gaplengths.txt", header=TRUE)

setwd(output_dir)

hist(gap_lengths$length,breaks=seq(from=0,to=10000, by=500),col="grey", xlab="Gap Length", main="SSPACE Run5 Gap Lengths",ylim=c(1,900))


