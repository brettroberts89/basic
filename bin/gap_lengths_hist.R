#TITLE: contig_length_hist.R
#AUTHOR: Brett Roberts
#DATE LAST MODIFIED: 11/07/12
#PURPOSE: creates a histogram showing the distribution of contig and scaffold lengths in the two assemblies

input_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/assembly_gaps/"
output_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/assembly_gaps/"

setwd(input_dir)
SGBv1 <-read.csv("SGBv1_gaplengths.txt", header=TRUE)
SGBv2 <-read.csv("SGBv2_gaplengths.txt", header=TRUE)
SGBv3 <-read.csv("SGBv3_gaplengths.txt", header=TRUE)
SGBv4 <-read.csv("SGBv4_gaplengths.txt", header=TRUE)
SGBv5 <-read.csv("SGBv5_gaplengths.txt", header=TRUE)
SGBv11<-read.csv("SGBv11_gaplengths.txt", header=TRUE)

hist(SGBv1$length,breaks=seq(from=0,to=6600, by=100),col="grey", xlab="Gap Length", main="SGBv1 Gap Lengths",ylim=c(1,10000))

hist(SGBv2$length,breaks=seq(from=0,to=6600, by=100),col="grey", xlab="Gap Length", main="SGBv2 Gap Lengths",ylim=c(1,10000))

hist(SGBv3$length,breaks=seq(from=0,to=6600, by=100),col="grey", xlab="Gap Length", main="SGBv3 Gap Lengths",ylim=c(1,10000))

hist(SGBv4$length,breaks=seq(from=0,to=6600, by=100),col="grey", xlab="Gap Length", main="SGBv4 Gap Lengths",ylim=c(1,10000))

hist(SGBv5$length,breaks=seq(from=0,to=6600, by=100),col="grey", xlab="Gap Length", main="SGBv5 Gap Lengths",ylim=c(1,10000))

hist(SGBv11$length,breaks=seq(from=0,to=4000, by=100),col="grey", xlab="Gap Length", main="SGBv11 Gap Lengths",ylim=c(1,10000))

