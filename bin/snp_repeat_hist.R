#TITLE: snp_repeat_hist.R
#AUTHOR: Brett Roberts
#DATE LAST MODIFIED: 08/20/12
#PURPOSE: creates a histogram showing the frequency of different types of repeats in the sequences 50 bp to either side of snps

input_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/snp_data/"
output_dir="/Users/robotien/Desktop/Brett/SG_Biofuels/dilara/snp_data/"

setwd(input_dir)
snps <-read.csv("snp_hist.csv", header=TRUE)
repeat_pct <-snps_repeats <-read.csv("snp_repeat_pct_hist.csv", header=TRUE,colClasses=c(rep("numeric",4)))

setwd(output_dir)

#create bar plot showing percentage of snps that have unique regions vs. repetitive regions
bp1 <- barplot(snps$Percentage,col="grey",xlab="Characteristics of Region +-100 bp of SNP",ylab="Proportion of SNPs",ylim=c(0,0.6),names.arg= c("Unique","Multiple Types of Repeats","Mononucleotide\n(>=5 repeats)","Dinucleotide\n(>=5 repeats)","Trinucleotide\n(>=5 repeats)"))

#create plot of average percent of region that is repetitive on both the left and right

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
	}


avgs = apply(snps_repeats, 2, mean)
sd = apply(snps_repeats, 2, sd)

 left_repeat <- avgs[3]
 right_repeat <- avgs[4]
 left_unique = 1-left_repeat
 right_unique = 1-right_repeat
 
 plot_values = c(left_repeat,left_unique,right_repeat,right_unique)
 sd_values = c(sd[3],sd[3],sd[4],sd[4])
 m1 <- matrix(plot_values,2)
 m2 <- matrix(sd_values,2)
 upper_arrow = m1+m2
 lower_arrow = m1-m2
 
 bp2 <- barplot(as.matrix(m1),ylab="Percentage of Sequence",col=rainbow(2),ylim = c(0,1),names.arg=c("Left","Right"),beside=TRUE)
 legend("topleft", c("Repetitive","Unique"), cex=0.6,bty="n",fill=rainbow(2))
error.bar(bp2,m1,m2)

 
plot(c(1,2,3,4),plot_values,xlab="",ylab="Percentage of Sequence",ylim=c(0,1),axes=FALSE)
axis(1,at=1:4, lab=c("Left Repetitive","Left Unique","Right Repetitive","Right Unique"),pos = -0.1)
axis(2, las=1,at=0.2*0:5)
arrows(c(1,2,3,4),lower_arrow,c(1,2,3,4),upper_arrow,length=.1,angle=90,code=3)







