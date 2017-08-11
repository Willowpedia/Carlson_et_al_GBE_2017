
require(ggplot2)
require(reshape2)

args <- commandArgs(TRUE)

data_file <- args[1]
image_file <- args[2]
out_file <- args[3]
ploidy_p1 <- args[4]
ploidy_p2 <- args[5]

data_file
image_file

data <- read.delim(data_file, header=TRUE, sep="\t")

FDR_cutoff <- 0.005

for (i in 1:nrow(data)) {
	if (data$parents.binom.qvalue[i]<=FDR_cutoff) {
		if (data$hybrids.binom.qvalue[i]<=FDR_cutoff) {
			if (data$fisher.exact.qvalue[i]<=FDR_cutoff) {
				ratio_in_p <- data$CPM_P1_IN_P1[i]/data$CPM_P2_IN_P2[i]
				ratio_in_h <- (data$CPM_P1_IN_H[i]/as.numeric(ploidy_p1))/(data$CPM_P2_IN_H[i]/as.numeric(ploidy_p2[1]))
				if (log2(ratio_in_p)/log2(ratio_in_h) > 1) {
					data$regulatory[i] <- "cis + trans"
				} else {
					data$regulatory[i] <- "cis x trans"
				}
			} else {
				data$regulatory[i] <- "cis only"
			}
		} else {
			if (data$fisher.exact.qvalue[i]<=FDR_cutoff) {
				data$regulatory[i] <- "trans only"
			} else {
				data$regulatory[i] <- "ambiguous"
			}
		}
	} else {
		if (data$hybrids.binom.qvalue[i]<=FDR_cutoff) {
			if (data$fisher.exact.qvalue[i]<=FDR_cutoff) {
				data$regulatory[i] <- "compensatory"
			} else {
				data$regulatory[i] <- "ambiguous"
			}
		} else {
			if (data$fisher.exact.qvalue[i]<=FDR_cutoff) {
				data$regulatory[i] <- "ambiguous"
			} else {
				data$regulatory[i] <- "conserved"
			}
		}
	}
}

#plot_data <- data[c(CPM_P1_IN_P1, CPM_P2_IN_P2, CPM_P1_IN_H, CPM_P2_IN_H, regulatory)]

data$log_ratio_p <- log2(data$CPM_P1_IN_P1/data$CPM_P2_IN_P2)
data$log_ratio_h <- log2(data$CPM_P1_IN_H/data$CPM_P2_IN_H)

png(paste0(image_file, ".png"), width=4000, height=4000, res=300)
#pdf(paste0(image_file, ".pdf"))

q <- ggplot()
size <- 2
alpha <- 0.5

group.colors <- c("cis only" = "black", "trans only" = "red", "ambiguous" = "grey", "conserved" = "yellow", "compensatory" = "orange", "cis + trans" = "purple", "cis x trans" = "green") 

group.colors

q <- q +
	geom_point(data=data, aes(x=log_ratio_p, y=log_ratio_h, group=regulatory, colour=regulatory), size=size, alpha=alpha) +
	ylim(-10,10) + xlim(-10,10) + ggtitle(data_file) +
	xlab("log_2(P1/P2) in parents") + ylab("log_2(P1/P2) in hybrids") +
	scale_colour_manual(values=group.colors) +
	theme(legend.justification=c(1,0), legend.position=c(1,0))

print(q)

#print(
#	ggplot()
#	+ geom_point(data=subset(data, regulatory=="ambiguous"), aes(x=log_ratio_p, y=log_ratio_h), colour="grey", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="conserved"), aes(x=log_ratio_p, y=log_ratio_h), colour="yellow", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="compensatory"), aes(x=log_ratio_p, y=log_ratio_h), colour="orange", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="cis + trans"), aes(x=log_ratio_p, y=log_ratio_h), colour="purple", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="cis x trans"), aes(x=log_ratio_p, y=log_ratio_h), colour="green", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="cis only"), aes(x=log_ratio_p, y=log_ratio_h), colour="black", size=size, alpha=alpha)
#	+ geom_point(data=subset(data, regulatory=="trans only"), aes(x=log_ratio_p, y=log_ratio_h), colour="red", size=size, alpha=alpha)
#	+ ylim(-10,10) + xlim(-10,10) + ggtitle(data_file) 
#	+ xlab("log_2(P1/P2) in parents") + ylab("log_2(P1/P2) in hybrids")
#)

#print(
#	ggplot()
# + geom_point(data=data, aes(x=log_ratio_p, y=log_ratio_h, group=regulatory, color=regulatory), shape=1, size=2)
# + ylim(-10,10) + xlim(-10,10)
# + xlab("log_2(P1/P2) in parents") + ylab("log_2(P1/P2) in hybrids")
#)

dev.off()

write.table(data, out_file, row.names=FALSE, sep="\t");


