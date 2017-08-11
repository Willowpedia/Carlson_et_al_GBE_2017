
require(qvalue)

args <- commandArgs(TRUE)

data_file <- args[1]
out_file <- args[2]
ploidy_p1 <- args[3]
ploidy_p2 <- args[4]

data_file
out_file
ploidy_p1
ploidy_p2


p_expect_in_progeny <- as.numeric(ploidy_p1)/(as.numeric(ploidy_p1)+as.numeric(ploidy_p2))

p_expect_in_progeny

data <- read.table(data_file, sep="\t", header=TRUE)

#data<-head(data, n=1000)

for (i in 1:nrow(data)) {
	data$parents.binom.pvalue[i] <- binom.test(round(data$NUM_RNA_READS_P1_IN_P1_NORM[i]),
												round(data$NUM_RNA_READS_P1_IN_P1_NORM[i])+round(data$NUM_RNA_READS_P2_IN_P2_NORM[i]),
												p=0.5)$p.value

	data$hybrids.binom.pvalue[i] <- binom.test(round(data$NUM_RNA_READS_P1_IN_H[i]),
												round(data$NUM_RNA_READS_P1_IN_H[i])+round(data$NUM_RNA_READS_P2_IN_H[i]),
												p=p_expect_in_progeny)$p.value

	table <- matrix(c(round(data$NUM_RNA_READS_P1_IN_P1_NORM[i]),
						round((as.numeric(ploidy_p2)/as.numeric(ploidy_p1))*data$NUM_RNA_READS_P2_IN_P2_NORM[i]),
						round(data$NUM_RNA_READS_P1_IN_H[i]),
						round(data$NUM_RNA_READS_P2_IN_H[i]) ),
						nr=2)
	fisher <- fisher.test(table)
	if (fisher$p.value >=1) {
		data$fisher.exact.pvalue[i] <- 1
	} else {
		data$fisher.exact.pvalue[i] <- fisher$p.value
	}
		
		
	#data$fisher.exact.pvalue[i] <- fisher.test(table)$p.value

}

#write.table(data, out_file, row.names=FALSE, sep="\t");

#p <- data$parents.binom.pvalue
#qobj <- qvalue(p)
#data$parents.binom.qvalue <- qobj$qvalues
#
#p <- data$hybrids.binom.pvalue
#qobj <- qvalue(p)
#data$hybrids.binom.qvalue <- qobj$qvalues
#
#p <- data$fisher.exact.pvalue
#qobj <- qvalue(p)
#data$fisher.exact.qvalue <- qobj$qvalues


data$parents.binom.qvalue <- qvalue(data$parents.binom.pvalue)$qvalues
data$hybrids.binom.qvalue <- qvalue(data$hybrids.binom.pvalue)$qvalues
data$fisher.exact.qvalue <- qvalue(data$fisher.exact.pvalue)$qvalues


#write.table(data, out_file, col.names=NA, row.names=TRUE, sep="\t");
#write.table(data, paste0(out_file,".tmp"), row.names=FALSE, sep="\t");
write.table(data, out_file, row.names=FALSE, sep="\t");

