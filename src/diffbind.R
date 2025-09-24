library(DiffBind)
library(tidyverse)
library(CancerSubtypes)

# read in peaksets
# column PeakCaller has been changed from narrow to bed because score is regarded as the 8th column for narrow and 5th column for bed, see pv.defaultScoreCol from https://rdrr.io/bioc/DiffBind/src/R/helper.R
samples <- read.csv('meta.csv')
dbObj <- dba(sampleSheet=samples)

# dbObj$peaks contains sample wise peaks
dbObj$totalMerged
dim(dbObj$merged)
dim(dbObj$allcalled)
dim(dbObj$called)
dim(dbObj$binding)
dbObj

# affinity binding matrix
# compute count for each of the peaks in the consensus set (appears in more than two of the samples)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits = FALSE)

dim(dbObj$called)
dim(dbObj$binding)
dim(dbObj$merged)
dbObj$totalMerged
dbObj

# exploratory data analysis
# see how well samples cluster with one another
pdf("dba.plotPCA.pdf")
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)
dev.off()

pdf("corr.heatmap.pdf")
plot(dbObj)
dev.off()

# output affinity matrix (TMM normalized as default in dba.count() step)
dba.peakset(dbObj, bRetrieve=TRUE, writeFile="BindingAffinityMatrix.csv")

# find the top 10,000 most variable features
data = read.csv('BindingAffinityMatrix.csv', sep="\t", header=FALSE)
rownames(data) = paste(data$V1,data$V2,data$V3,sep="-")
data = data.matrix(data[1:nrow(data),4:ncol(data)])
colnames(data) =  dbObj$samples$SampleID
dim(data)
# FSbyMAD from CancerSubtypes package, MAD: median absolute deviation
# for each row (feature), MAD is calulated by: ‘constant * cMedian(abs(x - center))’, where constant: 1.4826, so 1.4826 x median(|value-median(vector)|)
data1 = FSbyMAD(data, cut.type="topk", value=10000)
write.table(data1,file="mat.txt")

# keep non +/-1kbTSS region overlapped peeks
data = read.csv('BindingAffinityMatrix.csv', sep="\t", header=FALSE)
rownames(data) = paste(data$V1,data$V2,data$V3,sep="-")
data = data.matrix(data[1:nrow(data),4:ncol(data)])
dim(data)

annot = read.table("consensus_peaks.final.bed", header = F)
annot = annot[annot$V4 != 4,]
nonTSS = paste(annot$V1, annot$V2, annot$V3,sep="-")
data = data[rownames(data) %in% nonTSS,]
write.table(data,file="BindingAffinityMatrix.csv")
data1 = FSbyMAD(data, cut.type="topk", value=10000)
write.table(data1,file="mat.txt")
