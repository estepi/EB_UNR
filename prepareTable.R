library(edgeR)
setwd("~/Documents/UNR/clases/data/mouse/STAR")
###########################################
files<-dir(pattern="*\\.tab")
RG <- readDGE(files, header=F, skip=4) # be careful
write.table(RG$counts,"geneCountsMmu.tab", sep="\t", quote=F)
###########################################

counts<-counts0[,-1]
head(counts)
colnames(counts)<-gsub("_1", "",colnames(counts))

targets<-read.table("targets.txt", sep = "\t", header = F)
targets
colnames(counts)
ii<-match(targets$V1, colnames(counts))

group <- factor(targets$V2)
head(counts)
dim(counts)
head(group)
length(group)

delete.na <- function(df, n=0) {
  df[rowSums(is.na(df)) <= n,]
}

countsClea<-delete.na(counts)
dim(counts)
dim(countsClea)
y <- DGEList(counts=countsClea,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

png("MDS.png")
plotMDS(y, labels = group)
dev.off()
###########################################
