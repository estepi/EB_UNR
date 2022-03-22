knitr::opts_chunk$set(echo = FALSE)


library(RNAseqData.HNRNPC.bam.chr14)
?RNAseqData.HNRNPC.bam.chr14
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
bamfiles <- RNAseqData.HNRNPC.bam.chr14_BAMFILES

targets <-
     data.frame(
         row.names = names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
         bam = matrix(
             unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
             ncol = 1,
             byrow = TRUE
),
         condition = c(rep("CT", 4), rep("KD", 4))
     )

targets
write.table(targets, "targets.txt", quote = F, sep = "\t")
 
files = as.character(targets$bam)
reads <- lapply(files,
                 function(x) {
                     aln <- readGAlignments(x)
                 })

names(reads) <- rownames(targets)

getwd()


c14 <- makeTxDbFromGFF("../data/chr14.gtf")
saveDb(c14 ,"chr14.sqlite") 

#c14 <- loadDb("chr14.sqlite")

feature <- exonsBy(c14, by = "gen")
 hits <-
     lapply(reads, function(x) {
         countOverlaps(feature, x, ignore.strand = TRUE)
     })
       
 hits.ul <- do.call(cbind.data.frame, hits)

write.table(hits.ul, file = "genes.hits.txt" )
head(hits.ul) 

targets 
 

 counts <-
   read.table(
         "genes.hits.txt",
       row.names = 1,
       header = T,
         stringsAsFactors = F     )

 class(counts)
 dim(counts)
 head(counts)
 
targets <- read.table("../data/targets.txt", header = T, row.names = 1)

class(counts)
dim(counts)
head(counts)


library(edgeR)


 summary(counts)
 plot(density(rowMeans(counts)))


 par(mfrow=c(2,4))
 apply(counts, 2, boxplot)


cpms<-cpm(counts)
keep<-rowSums(cpms>1)>4


dim(counts)
countsf<-counts[keep,]
dim(countsf)


summary(countsf)
d<-DGEList(counts=countsf, group=targets$condition)
d<-calcNormFactors(d)
d
shortNames<-paste(targets$condition, rep(1:4, 2), sep=".")
targets<-cbind(targets,shortNames)
dev.off()
plotMDS(d, labels=targets$shortNames,
col=c("darkgreen","blue")[factor(targets$condition)])


d<-estimateCommonDisp(d, verbose=TRUE)

d<-estimateTagwiseDisp(d)

plotBCV(d)


de<-exactTest(d, pair=c("CT","KD"))


 tt <- topTags(de, n=nrow(d))


table(tt$table$FDR <0.05)
hist(tt$table$FDR, breaks=100 )

deg<-rownames(tt)[tt$table$FDR <.1]

plotSmear(d, de.tags=deg)
abline(h=c(-1,0,1))

write.csv(tt$table, file="toptags_edgeR.csv")

if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
cds<-  newCountDataSet(countData=counts, conditions=targets$condition)
cds<-estimateSizeFactors(cds)
sizeFactors(cds)
cdsB<-estimateDispersions(cds, method="blind")
vsd<-varianceStabilizingTransformation(cdsB)
p<-plotPCA(vsd)
#
cds<-estimateDispersions(cds)
plotDispEsts(cds)
res<-nbinomTest(cds, "CT","KD")
plotMA(res)
resSig<-res[which(res$padj < 0.1),]
table(res$padj < 0.1)
hist(res$pval, breaks=100)
