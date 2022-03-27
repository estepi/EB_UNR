## ----setup, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = FALSE)


## ----targets, echo=TRUE, eval=TRUE-----------------------------------
targets <- read.table("~/Documents/UNR/clases/data/HNRPC/targets.txt", header = T, row.names = 1)


## ----counts, echo=TRUE, eval=TRUE------------------------------------
counts <-  read.table( "~/Documents/UNR/clases/data/HNRPC/genes.hits.txt",   
                       row.names = 1, header = T, stringsAsFactors = F)


## ----descriptive, echo=TRUE, eval=TRUE-------------------------------
class(counts)
dim(counts)
head(counts)


## ----DESeq, echo=TRUE, eval=TRUE, message=FALSE----------------------
library(DESeq2)

coldata <- data.frame(sample=rownames(targets),
                      condition=targets$condition)
coldata$condition <- factor(coldata$condition)
coldata$type <- rep("SE")
condition <- factor(targets$condition)


## ----DESeqI, echo=TRUE, eval=TRUE------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              design   = ~condition,
                              colData = coldata)
dds


## ----DE, echo=TRUE, eval=TRUE----------------------------------------
dds <- DESeq(dds)
resultsNames(dds) 
res <- results(dds)
res


## ----orderDE, echo=TRUE, eval=TRUE-----------------------------------
resOrdered <- res[order(res$pvalue),]


## ----orderDE2, echo=TRUE, eval=TRUE----------------------------------
y <-as.data.frame(res)
resOrdered <- data.frame(res[order(abs(res$padj)),][1:20,])
head(resOrdered)


## ----vulcanoDE, echo=TRUE, eval=TRUE, warning=FALSE, message=FALSE----
library(ggplot2)
library(ggrepel)
y$gene_color <- rep("grey", nrow(y))
y$gene_color[y$log2FoldChange>1] <-"red"
y$gene_color[y$log2FoldChange< (-1)]<-"green"
y$imp_genes<-NA
ii <- match(rownames(resOrdered), rownames(y))
y$imp_genes[ii]<-rownames(y)[ii]

ggplot(y, aes(x=log2FoldChange,
              y=-log10(padj))) +
  geom_point(aes(col=gene_color), cex= 1.2) +
  scale_color_manual(values=c("dark green","dark grey", "dark red")) +
  labs(title="DEG", x="log2(FC)", y="-log10(FDR)") +
  geom_vline(xintercept= c(-1, 1), colour= 'black', linetype= 'dashed') +
  geom_hline(yintercept= 1.30103, colour= 'black', linetype= 'dashed') +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face="italic", hjust=0.4),
        axis.title.x = element_text(color = "black", size=12, hjust = 0.4),
        axis.title.y = element_text(size =12, hjust = 0.5)) +
  geom_text_repel(data=y,
                  aes(x=log2FoldChange, y=-log10(padj)),
                  label =y$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)


## ----tt2, echo=TRUE, eval=FALSE--------------------------------------
## write.table() o write.csv():
## write.csv(res, file="res_DESeq2.csv")

