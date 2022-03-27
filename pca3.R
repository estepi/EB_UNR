selector <- read.table("~/Documents/UNR/clases/data/mouse/selector.tsv", 
                       header = F, row.names = 1)

colnames(selector)
counts <-
    read.table(
        "~/Documents/UNR/clases/data/mouse/STAR/geneCountsMmu.tab",
        row.names = 1,
        header = T,
        stringsAsFactors = F
    )

dim(counts)
library(edgeR)
head(counts)

cpms<-cpm(counts)
summary(counts)
plot(density(rowMeans(counts)))


library(ggplot2)
dreads<-density(rowMeans(counts))
rmeans<-data.frame(rMeans=rowMeans(counts))
rmeans$color<-rep("gene", nrow(rmeans))

 ggplot(rmeans, aes(x=rMeans, fill=color)) + 
 geom_density()+
  theme_classic()+
  labs(title="Read density")


library(reshape2) 
df<- melt(counts)
ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")

#hay algun dato raro?
# como podemos identificarlo?

cpms<-cpm(counts)


## ----keep, echo=TRUE, eval=TRUE--------------------------------------
keep<-rowSums(cpms>1)>6


## ----keep2, echo=TRUE, eval=TRUE-------------------------------------
dim(counts)
countsf<-counts[keep,]
dim(countsf)
summary(countsf)
colnames(counts)
condition<-factor(rep(c("ENB", "LNB", "NSC"), each=3))

d<-DGEList(counts=counts, group=condition)
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]
d <- calcNormFactors(d)
d$samples
design <- model.matrix(~condition)
d <- estimateDisp(d,design)
fit <- glmFit(d,design)
lrt2 <- glmLRT(fit,coef=2)
topTags(lrt2)

lrt3 <- glmLRT(fit,coef=3)
topTags(lrt3)

## ----names, echo=TRUE, eval=TRUE-------------------------------------

plotMDS(d, labels=condition,
col=c("darkgreen","blue","red")[factor(condition)])

d<-estimateCommonDisp(d, verbose=TRUE)
d<-estimateTagwiseDisp(d)

plotBCV(d)

deLNB<-exactTest(d, pair=c("ENB","LNB"))
ttLNB <- topTags(deLNB)

deNSC<-exactTest(d, pair=c("ENB","NSC"))
ttNSC <- topTags(deNSC, n = nrow(deNSC))

df <- data.frame(
  exp="NSC",
  fdr=ttNSC$table$FDR)

dim(ttNSC)

ggplot(df, aes(x=fdr)) + 
geom_histogram(color="black", fill="white", binwidth=0.01)+
geom_vline(xintercept=0.05, linetype="dashed")+
theme_classic()+
labs(title="FDR distribution")

deg<-rownames(ttNSC)[ttNSC$table$FDR <.05 &   
                  abs(ttNSC$table$logFC )>1 ]
plotSmear(deNSC, de.tags=deg)
abline(h=c(-1,0,1))

y <-ttNSC$table
tt10 <- topTags(deNSC, n=20)

y$gene_color <- rep("grey", nrow(y))
y$gene_color[y$logFC>1] <-"red"   
y$gene_color[y$logFC< (-1)]<-"green"

y$imp_genes<-NA

ii <- match(rownames(tt10), rownames(y))
y$imp_genes[ii]<-rownames(y)[ii]

library(ggrepel)

ggplot(y, aes(x=logFC, y=-log10(FDR))) +
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
                  aes(x=logFC, y=-log10(FDR)), 
                  label =y$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)
logcpm <- cpm(d, log=TRUE)
#selecttipnar los top 100 de SD


###########################################

library(pheatmap)

fit <- glmFit(d,design)
lrt <- glmLRT(fit,coef=2:3)

fc<-data.frame(LNB=lrt$table$logFC.conditionLNB,
                   NSC=lrt$table$logFC.conditionNSC)
rownames(fc)<-rownames(lrt$table)
range<-apply(fc, 1, range)
range
fcToplot<-fc[sort(range, decreasing = T)[1:50],]
dim(fcToplot)

paletteLength <- 10
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(fcToplot), 0, 
                  length.out=ceiling(paletteLength/2) + 1), 
              seq(max(fcToplot)/paletteLength, 
                   max(fcToplot), 
                  length.out=floor(paletteLength/2)))

pheatmap(fcToplot, 
         fontsize_col= 2,
         fontsize_row=2, 
         myColor, 
         breaks=myBreaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
#######################################
