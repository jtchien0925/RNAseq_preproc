library("DESeq2")

#tuto2
directory<-"your dir"

#combine all files wiht "m" intiate
sampleFiles<-grep("your dataset name, e.g. Rzz<--you can use R in here",list.files(directory),value=TRUE)

#Given the batch
#e.g. -->sampleBatch <- c("TP3","TP6","TP3","TP6","TP3","TP6","TP3","TP6")
sampleBatch <- c("")
#Given the condition, in our exp is by monkey
#e.g. -->sampleCondition<-c("TP3","TP6","TP3","TP6","TP3","TP6","TP3","TP6")
sampleCondition<-c("")


sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition, Batch = sampleBatch)
#you can view sampleTable here
#sampleTable

#tell DEseq that data is from HT-seq
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)

#you can check ddsHTseq
#ddsHTseq

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("A","B")) # in here its youe comparison sets A and B

#ddsHTseq
#tuto2
dds<-DESeq(ddsHTSeq)
res<-results(dds)
# order in p-value adjustment
res<-res[order(res$padj),]
head(res)


#build a data frame that contains the results of DEseq
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="DATE-DESeq2-results-with-normalized.csv")

#take out outliners
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)

write.csv(as.data.frame(tab),file = 'DATE-DESeq2-replaceoutliers.csv')

#make a good looking table
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file = 'DATE-DESeq2-replaceoutliers-results.csv')

#upregulating
write.csv(as.data.frame(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ]),file = 'DATE-DESeq2upregulat.csv')
#downregulating
write.csv(as.data.frame(resSig[ order(resSig$log2FoldChange), ]),file = 'DATE-DESeq2downregulat.csv')

#find the significant expression
sum(resClean$padj < 0.1, na.rm=TRUE)
resSig <- subset(resClean, padj < 0.1)

#MA plot; goole MA plot
plotMA(dds,ylim=c(-10,10),main='DESeq2')

#check which is significant expressed
#using falseDiscoveryRate
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
#LogFoldChange
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)



#LFC csv
write.csv(as.data.frame(resLFC1),file = 'DATE-DESeq2-resLFC1.csv')




#topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
#with(resLFC1[topGene, ], {
#  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
#  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
#})

#topGene <- rownames(resLFC1)[which(resLFC1$padj<0.1)]
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=1)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
dev.copy(png,'deseq2_MAplot.png')
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)


# save normalized values
write.table(as.data.frame(assay(rld)),file='DATE-DESeq2-rlog-transformed-counts.csv', sep=',')
write.table(as.data.frame(assay(vsd)),file='DATE-DESeq2-vsd-transformed-counts.csv', sep=',')


#gene clustering, top variance genes, rLog trnasformation
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld))
heatmap(mat, annotation_col=df)

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"DATE-DESeq2_variance_stabilizing.png")
dev.off()



# clustering analysis
library("gplots")
library("RColorBrewer") 
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition,sampleFiles,sep=" : "))
hc <- hclust(distsRL)
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.copy(png,'deseq2_heatmaps_samplebysample.png')
dev.off()
#heatmap.2(mat, trace = "none", col = rev(vstcol), margin = c(10, 10))
#dev.copy(png, "DATE-DESeq2-clustering.png")
#dev.off()

# Principal components plot
# will show additional clustering of samples
# showing basic PCA function in R from DESeq2 package
# this lacks sample IDs and only broad sense of sample clustering
# its not nice - but it does the job
#print(plotPCA(rld, intgroup = c("Batch")))
#dev.copy(png, "DATE-DESeq2_PCA_initial_analysis.png")
#dev.off()


# or ggplot PCA plot
library("grDevices")
library('ggplot2')
library("genefilter")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vsd)[select,]))

condition <- c("")
#"Baseline","TP2","Baseline","TP2","Baseline","TP2","Baseline","TP2","Baseline","TP2"


scores <- data.frame(sampleFiles, pca$x, condition)

#better looking PCA plot, I have to tune some paramters to describe single dots
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
 + geom_point(size = 3)
 + ggtitle("Principal Components")
 + scale_colour_brewer(name = " ", palette = "Set1")
 + theme(
   plot.title = element_text(face = 'bold'),
   legend.key = element_rect(fill = 'NA'),
   legend.text = element_text(size = 3, face = "bold"),
   axis.text.y = element_text(colour = "Black"),
   axis.text.x = element_text(colour = "Black"),
   axis.title.x = element_text(face = "bold"),
   axis.title.y = element_text(face = 'bold'),
   panel.grid.major.x = element_blank(),
   panel.grid.major.y = element_blank(),
   panel.grid.minor.x = element_blank(),
   panel.grid.minor.y = element_blank(),
   panel.background = element_rect(color = 'black',fill = NA)
 ))



# heatmap of data
library("RColorBrewer")
library("gplots")
# 100 top expressed genes with heatmap.2, the number can be adjusted; see [1:100] part
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:100]

#making a heatmap
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Top100-DE-HeatMap")
dev.copy(png, "DATE-DESeq2-HEATMAP.png")
dev.off()


#### top 2000 genes based on row variance with heatmap3 based on VSD ###
library(heatmap3)
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=T)[seq_len(min(2000,length(rv)))]
my_palette <- colorRampPalette(c("blue", "white", "red"))(1024)
heatmap3(assay(vsd)[select,],col=my_palette,
         labRow = F,cexCol = 0.8,margins=c(6,6))
dev.copy(pdf, "DATE-DESeq2-heatmap3.pdf")
dev.off()




# I still want to make a heatmap here..
library('RColorBrewer')
library('gplots')
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10,10))
dev.copy(png,'DESeq2_heatmap1_raw counts')
dev.off()
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 10))

dev.copy(png,'DESeq2_heatmap2_rlog')
dev.off()


heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 10))
dev.copy(png,'DESeq2_heatmap3_vsd')
dev.off()




#########################################################
####needed to check
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj)
addmargins(tab)
res<-results(dds)
resClean <- results(ddsClean)
plotDispEsts(dds)

attr(resClean,'10')
plot(attr(res,'filterNumRej'),type='b', ylab='number of rejections')
dev.copy(png,'deseq2_filtering_treshold.png')
dev.off()
#?


W <- res$stat
maxCooks <- apply(assays(dds)[['cooks']],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab='rank of Wald statistic',
     ylab='maximum Cooks distance per gene',ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
dev.copy(png,'deseq2_cooksdist.png')
dev.off()




plot(res$baseMean+1, -log10(res$pvalue),
     log='x', xlab='mean of normalized counts',
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
dev.copy(png,'deseq2_indep_filt.png')
dev.off()


