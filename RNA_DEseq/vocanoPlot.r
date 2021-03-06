##volvanoPlot
##inport as txt, remember the first column header shold be 'gene'
res <- read.table(file.choose(), header=TRUE)
head(res)

# Make a basic volcano plot

with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)

with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))

with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(MASS)
library(calibrate)

with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(padj), labs=gene, cex=.4))
