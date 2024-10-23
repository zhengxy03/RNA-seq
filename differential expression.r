#pre-treatment
setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output/HTseq")
dataframe <- read.csv("merge.csv", header=TRUE, row.names =1)
countdata <- dataframe[-(1:5), ]
head(countdata)

row_names <- row.names(countdata)
new_row_names <- gsub("\\.\\w+","", row.names(countdata))
row.names(countdata) <- new_row_names

countdata <- countdata[rowSums(countdata) > 0,]

#download
# 使用bioconductor进行安装
source("http://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

# 安装包
biocLite("DESeq2")
biocLite("pheatmap")
biocLite("biomaRt")
biocLite("org.Rn.eg.db")
biocLite("clusterProfiler")

# 加载
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)


#create dds
setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output")
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )

head(coldata)
countdata <- countdata[row.names(coldata)]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds

#sample correlation(plot)
vsdata <- rlog(dds, blind=FALSE)
plotPCA(vsdata, intgroup="treatment") + ylim(-10, -10)

library("RColorBrewer")
gene_data_transform <- assay(vsdata)
sampleDists <- dist(t(gene_data_transform))
sampleDistsMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
# colnames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blue")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#results
dds$treatment <- factor(as.vector(dds$treatment), levels = c("control","treatment"))
dds <- DESeq(dds)

result <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
head(result)
result_order <- result[order(result$pvalue),]
head(result_order)

dir.create("../DESeq2") #output
write.csv(result, file="../DESeq2/results.csv", quote = F)

#analysis results
summary(result_order)

table(result_order$padj<0.05)

#different genes
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)
write.csv(diff_gene, file="../DESeq2/difference.csv", quote = F)

plotMA(result_order, ylim=c(-10,10))
