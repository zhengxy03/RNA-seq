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
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("biomaRt")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("clusterProfiler")

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
plotPCA(vsdata, intgroup="treatment") + ylim(-10, 10)

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

#biomaRt
#choose database & get symbols
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
ensembl_gene_id <- row.names(diff_gene)
rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)

#merge diff_gene & symbols
diff_gene$ensembl_gene_id <- ensembl_gene_id
diff_gene_dataframe <- as.data.frame(diff_gene)
diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by= c("ensembl_gene_id"))

#save data
write.table(result, "../stat/all_gene.tsv", sep="\t", quote = FALSE)
write.table(diff_gene_symbols, "../stat/diff_gene.tsv", row.names = F,sep="\t", quote = FALSE)

#count diff_gene(bash)

#plot
library(ggplot2)
data <- read.table("all_sampels.tsv", head=T)

pdf("samples_diff_gene_num.pdf")
  ggplot(data=data, aes(x=sample, y=num, fill=samples)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x= "samples", y= "num", title= "different gene number")
dev.off()