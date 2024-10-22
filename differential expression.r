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

