setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/annotation")
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("rn7.gff")

exons_gene <- exonsBy(txdb, by = "gene")

gene_len <- list()
for (i in names(exons_gene)){
    range_info = reduce(exons_gene[[i]])
    width_info = width(range_info)
    sum_len =sum(width_info)
    gene_len[[i]] = sum_len
}

gene_len <- lapply(exons_gene, function(x){sum(width(reduce(x)))})

data <- t(as.data.frame((gene_len)))

setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output/HTseq")
write.table(data, file="rn7_gene_len.tsv", row.names = TRUE, sep="\t", col.names = FALSE)

#RPKM
setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output/HTseq")
gene_len_file <- "rn7_gene_len.tsv" #annotation_exons_length
gene_len <- read.table(gene_len_file, header = FALSE, row.name = 1)
colnames(gene_len) <- c("length")

count_file <- "SRR2190795.count" #HEseq
count <- read.table(count_file, header = FALSE, row.name = 1)
colnames(count) <-c("count")
all_count <- sum(count["count"])

RPKM <- c()
for (i in row.names(count)){
    count_ = 0
    exon_kb = 1
    rpkm = 0
    count_ = count[i, ]
    exon_kb = gene_len[i, ] / 1000
    rpkm = (10 ^ 6 * count_ ) / (exon_kb * all_count )
    RPKM = c(RPKM, rpkm)
}
count["RPKM"] <- RPKM