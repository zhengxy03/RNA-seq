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

write.table(data, file="rn7_gene_len.tsv", row.names = TRUE, sep="\t", col.names = FALSE)