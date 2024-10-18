rm(list=ls())

setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output/HTseq")

files <- list.files(".", "*.count")
f_lists <- list()
for (i in files){
    prefix = gsub("(_\\w+)?\.count", "", i, perl=TRUE)
    f_lists[[prefix]]=i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for (i in id_list){
    count <- count + 1
    a <- read.table(f_lists[[i]], sep="\t", col.names=c("gene_id", i))
    data[[count]] <- a
}

data_merge <- data[[1]]
for (i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]], by="gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)