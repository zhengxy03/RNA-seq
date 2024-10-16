rm(list=ls())

setwd("~/project/rat/output/HTseq")

files <- list_files(".", "*.count")
f_lists <- list()
for (i in files){
    prefix = gsub("(_\\w+)?\.count", "", i, perl=TRUE)
    f_lists[[prefix]]=i
}
print(f_lists)