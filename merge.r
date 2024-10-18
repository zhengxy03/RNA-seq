rm(list=ls())

setwd("//wsl.localhost/Ubuntu/home/zxy0303/project")

files <- list.files(".", "*.count")
f_lists <- list()
for (i in files){
    prefix = gsub("(_\\w+)?\.count", "", i, perl=TRUE)
    f_lists[[prefix]]=i
}
print(f_lists)