ensembl_gene_id <- row.names(diff_gene)
rat_symbols <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"), filters='ensembl_gene_id', values=ensembl_gene_id, mart=mart)
diff_gene_ensembl_id <- rat_symbols$ensembl_gene_id

#GO
for(i in c("MF", "BP", "CC")){
    ego <- enrichGO(gene       = rat_symbols$entrezgene_id,
                    OrgDb      = org.Rn.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = i,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05)
    dotplot(ego, showCategory = 30, title = paste("The GO ", i, " enrichment analysis", sep = ""))
}

#KEGG
kk <- enrichKEGG(gene = gene, 
                 organism ='rno',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data = FALSE)