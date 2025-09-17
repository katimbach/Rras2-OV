library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

setwd("/your/parent/dir/")

heatmap_genes<- readRDS("./Results/Tumor_500_variable_genes.rds")
# or read in from "./Input_Files/Tumor_500_variable_genes.rds" if script 10 was not run

entrez_genes<- lapply(heatmap_genes, function(x){
  mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL") } )

#### Adeno ####
adeno_go<- enrichGO(gene= entrez_genes$Adeno,
      OrgDb        = org.Mm.eg.db,
      ont          = "BP",
      maxGSSize = 200,
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      readable = TRUE)


write.table(adeno_go@result, file= "./Results/GSEA/VariableGenes_Adenoma_BP.txt")

goplot(adeno_go, showCategory = 10)
plot<- dotplot(adeno_go, showCategory=10) + ggtitle("Adenoma Cluster \nVariable Genes BP Enrichment")

plot_dir<- "./Plots/Requested_Plots/"
pdf(file=paste0(plot_dir, "AdenoHeatmapClusterGenes_BPEnrich_DotPlot.pdf"), 
    width = 8, height=6)
plot
dev.off()

adeno_go <- pairwise_termsim(adeno_go)
emapplot(adeno_go)

# reactome <- enrichPathway(entrez_genes$Adeno,
#                        organism= "mouse",
#                        maxGSSize = 200,
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff  = 0.05,
#                        pAdjustMethod = "fdr", 
#                        readable = TRUE)
# 
# dotplot(reactome, showCategory=15) + ggtitle("dotplot for Adeno")


#### Clear cell ####
cc_go<- enrichGO(gene= entrez_genes$Clear_cell,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",
                    maxGSSize = 200,
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05,
                    readable = TRUE)

write.table(cc_go@result, file= "./Results/GSEA/VariableGenes_ClearCell_BP.txt")


goplot(cc_go, showCategory = 10)
plot2<- dotplot(cc_go, showCategory=10) + ggtitle("Clear Cell Cluster \nVariable Genes BP Enrichment")

pdf(file=paste0(plot_dir, "ClearCellHeatmapClusterGenes_BPEnrich_DotPlot.pdf"), 
    width = 8, height=6)
plot2
dev.off()

cc_go <- pairwise_termsim(cc_go)
emapplot(cc_go)

# reactome <- enrichPathway(entrez_genes$Clear_cell,
#                           organism= "mouse",
#                           maxGSSize = 200,
#                           pvalueCutoff = 0.05,
#                           qvalueCutoff  = 0.05,
#                           pAdjustMethod = "fdr", 
#                           readable = TRUE)
# 
# dotplot(reactome, showCategory=15) + ggtitle("dotplot for Adeno")
