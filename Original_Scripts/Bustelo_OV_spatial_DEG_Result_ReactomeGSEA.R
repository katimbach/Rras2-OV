library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(DOSE)

setwd("/home/bscuser/Documents/Spatial/Bustelo_OV")

sp.merged<- readRDS("./Robjs/AllSamp_seurat_regionsubset_merged.rds")

#plot genes for atretic follicles rel. to healthy follicles (consistenmt postive)
SpatialFeaturePlot(seurat, features = c("H2-D1", "Selenop", "Dcn"), images = "slice_OVA547")

#plot atretic-defining markers (relative to all other areas)
SpatialFeaturePlot(seurat, features = c("Ptprq", "Serpina5", "Angpt2"), images = "slice_OVA547")

#compare genes for evolution
plot_grid(SpatialFeaturePlot(seurat, features = c("4930486L24Rik", "Defb36", "Gatm"), images = "slice_OVA708", 
                   pt.size.factor = 2.5), 
          SpatialFeaturePlot(seurat, features = c("4930486L24Rik", "Defb36", "Gatm"), images = "slice_OVA818", 
                             pt.size.factor = 1.2), nrow = 2)


#### GSEA/GO analysis for gene lists ####

##### atretic spots from 547 ####

#look at atret v healthy follicles from 547
atret_genes<- read.table("./Results/DEGTestResults/OVA547_atretic_v_healthy_follicle_markers.txt")

atret_lfc<- atret_genes$avg_log2FC
names(atret_lfc)<- rownames(atret_genes)

gene_list = sort(atret_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

atret_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

# gse <- gseGO(geneList=gene_list, 
#              ont ="BP", 
#              keyType = "ENTREZID",
#              minGSSize = 50, 
#              maxGSSize = 200, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Mm.eg.db, 
#              pAdjustMethod = "fdr")
# 
dotplot(atret_res, showCategory=10, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)


#look at atret v all else
atret_genes<- read.table("./Results/SampleSpotIdentMarkers/OVA547_follicle_atretic_markers.txt")

atret_lfc<- atret_genes$avg_log2FC
names(atret_lfc)<- rownames(atret_genes)

gene_list = sort(atret_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

atret_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

# gse <- gseGO(geneList=gene_list, 
#              ont ="BP", 
#              keyType = "ENTREZID",
#              minGSSize = 50, 
#              maxGSSize = 200, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Mm.eg.db, 
#              pAdjustMethod = "fdr")
# 
# dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)

#GPCR is consistent, get information 
gpcr_genes<- reactome@result[["core_enrichment"]][which(reactome@result[["Description"]] ==  "GPCR ligand binding")]
library(stringr)
gpcr_genes<- str_split(gpcr_genes, "[/]", simplify = T)

names(gene_map[gene_map %in% gpcr_genes])
# [1] "Galr1"    "Cga"      "Gpr37"    "Tas2r118" "Uts2"     "Prok2"    "Chrm5"    "Lhcgr"    "Chrm1"    "Insl5"   
# [11] "Gabbr2"   "Taar2"    "Htr5a"    "Fshr"     "Oprl1"    "Sstr1"    "Gpr37l1"  "Ackr2"    "Taar1"    "Ccl20"   
# [21] "Grm7"     "Drd4"     "Grm6"     "Ucn"      "Apln"     "Ffar3"    "Hcrt"     "Cxcl1"    "Sct"      "Edn3"    
# [31] "Ccl22"    "Lpar4"    "Edn1"     "Gpr143"   "Gpha2"    "Ptger2"   "Lpar2"    "Penk"     "Pomc"     "Ghrhr"   
# [41] "Sstr3"    "Npff"     "Htr1d"    "F2rl3"    "C5ar2"    "Ednra"

gseaplot(reactome, by = "all", title = "GPCR ligand binding", 
         geneSetID = which(reactome@result[["Description"]] ==  "GPCR ligand binding"))


##### sertoli/clear cell vs adenoma #####
cc_v_ad_genes<- read.table("./Results/DEGTestResults/AllMutSamples_sertoli_cc_vs_adeno_markers.txt")

cc_v_ad_lfc<- cc_v_ad_genes$avg_log2FC
names(cc_v_ad_lfc)<- rownames(cc_v_ad_genes)

gene_list = sort(cc_v_ad_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)

#Look at just the DEGs from 708 only
cc_v_ad_genes<- read.table("./Results/DEGTestResults/OVA708_sertoli_cc_vs_adeno_markers.txt")

cc_v_ad_lfc<- cc_v_ad_genes$avg_log2FC
names(cc_v_ad_lfc)<- rownames(cc_v_ad_genes)

gene_list = sort(cc_v_ad_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

cc_v_ad_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)


#And DEGs from 818 only
cc_v_ad_genes<- read.table("./Results/DEGTestResults/OVA818_sertoli_cc_vs_adeno_markers.txt")

cc_v_ad_lfc<- cc_v_ad_genes$avg_log2FC
names(cc_v_ad_lfc)<- rownames(cc_v_ad_genes)

gene_list = sort(cc_v_ad_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

cc_v_ad_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)

##### granulosa tumor vs sertoli/clear cell and adenoma #####
gran_genes<- read.table("./Results/DEGTestResults/OVA818_granulosa_tumor_vs_sertoli_cc_and_adeno_markers.txt")

gran_lfc<- gran_genes$avg_log2FC
names(gran_lfc)<- rownames(gran_genes)

gene_list = sort(gran_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

gran_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)

##### innervated stroma #####
#Out of curiousity, look at innervated stroma group from 818 to see if it indeed has neuronal markers
innerv_genes<- read.table("./Results/SampleSpotIdentMarkers/OVA818_innerv_stroma_markers.txt")

innerv_lfc<- innerv_genes$avg_log2FC
names(innerv_lfc)<- rownames(innerv_genes)

gene_list = sort(innerv_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

innerv_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)

##### sertoli/clear cell phenotypes (more vs less clear) #####

#Look at just the DEGs from 708 only
cc_v_ad_genes<- read.table("./Results/DEGTestResults/OVA708_sertoli_cc_cluster9_vs_cluster2_markers.txt")

cc_v_ad_lfc<- cc_v_ad_genes$avg_log2FC
names(cc_v_ad_lfc)<- rownames(cc_v_ad_genes)

gene_list = sort(cc_v_ad_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]

cc_v_ad_res<- gseDO(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
) # no result

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)


#And DEGs from 818 only
cc_v_ad_genes<- read.table("./Results/DEGTestResults/OVA818_sertoli_cc_cluster5_vs_cluster0_markers.txt")

cc_v_ad_lfc<- cc_v_ad_genes$avg_log2FC
names(cc_v_ad_lfc)<- rownames(cc_v_ad_genes)

gene_list = sort(cc_v_ad_lfc, decreasing = TRUE)

gene_map<- mapIds(org.Mm.eg.db, names(gene_list), "ENTREZID", "SYMBOL")

names(gene_list)<- gene_map

sum(duplicated(names(gene_list)))
which(duplicated(names(gene_list)))

gene_list<- gene_list[! duplicated(names(gene_list))]


gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 50,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

reactome <- gsePathway(gene_list,
                       organism= "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr", 
                       verbose = FALSE)

dotplot(reactome, showCategory=5, split=".sign", font.size=10) + facet_grid(.~.sign)
gseaplot(reactome, by = "all", title = reactome$Description[1], geneSetID = 1)


#read in the clear cell results run in HPC 

ifelse(!dir.exists("./Plots/Spatial/GSEA"), dir.create("./Plots/Spatial/GSEA"),FALSE)

cc_708<- readRDS("./Results/GSEA/OVA708_adv_cc_vs_cc_deg_gsea.rds")
dotp<- dotplot(cc_708, showCategory=5, split=".sign", font.size=9,
               label_format=50) + facet_grid(.~.sign)

pdf(file="./Plots/Spatial/GSEA/OVA708_adv_clear_cell_gsea_dot.pdf", width=8, height=4)
dotp
dev.off()


num<- which(cc_708$Description == "ERK1 and ERK2 cascade")
gseap<- gseaplot(cc_708, geneSetID = num, 
                 by = "runningScore", title = cc_708$Description[num])

pdf(file="./Plots/Spatial/GSEA/OVA708_adv_clear_cell_gsea_ERK1ERK2.pdf", width=8, height=4)
gseap
dev.off()

cc_818<- readRDS("./Results/GSEA/OVA818_adv_cc_vs_cc_deg_gsea.rds")
dotp<- dotplot(cc_818, showCategory=5, split=".sign", font.size=9,
               label_format=50) + facet_grid(.~.sign)

pdf(file="./Plots/Spatial/GSEA/OVA818_adv_clear_cell_gsea_dot.pdf", width=6, height=4)
dotp
dev.off()

num<- which(cc_818$Description == "ERK1 and ERK2 cascade")
gseap<- gseaplot(cc_818, geneSetID = num, 
                 by = "runningScore", title = cc_818$Description[num])

pdf(file="./Plots/Spatial/GSEA/OVA818_adv_clear_cell_gsea_ERK1ERK2.pdf", width=8, height=4)
gseap
dev.off()
