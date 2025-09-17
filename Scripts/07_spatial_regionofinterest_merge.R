library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(cowplot)

setwd("/your/parent/dir/")

anno_path<- "./Input_Files/AllSamp_spot_annotations.txt"

anno_df<- read.table(anno_path)

anno_cols<- c("region_of_interest", "seurat_clusters", "spot_ident")

files<- list.files("./Robjs", pattern = "_seurat_unfiltered.rds", full.names = T)

objs<- lapply(files, readRDS)

objs<- lapply(objs, function(x){
  DefaultAssay(x)<- "SCT"
  return(x)
})

objs_sub<- lapply(objs, function(x){
  anno_sub<- anno_df %>% filter(orig.ident == unique(x$orig.ident))
  x@meta.data[anno_sub$barcode, anno_cols]<- anno_sub[,anno_cols]
  return(x[,x$region_of_interest == "TRUE"])
})

##### Merged subset object #####
objs_sub<- lapply(objs_sub, function(x){
  x$barcode<- colnames(x)
  names(x@images)<- paste0("slice_", unique(x$orig.ident))
  return(x)
})

objs_sub<- lapply(objs_sub, SCTransform, vst.flavor= "v2", method = "glmGamPoi", 
                          return.only.var.genes = F, assay = "Spatial")


sp.merged<- merge(objs_sub[[1]], objs_sub[c(2:length(objs_sub))])
sp.merged<- RenameCells(sp.merged,  
                        new.names = paste0(sp.merged$orig.ident, "_", sp.merged$barcode))


features <- SelectIntegrationFeatures(object.list = objs_sub, nfeatures = 3000)

VariableFeatures(sp.merged)<- features

sp.merged<- RunPCA(sp.merged)
#ElbowPlot(sp.merged, ndims = 50)

#DimPlot(sp.merged, group.by = "orig.ident") +
#  FeaturePlot(sp.merged, "nCount_Spatial", order = T)
#DimPlot(sp.merged, group.by = "orig.ident") +
#  FeaturePlot(sp.merged, "estimate", order = T)

sp.merged <- RunUMAP(sp.merged, reduction = "pca", dims = 1:50)

#DimPlot(sp.merged, group.by = "orig.ident") +
#  FeaturePlot(sp.merged, "estimate", order = T)
#Strong batch effect present- run CCA

sp.merged<- PrepSCTFindMarkers(sp.merged)

saveRDS(sp.merged, file= "./Robjs/AllSamp_seurat_regionsubset_merged.rds")