library(Seurat)
library(cowplot)
library(ggplot2)
library(hdf5r)
library(stringr)

setwd("/gpfs/scratch/bsc64/bsc64717/Bustelo_OV")

#### Create Seurat Objects ####

dirs<- c("Ovary_522_v2", "Ovary_547_v2",
         "Ovary_708", "Ovary_818_v2")
samps<- paste0("OVA", str_split_fixed(dirs, "_", n= 3)[,2])

##### Sample processing #####
for (i in 1:length(samps)){
  dir<- dirs[i]
  samp<- samps[i]
  seurat<- Load10X_Spatial(data.dir = paste0("./", dir, "/outs/"))
  seurat$orig.ident<- samp
  Idents(seurat)<- seurat$orig.ident
  
  seurat<- SCTransform(seurat, vst.flavor= "v2", method = "glmGamPoi", 
                          return.only.var.genes = F, assay = "Spatial")
  
  seurat<- RunPCA(seurat, assay = "SCT", verbose = FALSE)
  seurat<- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
  seurat<- FindClusters(seurat, verbose = FALSE)

  seurat<- RunUMAP(seurat, reduction = "pca", dims = 1:30)
  
  saveRDS(seurat, file= paste0("./Robjs/", samp, "_seurat_unfiltered.rds"))
}

##### QC/Cluster plots #####
files<- list.files("./Robjs", full.names = T)

objs<- lapply(files, readRDS)

#change point size depending on sample
pts<- c(3, 2, 2.5, 1.2)
names(pts)<- samps

for (i in 1:length(samps)){
  dir<- dirs[i]
  samp<- samps[i]
  
  seurat<- objs[[i]]
  
  Idents(seurat)<- seurat$orig.ident
  
  pt<- as.numeric(pts[samp])
  
  #Counts plots for QC
  plot1 <- VlnPlot(seurat, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() +
    ggtitle("Counts per spot") + 
    theme(plot.title = element_text(hjust = 0.5))
  plot2 <- SpatialFeaturePlot(seurat, features = "nCount_Spatial", 
                              alpha = c(0.5, 1),pt.size.factor = pt ) + 
    theme(legend.position = "bottom",  legend.text =  element_text(angle = 45,hjust=1)) + labs(fill= "Counts per spot")

  ifelse(!dir.exists("./Plots/Spatial/QC/"), dir.create("./Plots/Spatial/QC/", recursive = T),
         FALSE)
  
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_nCountSpatial.pdf"), width = 10, height = 7)
  print(plot_grid(plot1, plot2))
  dev.off()
  
  #Feature plots for QC 
  plot1 <- VlnPlot(seurat, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend() +
    ggtitle("Genes per spot") + 
    theme(plot.title = element_text(hjust = 0.5))
  plot2 <- SpatialFeaturePlot(seurat, features = "nFeature_Spatial", 
                              alpha = c(0.5, 1),pt.size.factor = pt ) + 
    theme(legend.position = "bottom",  legend.text =  element_text(angle = 45,hjust=1)) + 
    labs(fill= "Genes per spot")

  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_nFeatureSpatial.pdf"), width = 10, height = 7)
  print(plot_grid(plot1, plot2))
  dev.off()
  
  #Low feature spots highlight 
  plot1<- VlnPlot(seurat, features = "nFeature_Spatial") + NoLegend() + 
    geom_abline(intercept=1000, 
                slope = 0 , color = "red", linetype=3) + ggtitle("Genes per spot") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  spots_low<- colnames(seurat)[seurat$nFeature_Spatial < 1000]
  
  plot2<- SpatialDimPlot(seurat, cells.highlight = spots_low) +theme(legend.position = "bottom") + 
    labs(fill= "Spots gene count") + 
    scale_fill_manual(labels = c("<1000", ">=1000"), values = c("grey", "red"))
  
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_LowFeatureSpatial.pdf"), width = 10, height = 7)
  print(plot_grid(plot1, plot2))
  dev.off()
  
  #Top variable features- scatter & tissue location 
  var_top10 <- head(VariableFeatures(seurat), 10)
  var_p1 <- LabelPoints(plot = VariableFeaturePlot(seurat), points = seurat@assays[["SCT"]]@var.features[1:10], repel = TRUE)+
    ggtitle(paste0(samp, " Variable Genes")) + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_VarGeneScatter.pdf"), width = 6, height = 6)
  print(var_p1)
  dev.off()
  
  var_p2<- SpatialFeaturePlot(seurat, features = var_top10[1:10], ncol=5,
                              pt.size.factor = pt, alpha = 0.5, keep.scale = "all")
  
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_VarGeneSpatial.pdf"), width = 12, height = 6)
  print(var_p2)
  dev.off()
  
  #Cluster plot 
  Idents(seurat)<- seurat$seurat_clusters
  
  plot1<- SpatialDimPlot(seurat, label = T, pt.size.factor = pt)
  
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_SpatialClusters.pdf"), width = 6, height = 6)
  print(plot1)
  dev.off()
  
}

#FindMarkers, make heatmap 

files<- list.files("./Robjs", full.names = T)

objs<- lapply(files, readRDS)

samps<- as.character(lapply(objs, function(x){ return(unique(x$orig.ident))}))

markers<- lapply(objs, FindAllMarkers)

dir.create("./Results/Res0.8ClusterMarks")

save_marks<- function(df, samp){
  df[,"sample"]<- samp
  write.csv(df, paste0("./Results/Res0.8ClusterMarks/", samp, "_Res0.8_ClusterMarkers.csv"))
}

mapply(save_marks, markers, samps)

library(RColorBrewer)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"))(256)


marker_hm<- function(obj, mark){
  top10<- mark %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 3) %>%
    ungroup()
  plot<- DoHeatmap(obj, features = top10$gene, angle = 90, size =3, 
                   slot= "data") + 
   scale_fill_gradientn(colours = mapal) +
  guides(colour="none")
}

hms<- mapply(marker_hm, objs, markers)

mapply(hms, samps, FUN= function(hm, samp){
  pdf(file= paste0("./Plots/Spatial/QC/", samp, "_Res0.8Heatmap.pdf"), width = 10, height = 7)
  print(hm)
  dev.off()
})