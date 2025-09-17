library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(stringr)

setwd("/your/parent/dir/")

files<- list.files("./Robjs", pattern = "_seurat_unfiltered.rds", full.names = T)

objs<- lapply(files, readRDS)

set.seed(123)

#Get scaling value for point size per sample
pts<- lapply(objs, function(x){
  0.017*(1/x@images[[paste0("slice_", unique(x$orig.ident))]]@spot.radius)
})

samples<- lapply(objs, function(x){unique(x$orig.ident)} ) %>% 
  as.character()
names(pts)<- samples

# pulled from the Seurat source code--requires package RColorBrewer
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, 
                                                              name = "Spectral")))

# p<- SpatialFeaturePlot(sp.merged, features = c("Amh"),crop = F) &
#   theme( legend.position = "right", plot.title = element_blank())

#### Meta column plotting function ####
meta_plot<- function(obj, meta_cols, lims){
  samp<- unique(obj$orig.ident)
  ptsz<- pts[samp] %>% unlist
  
  p<- SpatialFeaturePlot(obj, meta_cols,  ncol = length(meta_cols), pt.size.factor = ptsz, stroke = NA, crop = T) &
    scale_fill_gradientn(limits=lims, colours=SpatialColors(n=100)) 
  
  #If meta values to be plotted are the UCD scores, change titles
  if(all(meta_cols %in% names(ucd_cols_titles))){
    for(i in 1:length(p)){
      title<- ucd_cols_titles[p[[i]]$labels$fill]
      p[[i]]<- p[[i]] + ggtitle(title) +theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
  } else{ 
    for(i in 1:length(p)){
      title<- p[[i]]$labels$fill
      p[[i]]<- p[[i]] + ggtitle(title) +theme(plot.title = element_text(hjust = 0.5, size = 10))
      }
  }
  
  p<- p& theme(legend.position = "right", legend.title = element_blank())
  
  p_leg<- get_legend(p[[1]])
  
  #If meta values to be plotted are the UCD scores, remove UCD for plot names
  if(all(meta_cols %in% names(ucd_cols_titles))){
    pdf(file = paste0(plot_dir, samp, "_", str_flatten(str_remove_all(meta_cols, "UCD_"), "_"), ".pdf"), 
        width = 8, height = 4)
    plot_grid(p & NoLegend(), p_leg,
              ncol = 2, rel_widths = c(1, 0.1)) %>% print()
    dev.off()
  } else{ #otherwise, just plot values normally
    pdf(file = paste0(plot_dir, samp, "_", str_flatten(str_remove_all(meta_cols, "UCD_"), "_"), ".pdf"), 
        width = 8, height = 4)
    plot_grid(p & NoLegend(), p_leg,
              ncol = 2, rel_widths = c(1, 0.1)) %>% print()
    dev.off()
  }
  
  plot<- plot_grid(p & NoLegend(), p_leg,
                   ncol = 2, rel_widths = c(1, 0.1))
  return(plot)
}


#### Deconvolution scores for cell types of interest ####
ucd_cols<- c("leydig.cell", "seminiferous.tubule.epithelial.cell","sertoli.cell",
             "granulosa.cell","follicular.cell.of.ovary",
             "meso.epithelial.cell", "cumulus.cell", "ovarian.surface.epithelial.cell")

#Add line breaks to make plotting cell names more compact
ucd_cols_titles<- str_replace_all(ucd_cols, "[.]", " ")
ucd_cols_titles<-str_replace_all(ucd_cols_titles, "seminiferous ", "seminiferous\n")
ucd_cols_titles<- str_replace_all(ucd_cols_titles, "follicular cell ", "follicular cell\n")
ucd_cols_titles<- str_replace_all(ucd_cols_titles, "meso epithelial ", "meso-epithelial\n")
ucd_cols_titles<- str_replace_all(ucd_cols_titles, "ovarian surface ", "ovarian surface\n")
names(ucd_cols_titles)<- paste0("UCD_", ucd_cols)

for (i in 1:length(objs)){
  samp<- objs[[i]]$orig.ident %>% unique()
  ucd_filt<- read.csv(paste0("./UCD_output/", samp, "_primary_propagated.csv"), 
                      sep = "\t")
  ucd_filt<- ucd_filt %>% tibble::column_to_rownames("X")
  #Omit cell and somatic cell 
  ucd_filt<- ucd_filt[,ucd_cols]
  
  cells<- intersect(rownames(ucd_filt), colnames(objs[[i]]))
  
  objs[[i]]@meta.data[cells, paste0("UCD_", ucd_cols)]<- ucd_filt[cells,]
  
}


plot_dir<- "./Plots/Requested_Plots/UCD_scores/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

##### Apply meta function to UCD columns ####
ucd_plots<- lapply(objs, meta_plot, meta_cols= paste0("UCD_", ucd_cols[1:4]), lims= c(0,1))

ucd_plots<- lapply(objs, meta_plot, meta_cols= paste0("UCD_", ucd_cols[5:8]), lims= c(0,1))

# plot_grid(plotlist =ucd_plots, nrow = 4)

#### Gene feature plots ####

plot_dir<- "./Plots/Requested_Plots/Genes/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)


genes<- c("Wt1", "Sox8", "Sox9", "Dmrt1", 
          "Amh", "Krt19", "Pcp4", "Sprr2f", 
          "Bcam", "Cdh2", "Foxl2", "Rras2", 
          "Pax8", "Wnt2b", "Upk1b", "Unc45b")

#get the max expression value from 
#all genes in all samples to keep scale constant in every plot
scale.max<- lapply(objs, function(obj){
  max(obj@assays$SCT@data[genes,])
}) %>% as.numeric() %>% max()

gene_plot<- function(obj, genes){
  DefaultAssay(obj)<- "SCT"
  samp<- unique(obj$orig.ident)
  ptsz<- pts[samp] %>% unlist
  
  p<- SpatialFeaturePlot(obj, genes,  ncol = length(genes), pt.size.factor = ptsz, stroke = NA, crop = T) &
    scale_fill_gradientn(limits=c(0, scale.max), colours=SpatialColors(n=100)) 
  
  for(i in 1:length(p)){
    title<- p[[i]]$labels$fill
    p[[i]]<- p[[i]] + ggtitle(title) +theme(plot.title = element_text(hjust = 0.5, size = 10))
  }
  
  p<- p& theme(legend.position = "right", legend.title = element_blank())
  
  p_leg<- get_legend(p[[1]])
  
  pdf(file = paste0(plot_dir, samp, "_", str_flatten(genes, "_"), ".pdf"), 
      width = 8, height = 4)
  plot_grid(p & NoLegend(), p_leg,
            ncol = 2, rel_widths = c(1, 0.1)) %>% print()
  dev.off()
  
  plot<- plot_grid(p & NoLegend(), p_leg,
                   ncol = 2, rel_widths = c(1, 0.1))
  return(plot)
}

#Plot the genes in group of 4
lapply(objs, gene_plot, genes= genes[1:4])
lapply(objs, gene_plot, genes= genes[5:8])
lapply(objs, gene_plot, genes= genes[9:12])
lapply(objs, gene_plot, genes= genes[13:16])

#### Pathway scores (from progeny) ####

plot_dir<- "./Plots/Requested_Plots/Progeny_scores/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

for (i in 1:length(objs)){
  samp<- objs[[i]]$orig.ident %>% unique()
  prog_df<- read.table(paste0("./Results/Progeny_Results/", samp, "_full_progeny_df.txt"), 
                         check.names= F) %>% t()
  
  objs[[i]]@meta.data[rownames(prog_df), colnames(prog_df)]<- prog_df
  
}

pathways<- c("MAPK", "PI3K")

scale.max<- lapply(objs, function(obj){
  max(obj@meta.data[,pathways])
}) %>% as.numeric() %>% max()

scale.min<- lapply(objs, function(obj){
  min(obj@meta.data[,pathways])
}) %>% as.numeric() %>% min()

proj_plots<- lapply(objs, meta_plot, meta_cols= pathways, lims= c(scale.min,scale.max))
