library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(msigdbr)
library(purrr)
library(tidyr)

setwd("/home/bscuser/Documents/Spatial/Bustelo_OV")

sp.merged<- readRDS("./Robjs/AllSamp_seurat_regionsubset_merged.rds")

subset<- sp.merged[,sp.merged$spot_ident %in% c("sertoli_like_cc", "sertoli_adeno_border", "adeno_like", "granulosa_tumor")]

h_gene_sets = msigdbr(species = "mouse", category = "H") %>% as.data.frame()
h_list<- with(h_gene_sets, split(gene_symbol, gs_name))

lapply(h_list, function(x){
  sum(x %in% rownames(subset)) / length(x)
})

subset<- AddModuleScore(subset, features = h_list, assay = "SCT", search = T, name = names(h_list))
colnames(subset@meta.data)[grepl("HALLMARK", colnames(subset@meta.data))]<- names(h_list)

df<- subset@meta.data[,c("spot_ident", names(h_list))]

dir.create("./Results/MSig_Hallmarks")

write.table(df, "./Results/MSig_Hallmarks/Hallmarks_Scores.txt")

VlnPlot(subset, features = names(h_list), group.by = "orig.ident")

h_median<- subset@meta.data %>% 
  summarise(across(all_of(names(h_list)), median))

h_median_cluster<- subset@meta.data %>% group_by(samp_cluster) %>%
  summarise(across(all_of(names(h_list)), median))

h_variance<- subset@meta.data %>% #group_by(samp_cluster) %>%
  summarise(across(all_of(names(h_list)), var))

which.max(h_variance)
# HALLMARK_CHOLESTEROL_HOMEOSTASIS 
# 9 

h_variance[order(h_variance, decreasing = T)]

#Test which pathways are different between adeno and clear cell areas? 
library(rstatix)
df_sub<- df[df$spot_ident %in% c("sertoli_like_cc", "adeno_like"),]
for(h in names(h_list)){
  if(which(names(h_list) == h) == 1){
    res<- df_sub %>% wilcox_test(as.formula(paste(h, " ~ spot_ident", sep = "")))
    effect<- df_sub %>% wilcox_effsize(as.formula(paste(h, " ~ spot_ident", sep = "")))
    res[,c("effsize", "magnitude")]<- effect[,c("effsize", "magnitude")]
  }
   
  if(which(names(h_list) == h) > 1){
    res2<- df_sub %>% wilcox_test(as.formula(paste(h, " ~ spot_ident", sep = "")))
    effect<- df_sub %>% wilcox_effsize(as.formula(paste(h, " ~ spot_ident", sep = "")))
    res2[,c("effsize", "magnitude")]<- effect[,c("effsize", "magnitude")]
    res<- rbind(res, res2)
  }
}

res$p.adj<- p.adjust(res$p, method = "fdr")
write.table(res, "./Results/MSig_Hallmarks/Hallmarks_ClearCell_Adeno_WilcoxRes.txt")

#plot the top 3 on the tissues (fourth doesn't have much expression)
dir.create("./Plots/Spatial/MSig_Hallmarks")

hallmarks<- res %>% top_n(-3, p.adj) %>% pull(.y.)

pts<- c(3, 2, 2.5, 1.2)

for (h in hallmarks){
  p<- SpatialFeaturePlot(subset, features =  h, 
                         keep.scale = "all", images = names(subset@images)[2:4]) 

  leg<- ggpubr::get_legend(p)
  p[[1]][["plot_env"]][["plot"]][["layers"]][[1]][["aes_params"]][["point.size.factor"]] <- pts[3]
  p[[2]][["plot_env"]][["plot"]][["layers"]][[1]][["aes_params"]][["point.size.factor"]] <- pts[3]
  p[[3]][["plot_env"]][["plot"]][["layers"]][[1]][["aes_params"]][["point.size.factor"]] <- pts[4]
  
  plot<- plot_grid(leg, p & NoLegend(), nrow=2,rel_heights = c(0.1,1))
  
  pdf(file= paste0("./Plots/Spatial/MSig_Hallmarks/MaligAreas_", h, "_score.pdf"), width = 12, height = 6)
  plot(plot)
  dev.off()
  
  v<- VlnPlot(subset[,subset$spot_ident %in% c("adeno_like", "sertoli_like_cc")], 
          features = h, group.by = "spot_ident") +  ggtitle(str_replace_all(h, "_", " ")) +
    ggpubr::stat_compare_means(label =  "p.signif", label.x = 1.5) & NoLegend() 
  pdf(file= paste0("./Plots/Spatial/MSig_Hallmarks/Adeno_ClearCell_", h, "_violin.pdf"), width = 6, height = 6)
  plot(v)
  dev.off() 
}

