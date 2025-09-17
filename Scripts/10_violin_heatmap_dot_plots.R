library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(ggpubr)

set.seed(123)

setwd("/your/parent/dir/")

sp.merged<- readRDS("./Robjs/AllSamp_seurat_regionsubset_merged.rds")

sp.merged@meta.data <- sp.merged@meta.data %>% mutate(sample_name = case_when(orig.ident == "OVA522" ~ "WT",
                                                                              orig.ident == "OVA547" ~ "HT",
                                                                              orig.ident == "OVA708" ~ "HM AdCre Rep1",
                                                                              orig.ident == "OVA818" ~ "HM AdCre Rep2"))

sp.merged$spot_ident_spec<- sp.merged$spot_ident

sp.merged$spot_ident_spec[sp.merged$spot_ident %in% c("rete", "rete1")]<- "normal_rete"
sp.merged$spot_ident_spec[sp.merged$spot_ident %in% c("rete2")]<- "cystic_rete"

sp.merged$spot_ident_spec[sp.merged$spot_ident %in% c("sertoli_like_cc")]<- "clear_cell"

sp.merged$spot_ident_spec[sp.merged$spot_ident %in% c("sertoli_like_cc") & sp.merged$SCT_snn_res.1 %in% c(9, 5)]<- "clear_cell_adv"

pts<- sapply(names(sp.merged@images), function(x){
  0.017*(1/sp.merged@images[[x]]@spot.radius)
})

samples<- unique(sp.merged$orig.ident)
pts_vector<- list(rep(pts[1] , sum(sp.merged$orig.ident == samples[1])), 
                  rep(pts[2] , sum(sp.merged$orig.ident == samples[2])),
                  rep(pts[3] , sum(sp.merged$orig.ident == samples[3])),
                 rep(pts[4] , sum(sp.merged$orig.ident == samples[4]))) 
names(pts_vector)<- names(sp.merged@images)

p<- SpatialFeaturePlot(sp.merged, features = c("Amh"),crop = F) &
  theme( legend.position = "right", plot.title = element_blank())


ucd_cols<- c("leydig.cell", "seminiferous.tubule.epithelial.cell","sertoli.cell",
             "granulosa.cell","follicular.cell.of.ovary",
             "meso.epithelial.cell", "cumulus.cell", "ovarian.surface.epithelial.cell")

sp.merged@meta.data[, paste0("UCD_", ucd_cols)]<- 0

for (i in 1:length(samples)){
  samp<- samples[i]
  ucd_filt<- read.csv(paste0("./Results/UCD_output/", samp, "_primary_propagated.csv"), 
                      sep = "\t")
  ucd_filt<- ucd_filt %>% tibble::column_to_rownames("X")
  rownames(ucd_filt)<- paste0(samp, "_", rownames(ucd_filt))
  
  #Keep only relevant cell types
  ucd_filt<- ucd_filt[rownames(ucd_filt) %in% colnames(sp.merged) ,ucd_cols]
  
  sp.merged@meta.data[rownames(ucd_filt), paste0("UCD_", ucd_cols)]<- ucd_filt
  
}

# pulled from the Seurat source code--requires package RColorBrewer
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))


p<- SpatialFeaturePlot(sp.merged, paste0("UCD_", ucd_cols[1:4]), images = "slice_OVA547", ncol = 4) &
  scale_fill_gradientn(limits=c(0, 1), colours=SpatialColors(n=100)) & 
  theme(legend.title.position = "top", legend.position = "bottom") &ggtitle(c(ucd_cols[1:4]))

for(i in 1:length(p)){
  p[[i]]<- p[[i]] + ggtitle(p[[i]]$labels$fill)
}

p<- p& theme(legend.position = "right", legend.title = element_blank())

p_leg<- get_legend(p[[1]])

plot<- plot_grid(p & NoLegend(), p_leg,
                 ncol = 2, rel_widths = c(1, 0.1))

#### Violin plot of UCD leydig, semin tubule and sertoli per sample & tumor type ####
tumor_cols<- c("adeno_like"= "#8EC7FF",
               "clear_cell" = "#E6797C", 
               "clear_cell_adv"= "#018C87")

sp.merged$groups_compare<- paste0(sp.merged$sample_name, " ", sp.merged$spot_ident_spec)
sp.merged$groups_compare[is.na(sp.merged$groups_compare)]<- " "

df<- sp.merged@meta.data %>% filter(spot_ident %in% c("sertoli_like_cc", "adeno_like")) 

labels<- df$sample_name
names(labels)<- df$groups_compare

#Factor the groups to have the order HT -> Adeno -> CC groups
sp.merged$groups_compare<- factor(sp.merged$groups_compare, 
                                  levels = c("HT clear_cell", "HM AdCre Rep1 adeno_like", "HM AdCre Rep1 clear_cell",    
                                             "HM AdCre Rep1 clear_cell_adv","HM AdCre Rep2 adeno_like", "HM AdCre Rep2 clear_cell",
                                             "HM AdCre Rep2 clear_cell_adv", " "))
#Specify the colors for the groups 
groups<- levels(sp.merged$groups_compare)

tumor_cols_exp<- rep("#8EC7FF", length(groups))

tumor_cols_exp[groups %>% endsWith(., "clear_cell")]<- "#E6797C"
tumor_cols_exp[groups %>% endsWith(., "clear_cell_adv")]<- "#018C87"
names(tumor_cols_exp)<- groups

#Add labels only for one of each tumor group
labels_legend<-  c( "HM AdCre Rep1 adeno_like" ="adenoma", #"HM AdCre Rep2 adeno_like" ="adenoma",
                    "HT clear_cell" = "early\nclear cell", #"HM AdCre Rep1 clear_cell" = "early\nclear cell", "HM AdCre Rep2 clear_cell" ="early\nclear cell",    
                    "HM AdCre Rep1 clear_cell_adv" ="advanced\nclear cell") #, "HM AdCre Rep2 clear_cell_adv"= "advanced\nclear cell")

#Make comparisons list for stat_compare_means
  #Compare all early clear cell, within tumor types per sample
comps<- c(lapply(c(3:4), function(x){ c(groups[2], groups[x])}),
          lapply(c(6:7), function(x){ c(groups[5], groups[x])}),
          lapply(c(3,6), function(x){ c(groups[1], groups[x])}))

ymnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                   symbols = c("****", "***", "**", "*", "ns"))

#No stat test for leydig cells, all very low
p<- lapply( paste0("UCD_", ucd_cols[1:3]), function(col){
  if(col == "UCD_leydig.cell"){
    VlnPlot(sp.merged, col, group.by =  "groups_compare",
            idents = c("sertoli_like_cc", "adeno_like"), cols = tumor_cols_exp)  &
      scale_y_continuous(breaks = seq(0, 1, 0.2)) &  # Only show ticks up to 1.0
      coord_cartesian(ylim = c(0, 1.6)) &  # Extend plotting space but keep tick labels fixed &
      theme(legend.position = "bottom") & 
      scale_x_discrete(labels= labels) &
      scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                        breaks = c("HM AdCre Rep1 adeno_like", "HT clear_cell",
                                   "HM AdCre Rep1 clear_cell_adv")) &
      guides(y = guide_axis(cap = "upper")) & 
      xlab("Sample") & ylab("Score")
  }else{
    VlnPlot(sp.merged, col, group.by =  "groups_compare",
            idents = c("sertoli_like_cc", "adeno_like"), cols = tumor_cols_exp)  &
      scale_y_continuous(breaks = seq(0, 1, 0.2)) &  # Only show ticks up to 1.0
      coord_cartesian(ylim = c(0, 1.6)) &  # Extend plotting space but keep tick labels fixed &
      theme(legend.position = "bottom") & 
      scale_x_discrete(labels= labels) &
      scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                        breaks = c("HM AdCre Rep1 adeno_like", "HT clear_cell",
                                   "HM AdCre Rep1 clear_cell_adv")) &
      guides(y = guide_axis(cap = "upper")) & 
      xlab("Sample") & ylab("Score") &
      stat_compare_means(comparisons = comps, hide.ns = T, symnum.args = ymnum.args, 
                         tip.length = 0.01, size= 3)
  }
})

require(patchwork)
p<- wrap_plots(p, nrow = 1)

#Plot with "typical" Violin from Seurat- cannot run stat_compare_means with split.by
# p<- VlnPlot(sp.merged, paste0("UCD_", ucd_cols[1:3]), group.by = "orig.ident",
#             split.by = "spot_ident_spec", 
#             idents = c("sertoli_like_cc", "adeno_like"), y.max = 1, cols = tumor_cols)  &
#   theme(legend.position = "bottom") & 
#   scale_fill_manual(values = tumor_cols, labels = c("adenoma",  
#                                                "early\nclear cell", "advanced\nclear cell")) &
#   xlab("Sample") & ylab("Score")


p_leg<- get_legend(p[[1]])

ucd_titles<- str_replace_all(ucd_cols, "[.]", " ")
ucd_titles[2]<- "seminiferous tubule\nepithelial cell"

for(i in 1:length(p)){
  p[[i]]$labels$title<- ucd_titles[i]
}

plot<- plot_grid(plot_grid(NULL, p_leg, NULL, ncol=3) ,p & NoLegend(), 
                 nrow = 2, rel_heights = c(0.1, 1))

plot_dir<- "./Plots/Requested_Plots/"
pdf(file=paste0(plot_dir, "TumorTypes_Leydig_SeminTubule_Sertoli_VlnPlot_revised.pdf"), 
    width = 12, height=7)
plot
dev.off()


#### Var genes between tumors - heatmap ####
subset<- subset(sp.merged, idents = c("sertoli_like_cc", "adeno_like"))

DefaultAssay(subset)<- "Spatial"
subset<- NormalizeData(subset)
subset<- FindVariableFeatures(subset, nfeatures = 500)

var_genes<- VariableFeatures(subset)

# var_gene_meta<- subset@assays$Spatial@meta.data
# var_gene_meta<- var_gene_meta %>% filter(vf_vst_counts_variable == TRUE)

var_genes<- var_genes[var_genes %in% rownames(subset@assays$SCT$scale.data)]


DefaultAssay(subset)<- "SCT"
# DoHeatmap(subset, var_genes, group.by = "spot_ident_spec")
meta<- subset@meta.data
meta$spot_id<- rownames(meta)
meta<- meta %>%  arrange(spot_ident_spec)
meta<- meta %>% mutate(ident_label = case_when(spot_ident_spec == "adeno_like" ~ "adenoma",
                                               spot_ident_spec == "clear_cell" ~ "early clear cell",
                                               spot_ident_spec == "clear_cell_adv" ~ "advanced clear cell")) %>% 
                       mutate(mut_status = case_when(orig.ident == "OVA547" ~ "WT/Q72L",
                                             orig.ident == "OVA708" ~ "Q72L/Q72L",
                                             orig.ident == "OVA818" ~ "Q72L/Q72L"))

meta$mut_status<- factor(meta$mut_status, levels = c("WT/Q72L", "Q72L/Q72L"))
meta$ident_label<- factor(meta$ident_label, levels= c("adenoma", "early clear cell", "advanced clear cell"))
meta$sample_name<- factor(meta$sample_name, levels = c("HT", "HM AdCre Rep1", "HM AdCre Rep2"))


mat<- as.matrix(subset@assays$SCT$scale.data[var_genes,meta$spot_id])
# mat[mat>2.5]<- 2.5
# mat[mat< -2.5]<- -2.5
#scaled_mat = t(scale(t(mat)))
#mat[mat> 3]<- 3
#mat[mat< -3]<- -3

samp_col<- c("HT" = "olivedrab", "HM AdCre Rep1" = "maroon", "HM AdCre Rep2" = "darkslateblue")

mut_col<- c("WT/Q72L" = "honeydew4", "Q72L/Q72L" = "firebrick4")

col_fun<-  circlize::colorRamp2(c(-2.5,0,2.5), c("blue4", "khaki2", "firebrick2"))

label_cols<- tumor_cols
names(label_cols)<- unique(meta$ident_label)

column_ha = ComplexHeatmap::HeatmapAnnotation(type = meta$ident_label,
                                       sample= meta$sample_name,
                                       Rras2= meta$mut_status,
                       col = list(type = label_cols, sample= samp_col, 
                                  Rras2= mut_col))

ht<- ComplexHeatmap::Heatmap(mat, show_row_names = FALSE, show_column_names = FALSE, 
                        cluster_rows = T, cluster_columns = F, 
                        name= "scaled\nexpression",
                        row_dend_reorder = T, 
                        show_row_dend = T,  
                        col = col_fun, top_annotation = column_ha)

#repeat, but specifying two clusters 
HM <- ComplexHeatmap::Heatmap(mat, km=2,,show_row_names = FALSE, show_column_names = FALSE,
                              cluster_rows = T, cluster_columns = F, 
                              name= "scaled\nexpression",
                              row_dend_reorder = T, 
                              show_row_dend = T,  
                              col = col_fun, top_annotation = column_ha)  #Make a heatmap, and have 3 clusters
HM <- ComplexHeatmap::draw(HM)

#Extract genes for GSEA
r.dend <- ComplexHeatmap::row_dend(HM)  #Extract row dendrogram
rcl.list <- ComplexHeatmap::row_order(HM)

lapply(rcl.list, length)

heatmap_genes<- lapply(rcl.list, function(x){
  rownames(mat)[x]
})

names(heatmap_genes)<- c("Adeno", "Clear_cell")

saveRDS(heatmap_genes, file = "./Results/Tumor_500_variable_genes.rds")

#### Selected genes between tumors - heatmap ####

#first, markers for cc vs adeno 
cc_vs_adeno<- read.table("./Results/DEGTestResults/AllMutSamples_sertoli_cc_vs_adeno_markers.txt")
cc_vs_adeno$gene<- rownames(cc_vs_adeno)


cc_ad_genes<- cc_vs_adeno %>% filter(avg_log2FC < 0) %>% 
  slice_min(p_val_adj, n= 15) %>% arrange(avg_log2FC) %>% pull(gene)

cc_ad_genes<- c(cc_ad_genes, cc_vs_adeno %>% filter(avg_log2FC > 0) %>% 
  slice_min(p_val_adj, n= 15) %>% arrange(avg_log2FC) %>% pull(gene))

#then, markers for adv vs early cc
adv_vs_early_818<- read.table("./Results/DEGTestResults/OVA818_sertoli_cc_cluster5_vs_cluster0_markers.txt")
adv_vs_early_818$gene<- rownames(adv_vs_early_818)

adv_early_genes_818<- adv_vs_early_818 %>% filter(avg_log2FC < 0) %>% 
  slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% pull(gene)

adv_early_genes_818<- c(adv_early_genes_818, 
                        adv_vs_early_818 %>% filter(avg_log2FC > 0) %>% 
                          slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% pull(gene))

adv_vs_early_708<- read.table("./Results/DEGTestResults/OVA708_sertoli_cc_cluster9_vs_cluster2_markers.txt")
adv_vs_early_708$gene<- rownames(adv_vs_early_708)

adv_early_genes_708<- adv_vs_early_708 %>% filter(avg_log2FC < 0) %>% 
  slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% pull(gene)

adv_early_genes_708<- c(adv_early_genes_708, 
                        adv_vs_early_708 %>% filter(avg_log2FC > 0) %>% 
                          slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% pull(gene))

hm_genes<- unique(c(cc_ad_genes, adv_early_genes_708, adv_early_genes_818))

mat<- as.matrix(subset@assays$SCT$scale.data[hm_genes,meta$spot_id])

#try again, this time grouping by tumor type and sample (avg normalized expr and scale)
mat<- as.matrix(subset@assays$SCT$data[hm_genes,meta$spot_id]) %>% t() %>% as.data.frame()
mat[,c("orig.ident","spot_ident_spec")]<- meta[, c("orig.ident","spot_ident_spec")]

mat<- mat %>% group_by(orig.ident, spot_ident_spec) %>% summarise_at(hm_genes, mean)
mat_rows<- paste0(mat$orig.ident, "_", mat$spot_ident_spec)
mat<- mat %>% ungroup() %>% select(all_of(hm_genes)) 
mat<- as.matrix(mat)
mat<- scale(mat)
mat<- t(mat)
colnames(mat)<- mat_rows

meta_small<- meta %>% select(sample_name, ident_label, mut_status, orig.ident)
meta_small$names<- paste0(meta$orig.ident, "_", meta$spot_ident_spec)
meta_small<- unique(meta_small)
rownames(meta_small)<- meta_small$names 

meta_small<- meta_small[colnames(mat),]

meta_small$mut_status<- factor(meta_small$mut_status, levels = c("WT/Q72L", "Q72L/Q72L"))
meta_small$ident_label<- factor(meta_small$ident_label, levels= c("adenoma", "early clear cell", "advanced clear cell"))

order<- meta_small %>% arrange(ident_label) %>% rownames()

meta_small<- meta_small[order,]

column_ha = ComplexHeatmap::HeatmapAnnotation(type = meta_small$ident_label,
                                              sample= meta_small$sample_name,
                                              Rras2= meta_small$mut_status,
                                              col = list(type = label_cols, sample= samp_col, 
                                                         Rras2= mut_col))


ComplexHeatmap::Heatmap(mat[,order], 
                        show_row_names = TRUE, show_column_names = FALSE, 
                        cluster_rows = T, cluster_columns = F, 
                        name= "scaled\nexpression",
                        row_dend_reorder = T, 
                        show_row_dend = F,  
                        col = col_fun, top_annotation = column_ha)


#### Vln plot of clear cell genes in different tumor types ####

##### Upregulated advanced clear cell #####
g1<- adv_vs_early_708 %>% filter(avg_log2FC > 0) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)
g2<- adv_vs_early_818 %>% filter(avg_log2FC > 0) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)

genes<- unique(c(g1, g2, "Pax8")) 
  #18 genes total

max.val<- ceiling(max(sp.merged@assays$SCT@data[genes, Idents(sp.merged) %in% c("sertoli_like_cc")]))

#Modify comps to remove any adeno comparisons
comps2<- c(list(groups[3:4], groups[6:7]),
          lapply(c(3,6), function(x){ c(groups[1], groups[x])}))


p<- lapply( genes, function(gene){
    VlnPlot(sp.merged, gene, group.by =  "groups_compare",
            idents = c("sertoli_like_cc"), cols = tumor_cols_exp)  &
      scale_y_continuous(breaks = seq(0, max.val, 0.5)) &  # Only show ticks up to 1.0
      coord_cartesian(ylim = c(0, max.val + 1.2)) &  # Extend plotting space but keep tick labels fixed &
      theme(legend.position = "bottom") & 
      scale_x_discrete(labels= labels) &
      scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                        breaks = c("HM AdCre Rep1 adeno_like", "HT clear_cell",
                                   "HM AdCre Rep1 clear_cell_adv")) &
      guides(y = guide_axis(cap = "upper")) & 
      xlab("Sample") & ylab("Expression") &
      stat_compare_means(comparisons = comps2, hide.ns = T, symnum.args = ymnum.args, 
                         tip.length = 0.01, size= 3)
  }
)

#Plot in groups of 4
range <- 1:18
group_size <- 3

# Split into groups of 4
result <- split(range, ceiling(seq_along(range) / group_size))

p_leg<- get_legend(p[[1]])


plot_dir<- "./Plots/Requested_Plots/Tumors_Gene_Violin/AdvCCGenes/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

lapply(result, function(x){
  p1<- wrap_plots(p[x], nrow = 1)
  plot<- plot_grid(plot_grid(NULL, p_leg, NULL, ncol=3) ,p1 & NoLegend(), 
                   nrow = 2, rel_heights = c(0.1, 1))
  plot_name<- paste0("ClearCell_", str_flatten(genes[x], "_"), "_VlnPlot.pdf")
  pdf(file=paste0(plot_dir, plot_name), 
      width = 12, height=7)
  print(plot)
  dev.off()
})

##### Upregulated in adeno vs clear cell (all samples) #####

genes<- cc_vs_adeno %>% filter(avg_log2FC < 0) %>% 
  slice_min(p_val_adj, n= 18) %>% arrange(avg_log2FC) %>% pull(gene)


p<- lapply( genes, function(gene){
  VlnPlot(sp.merged, gene, group.by =  "groups_compare",
          idents = c("sertoli_like_cc", "adeno_like"), cols = tumor_cols_exp)  &
    scale_y_continuous(breaks = seq(0, max.val, 0.5)) &  # Only show ticks up to 1.0
    coord_cartesian(ylim = c(0, max.val + 1.2)) &  # Extend plotting space but keep tick labels fixed &
    theme(legend.position = "bottom") & 
    scale_x_discrete(labels= labels) &
    scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                      breaks = c("HM AdCre Rep1 adeno_like", "HT clear_cell",
                                 "HM AdCre Rep1 clear_cell_adv")) &
    guides(y = guide_axis(cap = "upper")) & 
    xlab("Sample") & ylab("Expression") &
    stat_compare_means(comparisons = comps, hide.ns = F, symnum.args = ymnum.args, 
                       tip.length = 0.01, size= 3)
}
)

#Plot in groups of 4
range <- 1:length(p)
group_size <- 3

# Split into groups of 4
result <- split(range, ceiling(seq_along(range) / group_size))

p_leg<- get_legend(p[[1]])

plot_dir<- "./Plots/Requested_Plots/Tumors_Gene_Violin/AdenoGenes/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

lapply(result, function(x){
  p1<- wrap_plots(p[x], nrow = 1)
  plot<- plot_grid(plot_grid(NULL, p_leg, NULL, ncol=3) ,p1 & NoLegend(), 
                   nrow = 2, rel_heights = c(0.1, 1))
  plot_name<- paste0("ClearCell_", str_flatten(genes[x], "_"), "_VlnPlot.pdf")
  pdf(file=paste0(plot_dir, plot_name), 
      width = 12, height=7)
  print(plot)
  dev.off()
})

##### Genes in rete + early/late cc + adeno #####

#Pax8, Wt1, Sox9, Rras2

#Expand the groups to include rete 
idents_keep<- c("sertoli_like_cc", "adeno_like", "rete", "rete1", "rete2")

#Re-configure b/c it's gotten messed up from plotting
sp.merged$groups_compare<- paste0(sp.merged$sample_name, " ", sp.merged$spot_ident_spec)

df<- sp.merged@meta.data %>% filter(spot_ident %in% idents_keep) 

labels<- df$sample_name
names(labels)<- df$groups_compare


#Factor the groups to have the order HT -> Adeno -> CC groups
sp.merged$groups_compare<- factor(sp.merged$groups_compare, 
                                  levels = c("HT normal_rete", "HT clear_cell", 
                                             "HM AdCre Rep1 normal_rete", "HM AdCre Rep1 cystic_rete", 
                                             "HM AdCre Rep1 adeno_like", "HM AdCre Rep1 clear_cell",    
                                             "HM AdCre Rep1 clear_cell_adv","HM AdCre Rep2 adeno_like", "HM AdCre Rep2 clear_cell",
                                             "HM AdCre Rep2 clear_cell_adv", " "))
#Specify the colors for the groups 
groups<- levels(sp.merged$groups_compare)

tumor_cols_exp<- rep("#8EC7FF", length(groups))

tumor_cols_exp[groups %>% endsWith(., "normal_rete")]<- "#EFC950"
tumor_cols_exp[groups %>% endsWith(., "cystic_rete")]<- "#00A86A"
tumor_cols_exp[groups %>% endsWith(., "clear_cell")]<- "#E6797C"
tumor_cols_exp[groups %>% endsWith(., "clear_cell_adv")]<- "#018C87"
names(tumor_cols_exp)<- groups

#Add labels only for one of each tumor group
labels_legend<-  c( "HT normal_rete"= "normal rete",
                    "HM AdCre Rep1 cystic_rete" = "cystic rete",
                    "HM AdCre Rep1 adeno_like" ="adenoma", #"HM AdCre Rep2 adeno_like" ="adenoma",
                    "HT clear_cell" = "early\nclear cell", #"HM AdCre Rep1 clear_cell" = "early\nclear cell", "HM AdCre Rep2 clear_cell" ="early\nclear cell",    
                    "HM AdCre Rep1 clear_cell_adv" ="advanced\nclear cell") #, "HM AdCre Rep2 clear_cell_adv"= "advanced\nclear cell")

#Make comparisons list for stat_compare_means
#Compare all early clear cell, within tumor types per sample
comps<- c(lapply(c(2), function(x){ c(groups[1], groups[x])}),
          lapply(c(4:10), function(x){ c(groups[3], groups[x])}),
          lapply(c(5:10), function(x){ c(groups[4], groups[x])}))

genes<- c("Pax8", "Wt1", "Sox9", "Rras2", "Sox8", "Dmrt1")

ht_norm_rete<- read.table("./Results/SampleSpotIdentMarkers/OVA547_rete_markers.txt")
ht_norm_rete$gene<- rownames(ht_norm_rete)
g1<- ht_norm_rete %>% filter(avg_log2FC > 0) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)

hm_norm_rete<- read.table("./Results/SampleSpotIdentMarkers/OVA708_rete1_markers.txt")
hm_norm_rete$gene<- rownames(hm_norm_rete)
g2<- hm_norm_rete %>% filter(avg_log2FC > 0) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)


hm_cys_rete<- read.table("./Results/SampleSpotIdentMarkers/OVA708_rete2_markers.txt")
hm_cys_rete$gene<- rownames(hm_cys_rete)
g3<- hm_cys_rete %>% filter(avg_log2FC > 0) %>% slice_min(p_val_adj, n= 10) %>% pull(gene)

genes<- unique(c(genes, g1, g2, g3))

p<- lapply( genes, function(gene){
  VlnPlot(sp.merged, gene, group.by =  "groups_compare",
          idents = idents_keep, cols = tumor_cols_exp)  &
    scale_y_continuous(breaks = seq(0, max.val, 0.5)) &  # Only show ticks up to 1.0
    coord_cartesian(ylim = c(0, max.val + 1.2)) &  # Extend plotting space but keep tick labels fixed &
    theme(legend.position = "bottom") & 
    scale_x_discrete(labels= labels) &
    scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                      breaks = names(labels_legend)) &
    guides(y = guide_axis(cap = "upper")) & 
    xlab("Sample") & ylab("Expression") &
    stat_compare_means(comparisons = comps, hide.ns = F, symnum.args = ymnum.args, 
                       tip.length = 0.01, size= 3)
}
)

#Plot in groups of 4
range <- 1:length(genes)
group_size <- 3

# Split into groups of 4
result <- split(range, ceiling(seq_along(range) / group_size))

p_leg<- get_legend(p[[1]])

plot_dir<- "./Plots/Requested_Plots/Tumors_Gene_Violin/ReteGenes/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

lapply(result, function(x){
  p1<- wrap_plots(p[x], nrow = 1)
  plot<- plot_grid(plot_grid(NULL, p_leg, NULL, ncol=3) ,p1 & NoLegend(), 
                   nrow = 2, rel_heights = c(0.1, 1))
  plot_name<- paste0("Rete_", str_flatten(genes[x], "_"), "_VlnPlot.pdf")
  pdf(file=paste0(plot_dir, plot_name), 
      width = 12, height=7)
  print(plot)
  dev.off()
})


##### Epithelial genes: 708 epithelial + early/late cc + adeno #####
#Wnt2b, Upk1b, Unc45b (top epithelial markers)

idents_keep<- c("sertoli_like_cc", "adeno_like", "ov_epithelial")

sp.merged$groups_compare<- paste0(sp.merged$sample_name, " ", sp.merged$spot_ident_spec)

df<- sp.merged@meta.data %>% filter(spot_ident %in% idents_keep) 

df$groups_compare<- paste0(df$sample_name, " ", df$spot_ident_spec)

labels<- df$sample_name
names(labels)<- df$groups_compare

#Factor the groups to have the order HT -> Adeno -> CC groups
sp.merged$groups_compare<- factor(sp.merged$groups_compare, 
                                  levels = c("HT ov_epithelial", "HT clear_cell",
                                             "HM AdCre Rep1 adeno_like", "HM AdCre Rep1 clear_cell",    
                                             "HM AdCre Rep1 clear_cell_adv","HM AdCre Rep2 adeno_like", "HM AdCre Rep2 clear_cell",
                                             "HM AdCre Rep2 clear_cell_adv", " "))
sp.merged$groups_compare[is.na(sp.merged$groups_compare)]<- " "
#Specify the colors for the groups 
groups<- levels(sp.merged$groups_compare)

tumor_cols_exp<- rep("#8EC7FF", length(groups))

tumor_cols_exp[groups %>% endsWith(., "ov_epithelial")]<- "#7570B3"
tumor_cols_exp[groups %>% endsWith(., "clear_cell")]<- "#E6797C"
tumor_cols_exp[groups %>% endsWith(., "clear_cell_adv")]<- "#018C87"
names(tumor_cols_exp)<- groups

#Add labels only for one of each tumor group
labels_legend<-  c( "HT ov_epithelial" = "ovarian epithelium",
                    "HM AdCre Rep1 adeno_like" ="adenoma", #"HM AdCre Rep2 adeno_like" ="adenoma",
                    "HT clear_cell" = "early\nclear cell", #"HM AdCre Rep1 clear_cell" = "early\nclear cell", "HM AdCre Rep2 clear_cell" ="early\nclear cell",    
                    "HM AdCre Rep1 clear_cell_adv" ="advanced\nclear cell") #, "HM AdCre Rep2 clear_cell_adv"= "advanced\nclear cell")

#Make comparisons list for stat_compare_means
#Compare all early clear cell, within tumor types per sample
comps<- c(lapply(c(2:8), function(x){ c(groups[1], groups[x])}),
          lapply(c(4:5), function(x){ c(groups[3], groups[x])}),
          lapply(c(7:8), function(x){ c(groups[6], groups[x])}))


epi<- read.table("./Results/SampleSpotIdentMarkers/OVA547_ov_epithelial_markers.txt")
epi$gene<- rownames(epi)


genes<- epi %>% filter(avg_log2FC > 0) %>% 
  slice_min(p_val_adj, n= 12) %>% pull(gene)


p<- lapply( genes, function(gene){
  VlnPlot(sp.merged, gene, group.by =  "groups_compare",
          idents = idents_keep, cols = tumor_cols_exp)  &
    scale_y_continuous(breaks = seq(0, max.val, 0.5)) &  # Only show ticks up to 1.0
    coord_cartesian(ylim = c(0, max.val + 1.2)) &  # Extend plotting space but keep tick labels fixed &
    theme(legend.position = "bottom") & 
    scale_x_discrete(labels= labels) &
    scale_fill_manual(values = tumor_cols_exp, labels = labels_legend, 
                      breaks = names(labels_legend)) &
    guides(y = guide_axis(cap = "upper")) & 
    xlab("Sample") & ylab("Expression") &
    stat_compare_means(comparisons = comps, hide.ns = F, symnum.args = ymnum.args, 
                       tip.length = 0.01, size= 3)
}
)

#Plot in groups of 4
range <- 1:length(genes)
group_size <- 3

# Split into groups of 4
result <- split(range, ceiling(seq_along(range) / group_size))

p_leg<- get_legend(p[[1]])

plot_dir<- "./Plots/Requested_Plots/Tumors_Gene_Violin/EpithelialGenes/"
ifelse(!dir.exists(plot_dir), dir.create(plot_dir, recursive = T), FALSE)

lapply(result, function(x){
  p1<- wrap_plots(p[x], nrow = 1)
  plot<- plot_grid(plot_grid(NULL, p_leg, NULL, ncol=3) ,p1 & NoLegend(), 
                     nrow = 2, rel_heights = c(0.1, 1))
  plot_name<- paste0("Epithelial_", str_flatten(genes[x], "_"), "_VlnPlot.pdf")
  pdf(file=paste0(plot_dir, plot_name), 
        width = 12, height=7)
  print(plot)
  dev.off()
})


#### DotPlot of atretic vs healthy follicle ####

#Edit identities and levels for plotting in the correct order
Idents(sp.merged)<- sp.merged$spot_ident

follicles<- c("follicle_healthy", "follicle_atretic")

ident_levels<- c(follicles, unique(sp.merged$spot_ident)[!unique(sp.merged$spot_ident) %in% follicles])

Idents(sp.merged)<- factor(Idents(sp.merged), rev(ident_levels))
sp.merged$orig.ident<- factor(sp.merged$orig.ident, samples)

sp.merged$sample_ident<- paste0(sp.merged$sample_name, "_", sp.merged$spot_ident)

sp.merged$sample_ident_label<- str_replace_all(sp.merged$sample_ident, "_", " ")

sample_ident_label<- unique(sp.merged$sample_ident_label)

Idents(sp.merged)<- factor(sp.merged$sample_ident_label, 
                           levels= c(sample_ident_label[!grepl("follicle", sample_ident_label)],
                                     sample_ident_label[grepl("follicle atretic", sample_ident_label)],
                                     sample_ident_label[grepl("follicle healthy", sample_ident_label)]))

idents_keep<- unique(sp.merged$sample_ident_label)[grep("follicle", unique(sp.merged$sample_ident_label))]
idents_keep

#Get the gene lists for plotting
atretic_DEG<- read.table("./Results/DEGTestResults/OVA547_atretic_v_healthy_follicle_markers.txt")
atretic_DEG$gene<- rownames(atretic_DEG)

atretic_genes<- atretic_DEG %>% filter(avg_log2FC < 0) %>% slice_min(p_val_adj, n= 10) %>% arrange(avg_log2FC) %>% pull(gene)

atretic_genes2<- atretic_DEG %>% filter(avg_log2FC < 0) %>% 
  slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% pull(gene)

atretic_genes2<- c(atretic_DEG %>% filter(avg_log2FC > 0) %>% 
                     slice_min(p_val_adj, n= 5) %>% arrange(avg_log2FC) %>% 
                     pull(gene), atretic_genes2)

genes<- list(atretic_genes, atretic_genes2)

plots<- lapply(genes, function(g){
  DotPlot(sp.merged,  idents = idents_keep, 
          features = g, cols = "Spectral", scale = T) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
})

#Get the min/max scaled values from both plots to keep scale consistent
min.val<- floor(min(unlist(lapply(plots, function(x){ min(x[["data"]]$avg.exp.scaled) }))))

max.val<- ceiling(max(unlist(lapply(plots, function(x){ max(x[["data"]]$avg.exp.scaled) }))))

meta<- left_join(meta, sp.merged@meta.data)

meta_simple<- meta[,c("sample_name","spot_ident_spec", "spot_ident", 
                      "ident_label", "mut_status")] %>% unique


mut_col<- c(mut_col, "WT/WT" = "ivory3")

plots2 <- lapply(plots, function(p){
  p$data[["Sample"]]<- str_split_fixed(p$data$id, " follicle", n=2)[,1]
  p$data <- left_join(p$data, meta_simple[,c("sample_name", "mut_status")], 
                      by= join_by(Sample == sample_name), relationship = "many-to-many")
  
  p$data$mut_status<- factor(p$data$mut_status, levels = c("WT/WT","WT/Q72L", "Q72L/Q72L"))
  p$data$mut_status[p$data$Sample == "WT"]<- "WT/WT"
  
  p<- p+ scale_colour_gradientn(limits=c(min.val, max.val), colours=SpatialColors(n=100))+
    geom_tile(aes(y=id, x = length(unique(atretic_genes))+1, fill = mut_status), 
                   position = "identity") + scale_fill_manual(values=mut_col)+ labs(fill= "Rras2") + 
    theme(legend.title = element_text(size = 12), legend.text = element_text(size = 11))
  return(p)
})


plot_dir<- "./Plots/Requested_Plots/"

pdf(file=paste0(plot_dir, "AtreticGenes_DotPlot_AllNeg_revised.pdf"), 
    width = 10, height=5)
print(plots2[[1]])
dev.off()

pdf(file=paste0(plot_dir, "AtreticGenes_DotPlot_PosAndNeg_revised.pdf"), 
    width = 10, height=5)
print(plots2[[2]])
dev.off()
