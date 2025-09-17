library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(cowplot)

setwd("/gpfs/scratch/bsc64/bsc64717/Bustelo_OV")
#setwd("/mnt/beegfs/kimbach/Bustelo_OV")

#Best to use merged object for this, as we want cross-sample results as well

sp.merged<- readRDS("./Robjs/AllSamp_seurat_regionsubset_merged.rds")


#Need to make samp_spot_ident and samp_cluster groupings

sp.merged@meta.data[,"samp_spot_ident"]<- paste0(sp.merged$orig.ident, "_", sp.merged$spot_ident)
sp.merged@meta.data[,"samp_cluster"]<- paste0(sp.merged$orig.ident, "_", sp.merged$seurat_clusters)

#set idents to samp_spot_ident
Idents(sp.merged)<- sp.merged$samp_spot_ident

ifelse(!dir.exists("./Results/DEGTestResults"),
       dir.create("./Results/DEGTestResults"), FALSE)
#### DEG testing ####
##### Atretic v healthy follicles ####
#look at genes different from atretic v healthy follicles 

atret_v_healthy_547<- FindMarkers(sp.merged, ident.1="OVA547_follicle_atretic",
ident.2="OVA547_follicle_healthy", min.pct=0.1)

atret_v_healthy_547$sample<- "OVA547_only"
atret_v_healthy_547$test <- "atretic_vs_healthy_follicle"

#Problem is there is stroma mixed in w atretic spots b/c of Visium resolution
#Omit stroma positive markers from upregulated genes in atretic follicles ? 

marks_stroma_547<- read.table(file= "./Results/SampleSpotIdentMarkers/OVA547_ov_stroma_markers.txt")

#get genes upreg in stroma
pos_stroma<- rownames(marks_stroma_547)[marks_stroma_547$avg_log2FC > 0 & marks_stroma_547$p_val_adj < 0.05]

#get atretic upreg genes
pos_atretic<- rownames(atret_v_healthy_547)[atret_v_healthy_547$avg_log2FC > 0 & atret_v_healthy_547$p_val_adj< 0.05]

length(pos_atretic)
#[1] 391
length(pos_stroma)
#[1] 5641
sum(pos_atretic %in% pos_stroma)
#[1] 281

overlap<- pos_atretic[pos_atretic %in% pos_stroma]

#add these to table 
atret_v_healthy_547[,"stroma_pos"]<- "NA"
atret_v_healthy_547[pos_atretic,"stroma_pos"]<- "No"
atret_v_healthy_547[overlap,"stroma_pos"]<- "Yes"

write.table(atret_v_healthy_547, file= "./Results/DEGTestResults/OVA547_atretic_v_healthy_follicle_markers.txt")

#Out of curiosity, look at DEGs between 522 & 547 healthy follicles
atret_547_v_healthy_522<- FindMarkers(sp.merged, ident.1="OVA547_follicle_atretic",
ident.2="OVA522_follicle_healthy", min.pct=0.1)

atret_547_v_healthy_522$sample<- "OVA547_OVA522"
atret_547_v_healthy_522$test <- "atretic_vs_healthy_follicle"

#Repeat the same as previously 
pos_atretic<- rownames(atret_547_v_healthy_522)[atret_547_v_healthy_522$avg_log2FC > 0 & atret_547_v_healthy_522$p_val_adj< 0.05]


length(pos_atretic)
#[1] 472
length(pos_stroma)
#[1] 5641
sum(pos_atretic %in% pos_stroma)
#[1] 314

overlap<- pos_atretic[pos_atretic %in% pos_stroma]

#add these to table 
atret_547_v_healthy_522[,"stroma_pos"]<- "NA"
atret_547_v_healthy_522[pos_atretic,"stroma_pos"]<- "No"
atret_547_v_healthy_522[overlap,"stroma_pos"]<- "Yes"

write.table(atret_547_v_healthy_522, file= "./Results/DEGTestResults/OVA547_atretic_v_OVA522_healthy_follicle_markers.txt")

atret1<- rownames(atret_v_healthy_547)[atret_v_healthy_547$stroma_pos == "No"]
atret2<- rownames(atret_547_v_healthy_522)[atret_547_v_healthy_522$stroma_pos == "No"]

sum(atret1 %in% atret2)

true_atret<- atret1[atret1 %in% atret2]

atret_547<- atret_v_healthy_547[true_atret,]
atret_522<- atret_547_v_healthy_522[true_atret,]

true_atret_results<- merge(atret_547, atret_522, by= "row.names", suffixes= c("_OVA547healthy", "_OVA522healthy"))

write.table(true_atret_results, file= "./Results/DEGTestResults/OVA547_atretic_consistent_pos_markers.txt")

##### Clear cell vs adenoma ####

#First, within double mutant samples 

#OVA 708 
cc_v_adeno_708<- FindMarkers(sp.merged, ident.1="OVA708_sertoli_like_cc",
ident.2="OVA708_adeno_like", min.pct=0.1)

cc_v_adeno_708$sample<- "OVA708_only"
cc_v_adeno_708$test <- "sertoli_like_cc_vs_adeno_like"

write.table(cc_v_adeno_708, file= "./Results/DEGTestResults/OVA708_sertoli_cc_vs_adeno_markers.txt")

#OVA 818 
cc_v_adeno_818<- FindMarkers(sp.merged, ident.1="OVA818_sertoli_like_cc",
ident.2="OVA818_adeno_like", min.pct=0.1)

cc_v_adeno_818$sample<- "OVA818_only"
cc_v_adeno_818$test <- "sertoli_like_cc_vs_adeno_like"

write.table(cc_v_adeno_818, file= "./Results/DEGTestResults/OVA818_sertoli_cc_vs_adeno_markers.txt")

#Look at overlapping genes in each cell compartment across samples 
#cc first
cc_708<- cc_v_adeno_708[cc_v_adeno_708$avg_log2FC > 0 & cc_v_adeno_708$p_val_adj < 0.05, ]
nrow(cc_708)
	#1443 genes

cc_818<- cc_v_adeno_818[cc_v_adeno_818$avg_log2FC > 0 & cc_v_adeno_818$p_val_adj < 0.05, ]
nrow(cc_818)
	#1436 genes

sum(rownames(cc_708) %in% rownames(cc_818))
	#1063 genes overlap 

#adeno
ad_708<- cc_v_adeno_708[cc_v_adeno_708$avg_log2FC < 0 & cc_v_adeno_708$p_val_adj < 0.05, ]
nrow(ad_708)
	#2291 genes

ad_818<- cc_v_adeno_818[cc_v_adeno_818$avg_log2FC < 0 & cc_v_adeno_818$p_val_adj < 0.05, ]
nrow(ad_818)
	#3291 genes

sum(rownames(ad_708) %in% rownames(ad_818))
	#1836 genes overlap 

#Keep all the results that overlap between the two 
cc_overlap<- rownames(cc_708)[rownames(cc_708) %in% rownames(cc_818)]
ad_overlap<- rownames(ad_708)[rownames(ad_708) %in% rownames(ad_818)]

overlap<- c(cc_overlap, ad_overlap)

true_cc_v_ad_results<- merge(cc_v_adeno_708[overlap,], cc_v_adeno_818[overlap,], by= "row.names", suffixes= c("_OVA708", "_OVA818"))

write.table(true_cc_v_ad_results, file= "./Results/DEGTestResults/OVA708_OVA818_sertoli_cc_vs_adeno_consistent_markers.txt")

#Look at granulosa genes in 818 vs sertoli and adenoma 

gran_v_malig_818<- FindMarkers(sp.merged, ident.1="OVA818_granulosa_tumor",
                               ident.2=c("OVA818_adeno_like" ,"OVA818_sertoli_like_cc"), min.pct=0.1)

gran_v_malig_818$sample<- "OVA818_only"
gran_v_malig_818$test <- "granulosa_tumor_vs_sertoli_like_cc_and_adeno_like"

write.table(gran_v_malig_818, file= "./Results/DEGTestResults/OVA818_granulosa_tumor_vs_sertoli_cc_and_adeno_markers.txt")

#Look at general clear cell vs adenoma across all samples 
Idents(sp.merged)<- sp.merged$spot_ident

cc_v_adeno_all<- FindMarkers(sp.merged, ident.1="sertoli_like_cc",
ident.2="adeno_like", min.pct=0.1)

cc_v_adeno_all$sample<- "OVA547_OVA708_OVA818"
cc_v_adeno_all$test <- "sertoli_like_cc_vs_adeno_like"

write.table(cc_v_adeno_all, file= "./Results/DEGTestResults/AllMutSamples_sertoli_cc_vs_adeno_markers.txt")

#Compare different clear cell clusters within double mutant samples 
Idents(sp.merged)<- sp.merged$samp_cluster

#OVA 708- Cluster 9 (more advanced) vs cluster 2
clust9_v_2_708<- FindMarkers(sp.merged, ident.1=,"OVA708_9",
ident.2="OVA708_2", min.pct=0.1)

clust9_v_2_708$sample<- "OVA708"
clust9_v_2_708$test <- "sertoli_like_cc_cluster9_vs_cluster2"

write.table(clust9_v_2_708, file= "./Results/DEGTestResults/OVA708_sertoli_cc_cluster9_vs_cluster2_markers.txt")

#OVA 818- Cluster 5 (more advanced) vs cluster 0
clust5_v_0_818<- FindMarkers(sp.merged, ident.1=,"OVA818_5",
                             ident.2="OVA818_0", min.pct=0.1)

clust5_v_0_818$sample<- "OVA818"
clust5_v_0_818$test <- "sertoli_like_cc_cluster5_vs_cluster0"

write.table(clust5_v_0_818, file= "./Results/DEGTestResults/OVA818_sertoli_cc_cluster5_vs_cluster0_markers.txt")

##### Rete structures in 708 #####
Idents(sp.merged)<- sp.merged$samp_spot_ident

rete_708<- FindMarkers(sp.merged, ident.1="OVA708_rete1",
                       ident.2="OVA708_rete2", min.pct=0.1)

rete_708$sample<- "OVA708"
rete_708$test <- "rete1_vs_rete2"

write.table(rete_708, file= "./Results/DEGTestResults/OVA708_rete1_vs_rete2_markers.txt")

##### 708 cc rete vs actual clear cell area #####

rete_v_cc_708<- FindMarkers(sp.merged, ident.1="OVA708_rete2",
                       ident.2="OVA708_sertoli_like_cc", min.pct=0.1)

rete_v_cc_708$sample<- "OVA708"
rete_v_cc_708$test <- "rete2_vs_sertoli_like_cc"

write.table(rete_v_cc_708, file= "./Results/DEGTestResults/OVA708_rete2_vs_sertoli_like_cc_markers.txt")


##### 547 rete vs 708 rete1 #####
rete_547_v_rete_708<- FindMarkers(sp.merged, ident.1="OVA547_rete",
                            ident.2="OVA708_rete1", min.pct=0.1)

rete_547_v_rete_708$sample<- "OVA547_OVA708"
rete_547_v_rete_708$test <- "OVA547_rete_vs_OVA708_rete1"

write.table(rete_547_v_rete_708, file= "./Results/DEGTestResults/OVA547_rete_vs_OVA708_rete1_markers.txt")


