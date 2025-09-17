library(dplyr)
library(Seurat)
library(ggplot2)

setwd("/your/parent/dir/")

anno_path<- "./Input_Files/AllSamp_spot_annotations.txt"

anno_df<- read.table(anno_path)

files<- list.files("./Robjs", pattern = "_seurat_unfiltered.rds", full.names = T)

objs<- lapply(files, readRDS)

if(!dir.exists("./Results/SampleSpotIdentMarkers") dir.create("./Results/SampleSpotIdentMarkers")

for (i in 1:length(objs)){
        seurat<- objs[[i]]
        samp<- unique(seurat$orig.ident)
        anno_sub<- anno_df %>% filter(orig.ident == samp)
        seurat@meta.data[anno_sub$barcode, anno_cols]<- anno_sub[,anno_cols]
        Idents(seurat)<- seurat$spot_ident
        annos<- unique(seurat$spot_ident)
        annos<- annos[annos != "Other"]
        for (a in annos){
                file<- paste0("./Results/SampleSpotIdentMarkers/", samp, "_", a, "_markers.txt")
                #if(file.exists(file)){
                #       next
                #}else{
                marks<-  FindMarkers(seurat, ident.1=a, min.pct=0.1)
                write.table(marks, file= file)
                #}
        }
}