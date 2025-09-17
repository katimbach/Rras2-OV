library(Seurat)
library(stringr)

setwd("~/Documents/Spatial/Bustelo_OV")

robjs<- list.files("./Robjs", full.names = T)

obj_list<- lapply(robjs, readRDS)

samp_names<- str_split_fixed(list.files("./Robjs"), "_", n=3)[,1]

dir.create("./UCD_input")

for(i in 1:length(obj_list)){
  sample<- samp_names[i]
  write.csv(as.matrix(obj_list[[i]]@assays$SCT@data), paste0("./UCD_input/", sample, "_data_input.csv"))
  cells.keep<- colnames(obj_list[[i]])
  write.csv(obj_list[[i]]@images[[1]]@coordinates[cells.keep,], paste0("./UCD_input/", sample, "_coords_input.csv"))
}

