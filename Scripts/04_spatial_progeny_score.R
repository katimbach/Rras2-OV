library(Seurat)
library(decoupleR)

# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(cowplot)

setwd("/your/parent/dir/")

args <- commandArgs(trailingOnly = TRUE)

samp<- args[1]
print(paste0("the sample in this script is ", samp))

pts<- c(3, 2, 2.5, 1.2)

samps<- read.table("./samples.txt", header=F)
pt<- pts[which(samps$V1 == samp)]

data <- readRDS(paste0("./Robjs/", samp,"_seurat_unfiltered.rds"))

net <- get_progeny(organism = 'mouse', top = 100)
head(net)

# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$SCT@data)

# Run wmean
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', minsize = 5)
head(acts)

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"


progeny<- data@assays$pathwaysmlm@data

ifelse(!dir.exists("./Results/Progeny_Results"),
       dir.create("./Results/Progeny_Results"), FALSE)

write.table(progeny, file= paste0("./Results/Progeny_Results/", samp, "_full_progeny_df.txt"))
