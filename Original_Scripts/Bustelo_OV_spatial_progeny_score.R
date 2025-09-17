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

Idents(data)<- data$spot_ident

## ----"progeny", message=FALSE-----------------------------------------------------------------------------------------
net <- get_progeny(organism = 'mouse', top = 100)
head(net)

## ----"wmean", message=FALSE-------------------------------------------------------------------------------------------
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$SCT@data)

# Run wmean
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', minsize = 5)
head(acts)

## ----"new_assay", message=FALSE---------------------------------------------------------------------------------------
# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

# Save object with assay 
saveRDS(data, file= paste0("./Robjs/", samp,"_seurat_unfiltered.rds") )

## ----"projected_acts", message = FALSE, warning = FALSE, fig.width = 8, fig.height = 4--------------------------------
p1<- SpatialDimPlot(data, alpha=0) & NoLegend()

p2 <- SpatialFeaturePlot(data, features = c("Androgen", "Estrogen"), pt.size.factor = pt) & 
         scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') 

if(!dir.exists("./Plots/Spatial/Progeny_Results") dir.create("./Plots/Spatial/Progeny_Results")


pdf(file= paste0("./Plots/Spatial/Progeny_Results/", samp, "_Progeny_Andro_Estro_SpatialFeature.pdf"), width= 12, height=4)
plot_grid(p1, p2, nrow=1)
dev.off()

#Repeat, but for MAPK(ERK) and PI3K(AKT) pathways mentioned in paper
p2 <- SpatialFeaturePlot(data, features = c("MAPK", "PI3K"), pt.size.factor = pt) & 
         scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')

pdf(file= paste0("./Plots/Spatial/Progeny_Results/", samp, "_Progeny_MAPK_PI3K_SpatialFeature.pdf"), width= 12, height=4)
plot_grid(p1, p2, nrow=1)
dev.off()         

## ----"mean_acts", message = FALSE, warning = FALSE--------------------------------------------------------------------
# Extract activities from object as a long dataframe

progeny<- data@assays$pathwaysmlm@data

ifelse(!dir.exists("./Results/Progeny_Results"),
       dir.create("./Results/Progeny_Results"), FALSE)

write.table(progeny, file= paste0("./Results/Progeny_Results/", samp, "_full_progeny_df.txt"))

#keep only regions of interest
barcodes<- colnames(data)[data$region_of_interest == "TRUE"]

df <- t(as.matrix(data@assays$pathwaysmlm@data[,barcodes])) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)[barcodes]) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
pdf(file= paste0("./Plots/Spatial/Progeny_Results/", samp, "_regionofinterest_pheatmap.pdf"), width= 8, height=6)
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()

