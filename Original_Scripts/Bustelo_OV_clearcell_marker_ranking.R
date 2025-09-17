library(dplyr)
library(RankProd)

# Load the data of Golub et al. (1999). data(golub) 
# contains a 3051x38 gene expression
# matrix called golub, a vector of length called golub.cl 
# that consists of the 38 class labels,
# and a matrix called golub.gnames whose third column 
# contains the gene names.
data(golub)
subset <- c(1:4,28:30)
RP.out <- RP(golub[,subset],golub.cl[subset],rand=123) 

# class 2: label =1, class 1: label = 0
#pfp for identifying genes that are up-regulated in class 2 
#pfp for identifying genes that are down-regulated in class 2 
head(RP.out$pfp)

#For one class data, the label for each sample should be 1

#Get the clear cell marker genes for each sample 
setwd("~/Documents/Spatial/Bustelo_OV")
samps<- c("OVA547", "OVA708", "OVA818")

tbl<- lapply(samps, function(samp){
  x<- read.table(paste0("Results/SampleSpotIdentMarkers/", samp, "_sertoli_like_cc_markers.txt"))
  x<- x %>% filter(p_val_adj < 0.001)
  colnames(x)<- paste0(samp, "_", colnames(x))
  x$gene<- rownames(x)
  return(x)
  }
)

tbl<- tbl %>%
  Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by="gene"), .)

lfc_mat<- tbl %>% select(ends_with("avg_log2FC")) %>% as.matrix()
gnames<- tbl$gene
cl<- rep(1, ncol(lfc_mat))

RP.test <- RankProducts(data= lfc_mat,cl = cl, gene.names = gnames,rand=123, logged = T)

RPrank<- data.frame(Rank_product_rank = RP.test$RPrank[,2], gene= rownames(RP.test$RPrank), 
                    Rank_product_pval = RP.test$pval[,2])

out<- merge(tbl, RPrank, by= "gene")

membr_genes<- read.csv("./MEMBRANE.v2024.1.Hs.tsv", sep= "\t", skip = 16)[1,]
membr_genes<- membr_genes[,2] %>% strsplit(",") %>% unlist()

library(nichenetr)
out$Human_gene<- nichenetr::convert_mouse_to_human_symbols(out$gene)

out[,"Is_membrane_gene"]<- out$Human_gene %in% membr_genes


out<- out %>% relocate(Rank_product_rank, Rank_product_pval,  Human_gene, Is_membrane_gene, .after = gene)

write.table(out, "./Results/clearcell_genes_ranked.txt")
