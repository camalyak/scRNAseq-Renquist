library(Seurat)
library(sceasy)
library(reticulate)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(GEOquery)
library(scCustomize)
library(openxlsx)
library(RColorBrewer)


#set up environment
use_python("C:/Users/kaylamac/AppData/Local/Programs/Python/Python312/python.exe")


#get file to read the data
sceasy::convertFormat("C:/Users/kaylamac/Downloads/liver_3m_4m_6m.h5ad", from="anndata", to="seurat",  outFile='EC_liver.rds')
obj <- readRDS('EC_liver.rds')


#exclude the reverse diet, get only the chow and western diet
liver <- subset(x = obj, subset = diet == 'chow' | diet == 'Western')


#genes we care about
genes_of_interest <- c("Abat", "Aldh5a1", "Slc6a6", "Slc6a8", "Slc6a12", "Slc6a13", "Piezo1", "Flt4", "Flt1", "Vegfa", "Vegfb", "Vegfc", "Vegfd")


#get landmark genes for each zone
markers <- read.csv("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Zoning/all_markers_liver_ec.csv")


#rename columns
names(markers) <- c("gene", "zone", "control_p_val", "control_avg_log2FC", "control_pct.1", "control_pct.2", "control_p_val_adj", "max_pval", "min_pval")


#remove the first two rows
markers <- markers[-c(1, 2), ]


#gives us the genes that belong to each zone
genes_by_zone <- split(markers$gene, markers$zone)


#subset for only the marker genes
all_genes <- rownames(liver)
mark_genes <- intersect(all_genes, markers$gene)
liver_sub <- subset(liver, features = mark_genes)


#run normal processing
liver_sub <- ScaleData(liver_sub)
liver_sub <- FindVariableFeatures(liver_sub)

set.seed(402)
liver_sub <- RunPCA(liver_sub)


#gives us top 25 genes per zone
top_genes_by_zone <- lapply(genes_by_zone, function(genes) head(genes, 25))


#manual scoring of cells based on top genes in each zone
for (zone in names(top_genes_by_zone)) {
  zone_genes <- intersect(top_genes_by_zone[[zone]], rownames(liver_sub))
  
  #only calc scores if there are genes available
  if (length(zone_genes) > 0) {
    #calc the average expression for each cell
    avg_expr <- Matrix::colMeans(GetAssayData(liver_sub, slot = "data")[zone_genes, , drop = FALSE])
    
    #add the score to the metadata according to the average expression
    score_name <- paste0("Zone ", zone)
    liver_sub[[score_name]] <- avg_expr
  } else {
    cat("No genes found for zone:", zone, "/n")
  }
}


#get scores
zone_scores <- FetchData(liver_sub, vars = grep("Zone ", colnames(liver_sub@meta.data), value = TRUE))


#assign each cell to the zone with the highest score
liver_sub$zone <- apply(zone_scores, 1, function(x) {
  zones <- names(x)
  zones[which.max(x)]
})


#check the distribution of cells across zones
zone_distribution <- table(liver_sub$zone)
print(zone_distribution)


FeaturePlot(liver_sub, features = colnames(zone_scores), reduction = "pca")


combined_genes <- unique(c(mark_genes, genes_of_interest))
liver_combined <- subset(liver, features = combined_genes)

liver_combined$zone <- liver_sub$zone[colnames(liver_combined)]

#plot the zones
DotPlot_scCustom(seurat_object = liver_combined, features = genes_of_interest, x_lab_rotate = TRUE, group.by = "zone", 
                 split.by = "diet",  dot.scale = 8)
DotPlot_scCustom()(liver_combined, 
        features = genes_of_interest, 
        group.by = "zone",
        split.by = "diet",
        dot.scale = 8)


saveRDS(liver_combined, file = "C:/Users/kaylamac/Documents/Kayla's Code Stuff/Zoning/Liver_Zoning_Git.rds")
obj <- readRDS("C:/Users/kaylamac/Documents/EC_liver.rds")

DotPlot_scCustom(obj, features = genes_of_interest, x_lab_rotate = TRUE, group.by = "cell_type", 
                 split.by = "diet",  dot.scale = 8)