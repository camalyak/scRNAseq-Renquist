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
library(tibble)
library(ggrepel)


#virtual environment
use_python("C:/Users/kaylamac/AppData/Local/Programs/Python/Python312/python.exe")


#genes we care about
genes_of_interest <- c("Abat", "Aldh5a1", "Slc6a6", "Slc6a8", "Slc6a12", "Slc6a13", "Piezo1", "Flt1", "Flt4",
                       "Vegfa", "Vegfb", "Vegfc", "Vegfd")


#another way to read in matrix (easier)
nafld_ct_matrix <- Read10X("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/rawData_mouseNafld/countTable_mouseNafld", gene.column = 1)


#metadata for all cells
nafld_metadata <- read.csv("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/annot_mouseNafldAll.csv")
rownames(nafld_metadata) <- nafld_metadata$cell


#create seurat object that includes fibroblasts
nafld_all <- CreateSeuratObject(counts = nafld_ct_matrix, meta.data = nafld_metadata)


#save the seurat object of all cells
#31053 genes, 256345 cells
saveRDS(nafld_all, file = "C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all.rds")


#load in object
nafld_all <- readRDS("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all.rds")


#keep cells that have annotation (i.e. get rid of doublets and unannotated cells)
cells_to_keep <- rownames(nafld_all@meta.data)[!is.na(nafld_all@meta.data$annot)] 
nafld_all <- subset(nafld_all, cells = cells_to_keep)


#processing
nafld_all <- FindVariableFeatures(nafld_all)
nafld_all <- NormalizeData(nafld_all)
nafld_all <- ScaleData(nafld_all, features = VariableFeatures(nafld_all))


#get rid of rownames so it doesn't mess up when adding in umap coordinates
rownames(nafld_metadata) <- NULL


#plot umap according to coordiantes
umap_coords <- nafld_metadata %>%
  select(cell, UMAP_1, UMAP_2) %>%
  column_to_rownames(var = "cell")


umap_coords <- umap_coords[Cells(nafld_all), ]


#create dim reduction object -- umap
nafld_all@reductions$umap <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords),
  assay = "RNA",
  key = "UMAP_"
)


#change stromal cells to fibroblasts
nafld_all$annot <- ifelse(
  nafld_all$annot %in% c("Stromal cells"),
  "Fibroblasts",
  nafld_all$annot
)


#save the seurat object of filtered cells
#31053 genes, 121980 cells
saveRDS(nafld_all, file = "C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all.rds")


#load in filtered object
nafld_all <- readRDS("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all.rds")


#dimplot for all cells including immune cells
DimPlot(nafld_all, reduction = "umap", group.by = "annot", label = TRUE) + ggtitle("NAFLD All Cells")


#subset for cell types we are interested in (no immune cells)
nafld_all_sub <- subset(nafld_all, subset = annot %in% c("Hepatocytes", "Endothelial cells", "KCs", "Cholangiocytes", "Fibroblasts"))


#dimplot -- NO immune cells
png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_sub_dimplot.png", width = 2000, height = 1500, res = 300)
DimPlot(nafld_all_sub, reduction = "umap", group.by = "annot", label = TRUE) + ggtitle("Liver Atlas - NAFLD")
dev.off()


#create new metadata column that combines cell type and diet
nafld_all_sub$celltype_diet <- paste(
  nafld_all_sub$annot,
  nafld_all_sub$diet,
  sep = "_"
)


#dimplot by diet
png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_sub_diet_dimplot.png", width = 2000, height = 1500, res = 300)
DimPlot(nafld_all_sub, reduction = "umap", group.by = "diet", label = FALSE) + ggtitle("Liver Atlas - NAFLD by Diet")
dev.off()


#dotplot by cell type
png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_dotplot.png", width = 3000, height = 1000, res = 300)
DotPlot_scCustom(nafld_all_sub,
        features = genes_of_interest,
        group.by = "annot",
        x_lab_rotate = TRUE,
        colors_use = c("red", "green"),
        dot.min = 0.01,
        dot.scale = 10) + ggtitle("Mouse - Relative Cell Expression")
dev.off()


#dotplot by cell type and diet
png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_diet_dotplot.png", width = 3000, height = 1500, res = 300)
DotPlot_scCustom(nafld_all_sub,
        features = genes_of_interest,
        group.by = "celltype_diet",
        x_lab_rotate = TRUE,
        colors_use = c("red", "green"),
        dot.min = 0.01,
        dot.scale = 10) + ggtitle("Mouse - Relative Cell Expression by Diet")
dev.off()


#dotplot by cell type and diet -- just first 6 genes
DotPlot_scCustom(nafld_all_sub,
        features = genes_of_interest[1:6],
        group.by = "celltype_diet",
        x_lab_rotate = TRUE,
        colors_use = c("red", "blue"),
        dot.min = 0.01,
        dot.scale = 10) + ggtitle("Mouse - Relative Cell Expression by Diet")


#violin plot by cell type and diet -- just first 6 genes
VlnPlot(nafld_all_sub, features = genes_of_interest[1:6], group.by = "celltype_diet", ncol = 3)


#feature plot for all genes of interest
png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_featureplot.png", width = 6000, height = 3000, res = 300)
FeaturePlot(nafld_all_sub, features = genes_of_interest, ncol = 5)
dev.off()



#want to be specific with which fibroblasts are stellate cells
#need to add in meta data from other file

fibro_meta <- read.csv("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/annot_mouseNafldFibro.csv")

fibro_ct_matrix <- Read10X("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/rawData_mouseNafld/countTable_mouseNafld", gene.column = 1)


nafld_all_sub$annot.fibro <- ifelse(
  nafld_all_sub$annot == "Fibroblasts" & nafld_all_sub$cell %in% fibro_meta$cell,
  fibro_meta$annot[match(nafld_all_sub$cell, fibro_meta$cell)],
  nafld_all_sub$annot 
)

nafld_all_sub$annot.fibro[nafld_all_sub$annot.fibro == "Fibroblasts"] <- "Fibro"

saveRDS(nafld_all_sub, file = "C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_sub_fibro.rds")

View(nafld_all_sub@meta.data)

png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_fibro_dotplot.png", width = 3000, height = 1500, res = 300)
DotPlot_scCustom(nafld_all_sub,
        features = genes_of_interest,
        group.by = "annot.fibro",
        x_lab_rotate = TRUE,
        colors_use = c("red", "green"),
        dot.min = 0.01,
        dot.scale = 10) + ggtitle("Mouse - Relative Cell Expression")
dev.off()


nafld_all_sub$annot.fibro_diet <- paste(
  nafld_all_sub$annot.fibro,
  nafld_all_sub$diet,
  sep = "_"
)


png("C:/Users/kaylamac/Documents/Kayla's Code Stuff/Liver Atlas/NAFLD mouse/nafld_all_fibro_diet_dotplot.png", width = 3000, height = 1900, res = 300)
DotPlot_scCustom(nafld_all_sub,
        features = genes_of_interest,
        group.by = "annot.fibro_diet",
        x_lab_rotate = TRUE,
        colors_use = c("red", "green"),
        dot.min = 0.01,
        dot.scale = 10) + ggtitle("Mouse - Relative Cell Expression")
dev.off()