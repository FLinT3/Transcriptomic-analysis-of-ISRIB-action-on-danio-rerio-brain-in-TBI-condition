#####--------------- Load packages ---------------#####
library(ggplot2)
library(pheatmap)

#####--------------- Heatmap ---------------#####
select_genes <- dif_exp_genes$Geneid

mat <- assay(rld)[select_genes, ]
mat <- mat - rowMeans(mat)

Sy <- mapIds(org.Dr.eg.db, keys=rownames(mat), keytype="ENSEMBL", column="SYMBOL")
Sy <- as.character(Sy)
rownames(mat) <- Sy

# remove NA
mat <- mat[!is.na(rownames(mat)),]

ann_colors <- list(
  groups = c(C = "#3D4EE1", T = "#D52B18", TI = "#E88A0E", I = "#13A109"),
  condition_T = c("1" = "#E86DAF", "0" = "#15EBDC"),
  condition_I = c("1" = "#93E55D", "0" ="#DFF58B")
)

anno <- as.data.frame(colData(rld)[c("groups", "condition_T", "condition_I")])

# visualization
pheatmap(mat, annotation_colors = ann_colors, annotation_col = anno, fontsize = 6.9)
