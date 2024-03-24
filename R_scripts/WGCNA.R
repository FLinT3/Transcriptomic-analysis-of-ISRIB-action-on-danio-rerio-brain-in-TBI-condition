#####--------------- Load packages ---------------#####
library(WGCNA)
library(clusterProfiler)

#####--------------- WGCNA ---------------#####
t_rld <- t(assay(rld))

gsg <- goodSamplesGenes(t_rld, verbose = 4);

# if TRUE, then all OK. If not - filter FALSE
gsg$allOK

# Then a set of potential degrees for network analysis is determined and a soft threshold is searched using pickSoftThreshold, and the results are visualized.powers = c(c(1:15), seq(from = 16, to=30, by=2))
sft = pickSoftThreshold(t_rld, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

picked_power = 10

# Network creation
netKOF = blockwiseModules(t_rld, power = picked_power, networkType	= "unsigned",
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE, maxBlockSize = 4000,
                          saveTOMFileBase = "TBI", randomSeed = 999,
                          verbose = 3, loadTOM = TRUE)

module_df <- data.frame(
  gene_id = names(netKOF$colors),
  colors = labels2colors(netKOF$colors)
)

# Convert labels to colors for plotting
mergedColors = labels2colors(netKOF$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netKOF$dendrograms[[1]],
  mergedColors[netKOF$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

plotDendroAndColors(
  netKOF$dendrograms[[2]],
  mergedColors[netKOF$blockGenes[[2]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

plotDendroAndColors(
  netKOF$dendrograms[[3]],
  mergedColors[netKOF$blockGenes[[3]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

plotDendroAndColors(
  netKOF$dendrograms[[4]],
  mergedColors[netKOF$blockGenes[[4]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

plotDendroAndColors(
  netKOF$dendrograms[[5]],
  mergedColors[netKOF$blockGenes[[5]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)


sizeGrWindow(12, 9)
mergedColors = labels2colors(netKOF$colors)

moduleLabels = netKOF$colors
moduleColors = labels2colors(netKOF$colors)
MEs = netKOF$MEs;
geneTree = netKOF$dendrograms[[1]];

MEs0 = moduleEigengenes(t_rld, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# 


moduleTraitCor = cor(MEs, MetaDATA[-c(1:4)], use = "p", method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MetaDATA));
# moduleTraitPvalue[moduleTraitPvalue < 0.05] = NA
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


moduleTraitCor = cor(MEs, MEs, use = "p", method = "spearman");
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor, nrow(MetaDATA[-c(1:3)]));
moduleTraitPvalue[moduleTraitPvalue1 < 0.05] = NA
textMatrix1 = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1) = dim(moduleTraitCor)


# Correlation between modules and samples
par(mar = c(13, 9.5, 3, 3));
png(filename = "WGCNA_cor.png", width = 6000, height = 6000, res = 400)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = substring(names(MEs), 3),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

par(mar = c(8, 9.5, 3, 3));
png(filename = "WGCNA_cor_moduls.png", width = 6000, height = 6000, res = 400)
par(mar = c(6, 8.5, 3, 3));

# Module correlation
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               xSymbols = substring(names(MEs), 3),
               yLabels = names(MEs),
               ySymbols = substring(names(MEs), 3),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix1,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Module relationships"))
dev.off()

# Examine Expression Profiles
# pick out a few modules of interest here
modules_of_interest = c("turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized <- assay(rld)

subexpr = expr_normalized[submod$gene_id,]


submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

# Visualize gene expression in cluster between samples
submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

dev.off()

##### -----Over represented analysis----- #####
# Cluster extraction
color_list <- unique(module_df$colors)

# Directories creation
sapply(color_list, function(x) dir.create(paste0(path = "C:/Users/anton/Desktop/RNA-seq_TBI_v-2/WGCNA/modules/", x)))

auto_enrichment <- function(module){
  
  genes_of_interest = module_df %>% subset(colors %in% module)
  print(head(genes_of_interest))
  
  setwd(paste0("C:/Users/anton/Desktop/RNA-seq_TBI_v-2/WGCNA/modules/", module, "/"))
  # main analysis
  expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
  # print(head(expr_of_interest))
  
  genes <- as.data.frame(expr_of_interest)
  
  # Entrez gene ID
  ensemble <- rownames(genes)
  ENTREZ <- mapIds(org.Dr.eg.db, keys=ensemble, keytype="ENSEMBL", column="ENTREZID")
  
  print("Everything OK")
  
  ### GO over-representation analysis
  print("GO")
  tryCatch({
  ego_bp <- enrichGO(gene          = ensemble,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  ego1_bp <- as.data.frame(ego_bp)
  
  if (nrow(ego1_bp) > 0){
  
    p1 <- dotplot(ego_bp, showCategory=20, label_format = 40, font.size = 10) + ggtitle("dotplot for GO BP")
    ggsave("BP.png", plot = p1, device = "png")
    }
  }, error=function(e){})
  
  tryCatch({
  ego_mf <- enrichGO(gene          = ensemble,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  ego1_mf <- as.data.frame(ego_mf)
  
  if (nrow(ego1_mf) > 0){
    p2 <- dotplot(ego_mf, showCategory=20, label_format = 40, font.size = 10) + ggtitle("dotplot for GO MF")
    ggsave("MF.png", plot = p2, device = "png")}
  }, error=function(e){})
  
  tryCatch({
  ego_cc <- enrichGO(gene          = ensemble,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  p3 <- dotplot(ego_cc, showCategory=20, label_format = 40, font.size = 10) + ggtitle("dotplot for GO CC")
  ggsave("CC.png", plot = p3, device = "png")
  
  }, error=function(e){})
  
  
  
  ### KEGG pathway over-representation analysis
  tryCatch({
  kk <- enrichKEGG(gene         = ENTREZ,
                   organism     = 'dre',
                   pvalueCutoff = 0.05)
  
  p4 <- dotplot(kk, showCategory=20, label_format = 40, font.size = 10) + ggtitle("dotplot for KEGG")
  ggsave("KEGG.png", plot = p4, device = "png")
  
  }, error=function(e){})
  print("next")
}

sapply(color_list, function(x) auto_enrichment(x))



