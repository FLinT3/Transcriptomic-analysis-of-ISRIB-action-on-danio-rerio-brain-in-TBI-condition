#####--------------- Load packages ---------------#####
library(ggplot2)
library(ggrepel)

#####--------------- PCA Plot ---------------#####
select_genes <- rbind(head(re_1[order(re_1$log2FoldChange),], 250),
                      head(re_1[order(re_1$log2FoldChange, decreasing = TRUE),], 250)) %>% select(Geneid)
rld_1 <- rld[rownames(assay(rld)) %in% select_genes$Geneid]

pcaData <- plotPCA(rld_1, intgroup = c("groups", "condition_T", "condition_I"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = groups)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
