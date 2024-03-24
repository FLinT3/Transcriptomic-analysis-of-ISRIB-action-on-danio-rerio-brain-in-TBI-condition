#####--------------- Load packages ---------------#####
library(ggplot2)
library(ggrepel)

#####--------------- Volcano plot---------------#####
re_1$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
re_1$diffexpressed[re_1$log2FoldChange > 1 & re_1$padj < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
re_1$diffexpressed[re_1$log2FoldChange < -1 & re_1$padj < 0.05] <- "DOWN"

re_1$delabel <- NA
re_1$delabel[re_1$diffexpressed != "NO"] <- re_1$Symbol[re_1$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=re_1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
