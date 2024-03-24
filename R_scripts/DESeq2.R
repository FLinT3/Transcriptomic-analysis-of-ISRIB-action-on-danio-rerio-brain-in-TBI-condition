#####--------------- Load packages ---------------#####
library(AnnotationDbi)
library(org.Dr.eg.db)
library(Gviz)
library(DESeq2)
library(vsn) 

#####--------------- DESeq2 ---------------#####
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                design = ~ condition_T + condition_I + condition_T:condition_I)

# Counts filtering
keep <- rowSums(counts(dds) >= 10) > ncol(dds) * 0.3
dds <- dds[keep,]
DDS <- DESeq(dds)

result <- results(DDS, contrast = c("condition_T", 1, 0), alpha = 0.05)
summary(result)

result_1 <- results(DDS, contrast = list("condition_T_1_vs_0"), alpha = 0.05)
res_1 <- result_1[order(result_1$padj),]
summary(result_1)

result_2 <- results(DDS, contrast = list("condition_I_1_vs_0"), alpha = 0.05)
res_2 <- result_2[order(result_2$padj),]
summary(result_2)

result_3 <- results(DDS, contrast = list("condition_T1.condition_I1"), alpha = 0.05)
res_3 <- result_3[order(result_3$padj),]
summary(result_3)

res <- result[order(result),]
head(res)
plotMA(res)

# rlog functions to eliminate the variance dependency on the average
rld <- rlog(DDS, blind = FALSE)

re <- data.frame(res_1)
re_1 <- re %>% rownames_to_column("Geneid")

re <- data.frame(res_2)
re_2 <- re %>% rownames_to_column("Geneid")

re <- data.frame(res_3)
re_3 <- re %>% rownames_to_column("Geneid")

# Annotation
Symbol <- mapIds(org.Dr.eg.db, keys=re_1$Geneid, keytype="ENSEMBL", column="SYMBOL")
Symbol <- as.character(Symbol)
re_1$Symbol <- Symbol

# Filter genes by padj <= 0.05 and log2FoldChange > |1|
dif_exp_genes <- re_2 %>% 
  dplyr::filter((padj <= 0.05 & log2FoldChange > 1)|
                  (padj <= 0.05 & log2FoldChange < -1)) 
