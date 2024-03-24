#####--------------- Load packages ---------------#####
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)


# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#####--------------- Enrichmant analysis (GSEA) ---------------#####
organism = "org.Dr.eg.db"

dif_exp_genes <- re_1 %>% filter(padj < 0.05)
dif_exp_genes <- re_2 %>% filter(padj < 0.05)
dif_exp_genes <- re_3 %>% filter(padj < 0.05)

new_dif <- dif_exp_genes[order(dif_exp_genes$log2FoldChange, decreasing = TRUE),]

original_gene_list <- new_dif$log2FoldChange

names(original_gene_list) <- new_dif$Symbol
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Gene Ontology  
gseBP_1<- gseGO(geneList=gene_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.01, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")


gseMF_1 <- gseGO(geneList=gene_list, 
               ont ="MF", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.01, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")

gseCC_1 <- gseGO(geneList=gene_list, 
               ont ="CC", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.01, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")

# Visulization
dotplot(gseBP_1, 
        color = "pvalue", 
        showCategory=20, 
        split=".sign", 
        font.size = 14, 
        label_format = 60) + facet_grid(.~.sign)

dotplot(gseMF_1, 
        color = "pvalue", 
        showCategory=20, 
        split=".sign", 
        font.size = 10, 
        label_format = 65) + facet_grid(.~.sign)

dotplot(gseCC_1, 
        color = "pvalue", 
        showCategory=20, 
        split=".sign", 
        font.size = 10, 
        label_format = 65) + facet_grid(.~.sign)

cnetplot(gseBP_1, categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = categories_TI_BP, 
         cex_label_category = 0.8, 
         cex_label_gene = 0.5)

cnetplot(gseMF_1, categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 10, 
         cex_label_category = 0.8, 
         cex_label_gene = 0.5)

cnetplot(gseCC_1, categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 10, 
         cex_label_category = 0.8, 
         cex_label_gene = 0.5)

# KEGG
# Transfer the list of genes to ENTREZ ID
ENTREZ <- mapIds(org.Dr.eg.db, keys=new_dif$Geneid, keytype="ENSEMBL", column="ENTREZID")
new_dif$ENTREZ <- ENTREZ

# we want the log2 fold change 
kegg_gene_list <- new_dif$log2FoldChange

# name the vector
names(kegg_gene_list) <- new_dif$ENTREZ
kegg_df <- na.omit(data.frame(as.numeric(new_dif$ENTREZ), as.numeric(kegg_gene_list)))

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "dre"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

# Visulization
dotplot(kk2, 
        color = "pvalue", 
        showCategory = 20, 
        title = "Enriched Pathways", 
        split=".sign", font.size = 14, 
        label_format = 60) + facet_grid(.~.sign)

cnetplot(kk2, 
         categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 8, 
         color_gene = "red")

# Produce the native KEGG plot (PNG) by pathview 
kegg_path <- function(pathway_list){
  for (i in pathway_list)
    pathview(gene.data=kegg_gene_list, pathway.id = i, species = kegg_organism)
}

# exclude the wrong pathways 
pathway_list <- kegg_res$ID[kegg_res$ID != "dre01200" 
                            & kegg_res$ID != "dre01232"
                            & kegg_res$ID != "dre01250"
                            & kegg_res$ID != "dre00511"
                            & kegg_res$ID != "dre00563"
                            & kegg_res$ID != "dre00531"
                            & kegg_res$ID != "dre01230"
                            & kegg_res$ID != "dre01210"
                            & kegg_res$ID != "dre01240"
                            & kegg_res$ID != "dre00510"
                            & kegg_res$ID != "dre00534"] 

kegg_path(pathway_list)

