#####--------------- Load packages ---------------#####
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)
library(org.Dr.eg.db)

#####--------------- Scaling of Log2FoldChange in range [-1, 1] ---------------#####
re_1$scaled_log2FoldChange <- sapply(re_1$log2FoldChange, function(gene) ifelse(gene > 0, gene/max(re_1$log2FoldChange), gene/-min(re_1$log2FoldChange)))
re_2$scaled_log2FoldChange <- sapply(re_2$log2FoldChange, function(gene) ifelse(gene > 0, gene/max(re_2$log2FoldChange), gene/-min(re_2$log2FoldChange)))
re_3$scaled_log2FoldChange <- sapply(re_3$log2FoldChange, function(gene) ifelse(gene > 0, gene/max(re_3$log2FoldChange), gene/-min(re_3$log2FoldChange)))

#####--------------- Main part ---------------#####
# Target genes
genes <- c("ENSDARG00000062139", "ENSDARG00000068729",
           "ENSDARG00000104276", "ENSDARG00000093182",
           "ENSDARG00000101061", "ENSDARG00000111939",
           "ENSDARG00000069135", "ENSDARG00000068128")

re_filtered_1 <- re_1[re_1$Geneid %in% genes, c("Geneid", "log2FoldChange", "padj", "scaled_log2FoldChange")] 
re_filtered_2 <- re_2[re_2$Geneid %in% genes, c("Geneid", "log2FoldChange", "padj", "scaled_log2FoldChange")] 
re_filtered_3 <- re_3[re_3$Geneid %in% genes, c("Geneid", "log2FoldChange", "padj", "scaled_log2FoldChange")]  

# annotation
Symbol <- mapIds(org.Dr.eg.db, keys=re_filtered_1$Geneid, keytype="ENSEMBL", column="SYMBOL")
Symbol <- as.character(Symbol)
re_filtered$Symbol <- Symbol

shape_1 = c("Other", "Kinase", "Other", "Other", "Kinase", "Other", "Kinase", "Kinase")
shape_2 = c("Kinase", "Other", "Other", "Kinase", "Other", "Other", "Kinase", "Kinase") 
shape_3 = c("Kinase", "Other", "Other", "Kinase", "Other", "Other", "Kinase", "Kinase") 

re_filtered_1$shape <- as.factor(shape_1)
re_filtered_2$shape <- as.factor(shape_2)
re_filtered_3$shape <- as.factor(shape_3)

# graph building
g <- graph(edges = c("eif2a", "atf4a",
                     "eif2ak1", "eif2a",
                     "eif2ak2", "eif2a",
                     "eif2ak3", "eif2a",
                     "eif2ak4", "eif2a",
                     "atf4a", "ppp1r15a",
                     "atf4a", "ppp1r15b"), 
           directed = TRUE)

V(g)$shape <- re_filtered$shape[match(V(g)$name, re_filtered$Symbol)]
V(g)$log2FoldChange <- re_filtered$scaled_log2FoldChange[match(V(g)$name, re_filtered$Symbol)]
V(g)$padj <- re_filtered$padj[match(V(g)$name, re_filtered$Symbol)]

# creating ggraph object
set.seed(42)  
ggraph(g, layout = "fr") +  
  geom_edge_arc(aes(start_cap = label_rect(node1.name), 
                    end_cap = label_rect(node2.name)), 
                arrow = arrow(length = unit(3, 'mm')),  
                lineend = "round", 
                color = "grey40",  
                width = 0.9,  
                curvature = 0.3) +  
  geom_node_point(aes(fill = V(g)$log2FoldChange,
                      shape = V(g)$shape), 
                  color = "black", size = 12, stroke = 1.2) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, vjust = 1.5) +  
  geom_node_text(aes(label = sprintf("padj = %.2f", padj)), 
                 nudge_y = -0.1,  
                 size = 3, 
                 color = "black") +
  geom_node_text(aes(label = V(g)$name), repel = TRUE, size = 3, vjust = 1.5) + 
  scale_fill_identity(name = "Log2FoldChange") +  
  scale_shape_manual(values = c(22, 21)) + 
  scale_fill_gradient2(name = "Log2FoldChange", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_void() + 
  theme(legend.position = "right") 

# save image
ggsave("graph.pdf", plot = p, device = "pdf", width = 8, height = 6)

