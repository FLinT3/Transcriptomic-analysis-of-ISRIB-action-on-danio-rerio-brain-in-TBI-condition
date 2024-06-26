#####--------------- Load packages ---------------#####
library(GenomicFeatures)
library(tximportData)
library(tximport)
library(readxl)
library(xfun)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(ggbeeswarm)
library(tximeta)
library(AnnotationDbi)

#####--------------- Data preparation ---------------#####
# MetaData loading
MetaDATA <- read_xlsx("Metada.xlsx")
MetaDATA <- MetaDATA[order(MetaDATA$group, decreasing = FALSE),]
MetaDATA$condition <- NULL

# remove C3 sample
MetaDATA <- subset(MetaDATA, MetaDATA$sample != "C3")

# add factors
MetaDATA %>% mutate(T = ifelse(Condition == "T" | Condition == "TI", 1, 0),
                    I = ifelse(Condition == "I" | Condition == "TI", 1, 0)) -> MetaDATA

SRR_list <- as.integer(MetaDATA$group)
sample_list <- MetaDATA$sample

# Annotation
# txdb <- makeTxDbFromGFF("Danio_rerio.GRCz11.110.gtf")
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- read.table("D:/Ubunta files/RNAseq/ZEBRA/Annotation.txt", sep = "\t", header = T)

# verification
length(unique(tx2gene$GENEID))

# Create RefSeq Transcript List
files <- file.path("counts", SRR_list, "quant.sf")
names(files) <- sample_list
all(file.exists(files))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi.salmon)

# Total number of reads per library
barplot(colSums(txi.salmon$counts)/1e+6,las=2,col='steelblue', ylab="Expression, mln reads")

# sample table creation
sampleTable <- data.frame(condition_T = as.factor(MetaDATA$T), condition_I = as.factor(MetaDATA$I),
                          groups = as.factor(MetaDATA$Condition))
 
