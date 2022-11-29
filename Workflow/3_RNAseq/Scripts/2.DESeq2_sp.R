# Author: M, de Celis Rodriguez
# Date: 13/09/2022
# Project: Wineteractions - RNAseq analysis

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(reshape2)
library(DESeq2)
library(ggplot2)
library(limma)
library(ggforce)
library(clusterProfiler)
library(goseq)

#load("DESeq_sp.Rdata")

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

## COUNT DATA
count.list <- list()
i <- 0

for (count in list.files("Inputs/RNA_data/Count/UniprotKB/")) {
  if (grepl(".summary", count, fixed = TRUE) == FALSE) {
    
    i <- i+1
    
    table <- read.table(paste("Inputs/RNA_data/Count/UniprotKB/", count, sep = "/"))
    colnames(table) <- table[1,]
    table <- table[-1,]
    colnames(table)[ncol(table)] <- gsub(".txt", "", count)
    count.list[[i]] <- table[,c(1,6,7)]
    names(count.list)[i] <- gsub(".txt", "", count)
    
  }
  
}

count_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "Geneid"), count.list)

## TAX DATA 
tax_sacc <- read.table("Inputs/tax_sacc.txt", header = TRUE)
tax_lach <- read.table("Inputs/tax_lach.txt", header = TRUE)

## COUNT DATA
countData <- cbind.data.frame(Geneid = count_df[,1],
                              mapply(count_df[,seq(3,119,2)], FUN = as.numeric))

countData_sacc <- countData[,colnames(countData) %in% subset(sample_df, Genus == "Saccharomyces")$Sample_ID]
countData_sacc <- cbind.data.frame(Geneid = countData[,1], countData_sacc)
countData_sacc <- countData_sacc[countData_sacc$Geneid %in% tax_sacc$Geneid,]
countData_sacc[is.na(countData_sacc)] <- 0

countData_lach <- countData[,colnames(countData) %in% subset(sample_df, Genus == "Lachancea")$Sample_ID]
countData_lach <- cbind.data.frame(Geneid = countData[,1], countData_lach)
countData_lach <- countData_lach[countData_lach$Geneid %in% tax_lach$Geneid,]
countData_lach[is.na(countData_lach)] <- 0

countData_hans <- countData[,colnames(countData) %in% subset(sample_df, Genus == "Hanseniaspora")$Sample_ID]
countData_hans <- cbind.data.frame(Geneid = countData[,1], countData_hans)
countData_hans <- countData_hans[rowSums(is.na(countData_hans[,-1])) < ncol(countData_hans[,-1])*0.5,]
countData_hans[is.na(countData_hans)] <- 0

## GENE LENGTH
length_df <- cbind.data.frame(Geneid = count_df[,1],
                              length = rowMeans(mapply(count_df[,seq(2,119,2)], FUN = as.numeric), na.rm = TRUE))

## LOAD DESeq2 OBJECT
dds_sacc <- DESeqDataSetFromMatrix(countData = countData_sacc, 
                                   colData = sample_df[sample_df$Sample_ID %in% colnames(countData_sacc),], 
                                   design = ~Condition, tidy = TRUE)

dds_lach <- DESeqDataSetFromMatrix(countData = countData_lach, 
                                   colData = sample_df[sample_df$Sample_ID %in% colnames(countData_lach),], 
                                   design = ~Condition, tidy = TRUE)

dds_hans <- DESeqDataSetFromMatrix(countData = countData_hans, 
                                   colData = sample_df[sample_df$Sample_ID %in% colnames(countData_hans),], 
                                   design = ~Condition, tidy = TRUE)

#
#### DESeq2 ####
dds_sacc <- DESeq(dds_sacc)

dds_lach <- DESeq(dds_lach)

dds_hans <- DESeq(dds_hans)

#
#### GLOBAL PCA ####
vst <- vst(dds_sacc, blind = TRUE)

pcaData <- plotPCA(vst, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
ggplot(pcaData, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#1b9e77", "#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "gray")) +
  labs(color = "Dominant Genus")

#


#### SAVE DATA ####
save.image("DESeq_sp.Rdata")
