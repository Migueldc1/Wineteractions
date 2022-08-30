# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - RNAseq analysis

# Set the project location as working directory
setwd("C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(DESeq2)
library(ggplot2)
library(reshape2)

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df <- cbind.data.frame(Sample_ID = sample_df[,1],
                              colsplit(sample_df[,1], "-",
                                       names = c("Origin", "Farming", "Condition")),
                              sample_df[,-1])
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

## Bracken Taxonomy
asv_RNA <- readRDS("Outputs/ASV_t1-RNA.rds")
asv_RNA <- apply(asv_RNA, 1, function(x) x/sum(x))

tax_RNA <- readRDS("Outputs/tax_t1-RNA.rds")
tax_RNA <- cbind.data.frame(tax_RNA, Id = row.names(tax_RNA))

asv_RNA <- melt(asv_RNA)
colnames(asv_RNA) <- c("Id", "Sample_ID", "value")
asv_RNA <- merge(asv_RNA, tax_RNA[,c(1,3)])

## Dominant Genus
dom_sp <- aggregate(asv_RNA$value, list(asv_RNA$Sample_ID), max)
colnames(dom_sp) <- c("Sample_ID", "value")
dom_sp <- merge(dom_sp, asv_RNA[,2:4], by = c("Sample_ID", "value"))

sample_df <- merge(sample_df, dom_sp[,c(1,3)], by = "Sample_ID")

## COUNT DATA
countData <- read.csv("Outputs/count_df.csv", header = TRUE, sep = ",")
colnames(countData) <- gsub("\\.", "-", colnames(countData))
countData[is.na(countData)] <- 0

## LOAD DESeq2 OBJECT
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = sample_df, 
                              design = ~Genus, tidy = TRUE)

#
#### GLOBAL PCA ####
vsdata <- vst(dds, blind = FALSE)

plotPCA(vsdata, intgroup = "Genus")









