# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(reshape2)
library(DESeq2)
library(clusterProfiler)
library(goseq)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - Wineteractions/GitHub/Wineteractions/")

#
#### CUSTOM FUNCTION ####

## Load Count data + Emapper Annotations to get KEGG Ortholog counts
read.emapper <- function(path.count, path.orth) {
  
  count.list <- NULL
  
  for (smpl in list.files(path.count)) {
    sample <- gsub(".txt", "", smpl)
    
    count_df <- read.table(paste(path.count, smpl, sep = "/"), header = TRUE)[,c(1,7)]
    colnames(count_df) <- c("seed_ortholog", sample)
    
    orth_file <- list.files(paste(path.orth, sample, sep = "/"), full.names = TRUE)
    orth_df <- read.delim(orth_file, sep = "\t", skip = 4)
    orth_df <- unique(orth_df[1:(nrow(orth_df)-3), c(2,12)])
    
    count.deff_df <- merge(count_df, orth_df, by = "seed_ortholog")
    count.deff_df <- subset(count.deff_df, KEGG_ko != "-")
    
    count.deff_df <- aggregate(count.deff_df[,2], list(count.deff_df$KEGG_ko), sum)
    colnames(count.deff_df) <- c("KEGG_ko", sample)
    
    count.list[[sample]] <- count.deff_df
    
  }
  
  ko_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), count.list)
  ko_df[is.na(ko_df)] <- 0
  
  return(ko_df)
  
}

## Load Gene Ontology Data
read.GOdata <- function(file.ext, path.annot, cols) {
  
  kegg.go_df <- NULL
  
  for (smpl in list.files(path.annot)) {
    
    #Open Annotation files and select KO and GO columns
    annot_df <- read.delim(list.files(paste(path.annot, smpl, sep = "/"), full.names = TRUE), sep = "\t", skip = 4)[,cols]
    annot_df <- unique(annot_df[annot_df$KEGG_ko != "" & annot_df$KEGG_ko != "-", ])
    
    #Bind all the annotation files in a single table, removing duplicated info
    kegg.go_df <- rbind.data.frame(kegg.go_df, annot_df)
    kegg.go_df <- unique(kegg.go_df)
    
  }
  
  #Format the table into a useful form
  go_df <- NULL
  for (n in 1:nrow(annot_df)) {
    
    go_df <- rbind.data.frame(go_df, cbind(gene_id = annot_df$KEGG_ko[n],
                                           category = unlist(strsplit(as.character(annot_df$GOs[n]),",", fixed = TRUE))))
    
    go_df <- unique(go_df)
    
  }
  
  return(go_df)
  
}

#
#### LOAD DATA ####

## KO TABLE
ko_df <- read.emapper(path.count = "Data/Meta-transcriptomics/Count/eggnog", 
                      path.orth = "Data/Meta-transcriptomics/Annotation/eggnog")
ko_df.f <- ko_df[!grepl(",", ko_df$KEGG_ko),]


## ANNOTATION TABLES
kegg_df <- read.table("Data/Meta-transcriptomics/Annotation/KEGG_names.txt", sep = "\t", header = TRUE, quote = "")
go_df <- read.GOdata(".emapper.annotations", "Data/Meta-transcriptomics/Annotation/eggnog", c("KEGG_ko", "GOs"))


## SAMPLE DATA
sample_sgm <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_sgm$Condition <- factor(sample_sgm$Condition, levels = c("Control", "18C", "NH4", "SO2"))


## BRACKEN REPORT
#RNAseq derived data
rna_sgm <- read.table("Data/Meta-transcriptomics/bracken.report.txt", check.names = FALSE)


#
#### GET DOMINANT GENUS ####
rna_sgm.genus <- cbind.data.frame(Genus = colsplit(row.names(rna_sgm), pattern = " ", names = c("Genus", ""))[1], rna_sgm)
rna_sgm.genus <- aggregate(. ~ Genus, rna_sgm.genus, sum)
rna_sgm.genus <- melt(rna_sgm.genus)
rna_sgm.genus <- merge(sample_sgm[,1:4], rna_sgm.genus, by.x = "Sample_ID", by.y = "variable")

dom.gen_sgm <- aggregate(value ~ Sample_ID, rna_sgm.genus, max)
dom.gen_sgm <- merge(rna_sgm.genus, dom.gen_sgm, by = c("Sample_ID", "value"))
dom.gen_sgm$Genus <- ifelse(dom.gen_sgm$value >= 0.9, dom.gen_sgm$Genus, "Other")

sample_sgm <- merge(dom.gen_sgm[,c(1,3:6)], sample_sgm, by = c("Sample_ID", "Origin", "Farming", "Condition"))
sample_sgm$Genus <- ifelse(sample_sgm$Genus %in% c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"), sample_sgm$Genus, "Other")

#
#### DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS - GLOBAL ####

dds_sgm <- DESeqDataSetFromMatrix(countData = ko_df.f, 
                                  colData = sample_sgm, 
                                  design = ~ Condition + Genus + Origin, tidy = TRUE)

dds_sgm$Genus <- relevel(dds_sgm$Genus, "Saccharomyces")

dds_sgm <- DESeq(dds_sgm, fitType = "local")
resultsNames(dds_sgm)

#
#### DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS - YEASTS ####

# KEEP SAMPLES DOMINATED BY ONE OF SAID SPECIES
sample_gen <- sample_sgm[sample_sgm$Genus %in% c("Saccharomyces", "Lachancea", "Hanseniaspora"),]
ko_df.gen <- ko_df.f[,c("KEGG_ko", sample_gen$Sample_ID)]

dds_gen <- DESeqDataSetFromMatrix(countData = ko_df.gen, 
                                  colData = sample_gen, 
                                  design = ~ Condition + Origin + Genus, tidy = TRUE)

dds_gen$Genus <- relevel(dds_gen$Genus, "Saccharomyces")
keep <- rowSums(counts(dds_gen) >= 10) >= 5
dds_gen <- dds_gen[keep,]

dds_gen <- DESeq(dds_gen, fitType = "local", betaPrior = FALSE)
resultsNames(dds_gen)

## GET DIFFERENTIALLY EXPRESSED ORTHOLOGS
# Saccharomyces v Lachancea
res_lt.sc <- results(dds_gen, contrast = c("Genus", "Lachancea", "Saccharomyces"), alpha = 0.05)
res_lt.sc <- lfcShrink(dds_gen, res = res_lt.sc, type = "ashr")
summary(res_lt.sc)

ress_lt.sc <- as.data.frame(res_lt.sc)
ress_lt.sc$KEGG_ko <- row.names(ress_lt.sc)

ress_lt.sc <- merge(ress_lt.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.sc <- ress_lt.sc[order(ress_lt.sc$padj),]
ress_lt.sc$DEO <- ifelse(abs(ress_lt.sc$log2FoldChange) > 1 & ress_lt.sc$padj <= 0.05, 1, 0)

# Saccharomyces v Hanseniaspora
res_hs.sc <- results(dds_gen, contrast = c("Genus", "Hanseniaspora", "Saccharomyces"), alpha = 0.05)
res_hs.sc <- lfcShrink(dds_gen, res = res_hs.sc, type = "ashr")
summary(res_hs.sc)

ress_hs.sc <- as.data.frame(res_hs.sc)
ress_hs.sc$KEGG_ko <- row.names(ress_hs.sc)

ress_hs.sc <- merge(ress_hs.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.sc <- ress_hs.sc[order(ress_hs.sc$padj),]
ress_hs.sc$DEO <- ifelse(abs(ress_hs.sc$log2FoldChange) > 1 & ress_hs.sc$padj <= 0.05, 1, 0)

#
#### BIOLOGICAL ENRICHMENT - GENE ONTOLOGY (BP) - YEASTS ####
count.bias <- rowSums(ko_df.gen[,-1])
names(count.bias) <- ko_df.gen[,1]

DEG_lt.sc <- as.integer(ress_lt.sc$padj < 0.05 & abs(ress_lt.sc$log2FoldChange) >= 1)
names(DEG_lt.sc) <- ress_lt.sc$KEGG_ko

pwf_lt.sc <- nullp(DEgenes = DEG_lt.sc, bias.data = count.bias[names(DEG_lt.sc)])
go_lt.sc <- goseq(pwf_lt.sc, gene2cat = go_df, test.cats = c("GO:BP"))
go_lt.sc <- subset(go_lt.sc, ontology == "BP" & numDEInCat > 5)

enrichGO_lt.sc <- go_lt.sc[p.adjust(go_lt.sc$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_hs.sc <- as.integer(ress_hs.sc$padj < 0.05 & abs(ress_hs.sc$log2FoldChange) >= 1)
names(DEG_hs.sc) <- ress_hs.sc$KEGG_ko

pwf_hs.sc <- nullp(DEgenes = DEG_hs.sc, bias.data = jitter(rep(1000, length(DEG_hs.sc))))
go_hs.sc <- goseq(pwf_hs.sc, gene2cat = go_df)
go_hs.sc <- subset(go_hs.sc, ontology == "BP" & numDEInCat > 5)

enrichGO_hs.sc <- go_hs.sc[p.adjust(go_hs.sc$over_represented_pvalue, method = "fdr") < 0.05,]

#
#### DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS - SACCHAROMYCES ####
sample_sc <- sample_sgm[sample_sgm$Genus == "Saccharomyces",]
sample_sc <- sample_sc[-c(5,6),]

ko_sc <- ko_df.f[,c("KEGG_ko", sample_sc$Sample_ID)]

dds_sc <- DESeqDataSetFromMatrix(countData = ko_sc, 
                                 colData = sample_sc, 
                                 design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc)) >= 10
dds_sc <- dds_sc[keep,]

dds_sc <- DESeq(dds_sc, fitType = "local")
resultsNames(dds_sc)

## DIFFERENTIAL EXPRESSION

# Low Temperature
res_sc.18C <- results(dds_sc, contrast = c("Condition", "18C", "Control"), alpha = 0.05)
res_sc.18C <- lfcShrink(dds_sc, res = res_sc.18C, contrast = c("Condition", "18C", "Control"), type = "normal")
summary(res_sc.18C)

ress_sc.18C <- as.data.frame(res_sc.18C)
ress_sc.18C$KEGG_ko <- row.names(ress_sc.18C)

ress_sc.18C <- merge(ress_sc.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.18C <- ress_sc.18C[order(ress_sc.18C$padj),]
ress_sc.18C$DEO <- ifelse(abs(ress_sc.18C$log2FoldChange) > 1 & ress_sc.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_sc.NH4 <- results(dds_sc, contrast = c("Condition", "NH4", "Control"), alpha = 0.05)
res_sc.NH4 <- lfcShrink(dds_sc, res = res_sc.NH4, contrast = c("Condition", "NH4", "Control"), type = "normal")
summary(res_sc.NH4)

ress_sc.NH4 <- as.data.frame(res_sc.NH4)
ress_sc.NH4$KEGG_ko <- row.names(ress_sc.NH4)

ress_sc.NH4 <- merge(ress_sc.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.NH4 <- ress_sc.NH4[order(ress_sc.NH4$padj),]
ress_sc.NH4$DEO <- ifelse(abs(ress_sc.NH4$log2FoldChange) > 1 & ress_sc.NH4$padj <= 0.05, 1, 0)

# High Sulfites
res_sc.SO2 <- results(dds_sc, contrast = c("Condition", "SO2", "Control"), alpha = 0.05)
res_sc.SO2 <- lfcShrink(dds_sc, res = res_sc.SO2, contrast = c("Condition", "SO2", "Control"), type = "normal")
summary(res_sc.SO2)

ress_sc.SO2 <- as.data.frame(res_sc.SO2)
ress_sc.SO2$KEGG_ko <- row.names(ress_sc.SO2)

ress_sc.SO2 <- merge(ress_sc.SO2, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.SO2 <- ress_sc.SO2[order(ress_sc.SO2$padj),]
ress_sc.SO2$DEO <- ifelse(abs(ress_sc.SO2$log2FoldChange) > 1 & ress_sc.SO2$padj <= 0.05, 1, 0)

#
#### BIOLOGICAL ENRICHMENT - GENE ONTOLOGY (BP) - SACCHAROMYCES ####
count.bias <- rowSums(ko_sc[,-1])
names(count.bias) <- ko_sc[,1]

DEG_sc.18C <- as.integer(ress_sc.18C$padj < 0.05 & abs(ress_sc.18C$log2FoldChange) >= 1)
names(DEG_sc.18C) <- ress_sc.18C$KEGG_ko

pwf_sc.18C <- nullp(DEgenes = DEG_sc.18C, bias.data = count.bias[names(DEG_sc.18C)])
go_sc.18C <- goseq(pwf_sc.18C, gene2cat = go_df, test.cats = c("GO:BP"))
go_sc.18C <- subset(go_sc.18C, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.18C <- go_sc.18C[p.adjust(go_sc.18C$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_sc.NH4 <- as.integer(ress_sc.NH4$padj < 0.05 & abs(ress_sc.NH4$log2FoldChange) >= 1)
names(DEG_sc.NH4) <- ress_sc.NH4$KEGG_ko

pwf_sc.NH4 <- nullp(DEgenes = DEG_sc.NH4, bias.data = count.bias[names(DEG_sc.NH4)])
go_sc.NH4 <- goseq(pwf_sc.NH4, gene2cat = go_df)
go_sc.NH4 <- subset(go_sc.NH4, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.NH4 <- go_sc.NH4[p.adjust(go_sc.NH4$over_represented_pvalue, method = "fdr") < 0.05,]

#
#### EXPORT TABLES ####
## SAMPLE DATA
saveRDS(sample_sgm, "Data/Meta-transcriptomics/sample_sgm.rds")


## ORTHOLOGS DATA
write.table(go_df, "Data/Meta-transcriptomics/Annotation/go_df.txt", sep = "\t", row.names = FALSE)
saveRDS(counts(dds_sgm, normalized = TRUE), "Data/Meta-transcriptomics/ko.n_df.rds")


## DESeq2 object (All samples)
saveRDS(dds_sgm, "Data/Meta-transcriptomics/dds_sgm.rds")
saveRDS(dds_gen, "Data/Meta-transcriptomics/dds_gen.rds")
saveRDS(dds_sc, "Data/Meta-transcriptomics/dds_sc.rds")


## DEO TABLES
# YEASTS
write.table(ress_lt.sc, "Data/Meta-transcriptomics/ress_lt.sc.txt", sep = "\t", row.names = FALSE)
write.table(ress_hs.sc, "Data/Meta-transcriptomics/ress_hs.sc.txt", sep = "\t", row.names = FALSE)

# SACCHAROMYCES
write.table(ress_sc.18C, "Data/Meta-transcriptomics/ress_sc.18C.txt", sep = "\t", row.names = FALSE)
write.table(ress_sc.NH4, "Data/Meta-transcriptomics/ress_sc.NH4.txt", sep = "\t", row.names = FALSE)
write.table(ress_sc.SO2, "Data/Meta-transcriptomics/ress_sc.SO2.txt", sep = "\t", row.names = FALSE)


## BIOLOGICAL ENRICHMENT TABLES
# YEASTS
write.table(enrichGO_lt.sc, "Data/Meta-transcriptomics/enrichGO_lt.sc.txt", sep = "\t", row.names = FALSE)
write.table(enrichGO_hs.sc, "Data/Meta-transcriptomics/enrichGO_hs.sc.txt", sep = "\t", row.names = FALSE)

# SACCHAROMYCES
write.table(enrichGO_sc.18C, "Data/Meta-transcriptomics/enrichGO_sc.18C.txt", sep = "\t", row.names = FALSE)
write.table(enrichGO_sc.NH4, "Data/Meta-transcriptomics/enrichGO_sc.NH4.txt", sep = "\t", row.names = FALSE)

#