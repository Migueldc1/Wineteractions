# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggforce)
library(clusterProfiler)
library(reshape2)
library(goseq)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

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

ko_df <- read.emapper(path.count = "Data/Meta-transcriptomics/Count/eggnog", 
                      path.orth = "Data/Meta-transcriptomics/Annotation/eggnog")

ko_df.f <- ko_df[!grepl(",", ko_df$KEGG_ko),]

## ANNOTATION
kegg_df <- read.table("Data/Meta-transcriptomics/Annotation/KEGG_names.txt", 
                      sep = "\t", header = TRUE, quote = "")

go_df <- read.GOdata(".emapper.annotations", "Data/Meta-transcriptomics/Annotation/eggnog", c("KEGG_ko", "GOs"))

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

tax_max <- read.table("Data/Meta-transcriptomics/bracken.report.txt", check.names = FALSE)
tax_max <- cbind.data.frame(Genus = colsplit(row.names(tax_max), pattern = " ", names = c("Genus", "to.rm"))[1], tax_max)
tax_max <- aggregate(. ~ Genus, tax_max, sum)
tax_max <- cbind.data.frame(Sample_ID = colnames(tax_max[,-1]), Genus = tax_max$Genus[apply(tax_max[,-1], 2, which.max)],
                            value = apply(tax_max[,-1], 2, max))
tax_max$Genus <- ifelse(tax_max$value >= 0.9, tax_max$Genus, "Other")
tax_max$Genus <- ifelse(tax_max$Genus %in% c("Saccharomyces", "Hanseniaspora", "Lachancea"), tax_max$Genus, "Other")

sample_df <- merge(sample_df, tax_max[,-3], by = "Sample_ID")


## COLORS
col_cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")

#
#### DIFFERENTIAL ANALYSIS - SACCHAROMYCES ####
sample_sc <- sample_df[sample_df$Genus == "Saccharomyces",]
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

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_sc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                     list(ress_sc.18C[,c(1,10)], ress_sc.NH4[,c(1,10)], ress_sc.SO2[,c(1,10)]))
colnames(venn.df_sc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sc[is.na(venn.df_sc)] <- 0

venn.plot_sc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 0)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,2] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,3] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,4] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,4] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,3] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,2] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 3)))

venn.plot_sc <- cbind.data.frame(venn.plot_sc, 
                                 x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                 y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col_cond[-1]) +
  labs(fill = NULL)

gg.venn_sc

# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                   cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col_cond[-1]) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 14, color  = "black"),
        axis.title.x = element_text(size = 16, color  = "black"),
        axis.title.y = element_text(size = 16, color  = "black"),
        legend.text = element_text(size = 15, color  = "black"),
        legend.title = element_text(size = 16, color  = "black"),
        axis.text.x = element_text(size = 14, color  = "black")) 

gg.hist_sc




#
#### EXPORT FIGURE 5 ####

gg.figure5 <- plot_grid(plot_grid(gg.venn_sc, gg.hist_sc, labels = c("A", "B"), label_size = 18), 
                        gg.go_sc, rel_heights = c(1, 1.5), labels = c(NA, "C"), label_size = 18, ncol = 1)
gg.figure5

ggsave("Figures/Figure_5.png", gg.figure5, bg = "white", width = 12.6, height = 11)

#