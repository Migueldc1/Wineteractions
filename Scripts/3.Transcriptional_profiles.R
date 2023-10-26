# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Transcriptional profiles of fermenting yeast communities

library(reshape2)
library(DESeq2)
library(ggplot2)
library(ggforce)
library(cowplot)
library(vegan)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - Wineteractions/GitHub/Wineteractions/")

#
#### LOAD DATA ####

## SAMPLE DATA
sample_sgm <- readRDS("Data/Meta-transcriptomics/sample_sgm.rds")
row.names(sample_sgm) <- sample_sgm$Sample_ID


## DESeq2 OBJECTS
dds_sgm <- readRDS("Data/Meta-transcriptomics/dds_sgm.rds")


## DEO TABLES
ress_lt.sc <- read.table("Data/Meta-transcriptomics/ress_lt.sc.txt", header = TRUE, sep = "\t")
ress_hs.sc <- read.table("Data/Meta-transcriptomics/ress_hs.sc.txt", header = TRUE, sep = "\t")


## BIOLOGICAL ENRICHMENT
enrichGO_lt.sc <- read.table("Data/Meta-transcriptomics/enrichGO_lt.sc.txt", header = TRUE, sep = "\t")
enrichGO_hs.sc <- read.table("Data/Meta-transcriptomics/enrichGO_hs.sc.txt", header = TRUE, sep = "\t")


## COLORS
col_cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")

#
################################################################################ FIGURE 3 ####
#### PCA GLOBAL - GENUS ####

vst_sgm <- vst(dds_sgm, blind = TRUE)

pcaData_sgm <- plotPCA(vst_sgm, intgroup = "Genus", returnData = TRUE, ntop = 5000)
pcaData_sgm$Genus <- factor(pcaData_sgm$Genus, levels = c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"))
percentVar_sgm <- round(100 * attr(pcaData_sgm, "percentVar"), 2)

pcaData_sgm <- merge(pcaData_sgm, sample_sgm[,1:4], by.x = "name", by.y = "Sample_ID")

gg.pca_gen <- ggplot(pcaData_sgm, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#cc3939", "#8da0cb", "#1b9e77", "gray70")) + 
  xlab(paste0("PC1: ", percentVar_sgm[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sgm[2],"% variance")) + 
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.margin = margin(t = -0.25, r = 0, b = 0, l = 0, unit = "cm"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 14, color = "black"))

gg.pca_gen

#
#### PERMANOVAs ####

norm_smg.dis <- as.matrix(dist(t(assay(vst_sgm)), method = "euclidean"))

adonis2(norm_smg.dis ~ Genus*Condition, sample_sgm[row.names(norm_smg.dis),])

#
#### DIFFERENTIAL EXPRESSION - VENN DIAGRAM ####

venn.df_gen <- merge(ress_lt.sc[,c(1, 9)], ress_hs.sc[,c(1, 9)], by = "KEGG_ko")
colnames(venn.df_gen) <- c("KEGG_ko", "Lt.Sc", "Hs.Sc")

venn.plot_gen <- rbind.data.frame(cbind(Lt.Sc = 0, Hs.Sc = 0, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 0)),
                                  cbind(Lt.Sc = 0, Hs.Sc = 1, 
                                        Counts = sum(venn.df_gen[,2] == 0 & venn.df_gen[,3] == 1)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 0, 
                                        Counts = sum(venn.df_gen[,2] == 1 & venn.df_gen[,3] == 0)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 1, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 2)))

venn.plot_gen <- cbind.data.frame(venn.plot_gen, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out2 <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("Lachancea", "Hanseniaspora"))
venn.out2$labels <- factor(venn.out2$labels, levels = c("Lachancea", "Hanseniaspora"))

gg.venn_gen <- ggplot() +
  geom_circle(data = venn.out2, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_gen, aes(x = x, y = y, label = Counts), size = 6) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 13, color  = "black"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5)) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = NULL)

gg.venn_gen

#
#### DIFFERENTIAL EXPRESSION - HISTOGRAM ####

hist_cond <- rbind(cbind(subset(ress_lt.sc, DEO == 1), Genus = "Lachancea"),
                   cbind(subset(ress_hs.sc, DEO == 1), Genus = "Hanseniaspora"))

gg.hist_gen <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Genus)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30", show.legend = FALSE) +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DEO") +
  scale_fill_manual(values = c("#cc3939", "#8da0cb")) +
  theme(axis.text.y = element_text(size = 14, color  = "black"),
        axis.title.x = element_text(size = 15, color  = "black"),
        axis.title.y = element_text(size = 15, color  = "black"),
        axis.text.x = element_text(size = 14, color  = "black")) 

gg.hist_gen

#
#### BIOLOGICAL ENRICHMENT vs SACCHAROMYCES ####

go_gen <- rbind.data.frame(cbind(enrichGO_lt.sc, Comparison = "Lachancea"),
                           cbind(enrichGO_hs.sc, Comparison = "Hanseniaspora"))

go_gen$Comparison <- factor(go_gen$Comparison, levels = c("Lachancea", "Hanseniaspora"))
go_gen$comp.cat <- ifelse(duplicated(go_gen$term) | duplicated(go_gen$term, fromLast = TRUE), "Both", go_gen$Comparison)

go_gen$term <- factor(go_gen$term, 
                      levels = unique(go_gen$term[order(go_gen$comp.cat, go_gen$numDEInCat, decreasing = TRUE)]))

gg.go_gen <- ggplot(go_gen) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_gen

#
#### EXPORT FIGURE 3 ####

gg.legend3 <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Genus)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = c("#cc3939", "#8da0cb")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15, color  = "black"),
        legend.title = element_text(size = 15, color  = "black"),
        legend.margin = margin(-0.25,-0.25,-0.25,-0.25)) 

gg.legend3 <- get_legend(gg.legend3)

gg.figure3 <- plot_grid(gg.pca_gen,
                        plot_grid(gg.venn_gen, gg.hist_gen, rel_widths = c(0.75, 1), labels = c("B", "C"), label_size = 18),
                        gg.legend3, gg.go_gen, rel_heights = c(1.2, 0.8, 0.10, 2), ncol = 1, labels = c("A", "", "", "D"), label_size = 18)
gg.figure3

ggsave("Figures/Figure_3.png", gg.figure3, bg = "white", width = 8, height = 15)

#
################################################################################ FIGURE 3 ####
#### PCA GLOBAL - CONDITION ####

gg.pca_con <- ggplot(pcaData_sgm, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_cond) + 
  xlab(paste0("PC1: ", percentVar_sgm[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sgm[2],"% variance")) + 
  guides(color = guide_legend(nrow = 2)) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.margin = margin(t = -0.25, r = 0, b = 0, l = 0, unit = "cm"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 14, color = "black"))

gg.pca_con

#
#### EXPORT SUPPLEMENTARY FIGURE S6 ####

gg.figureS6 <- gg.pca_con
gg.figureS6

ggsave("Figures/Figure_S6.png", gg.figureS6, bg = "white", width = 8, height = 8)

#