# Author: M, de Celis Rodriguez
# Date: 14/11/2023
# Project: Wineteractions - Fermented grape must metabolite and fungal diversity

library(reshape2)
library(ggplot2)
library(vegan)
library(cowplot)


rm(list=ls()) #Clear R environment

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - Wineteractions/GitHub/Wineteractions/")

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_GM.txt", sep = "\t", header = TRUE, encoding = "UTF-8")
sample_df <- cbind.data.frame(Sample_ID = paste(sample_df$Origin, sample_df$Farming, sample_df$Condition, sep = "-"), sample_df)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha", "Madrid", "Rioja"))
sample_df$Farming <- ifelse(sample_df$Farming == "ECO", "Organic", "Conventional")

sample_df <- subset(sample_df, Stage == "2_final")
sample_df <- sample_df[!is.na(sample_df$Seq_ID),]
row.names(sample_df) <- sample_df$Seq_ID

## METABOLITE DATA
must_df <- read.table("Data/Metadata/sample_GM.end.txt", sep = "\t", header = TRUE, encoding = "UTF-8")
row.names(must_df) <- must_df$Sample_ID

## COMMUNITY DATA
asv_GM <- readRDS("Data/Sequencing/Outputs/ASV_GM.rds")
asv_GM <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM <- asv_GM[, colSums(asv_GM != 0) > 0]

asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Data/Sequencing/Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]
tax_GM[is.na(tax_GM)] <- "Unidentified"


## COLORS
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")
names(col_origin) <- c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C")

#
################################################################################ SUPPLEMENTARY FIGURE S4 ####
#### TAXONOMIC EXPLORATION ####
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1:5)], by = "Seq_ID")

## Genus
asv.t_plot <- aggregate(value ~ Sample_ID + Origin + Farming + Condition + Genus, asv.t_GM.p, sum)
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(value ~ Sample_ID + Origin + Farming + Condition + Genus, asv.t_plot, sum)

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

gg.gen <- ggplot(asv.t_plot, 
                 aes(x = Condition, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ffff33", "#80b1d3", "#e6f5c9", "#ff7f00", "#cab2d6", "#8da0cb", "#1b9e77", "#6a3d9a")) +
  guides(fill = guide_legend(nrow = 2)) + 
  facet_grid(Farming ~ Origin) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.gen

#
#### PCA - FERMENTED GRAPE MUST ####

gm.pca_df <- must_df[, c("Sugars", "Ethanol", "Glycerol",  "pH", 
                         "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                         "Citric_acid", "Succinic_acid", "EEFA", "SCFA", "MCFA", "Fusel.alcohols", 
                         "Fusel.alcohol.acetates", "Ethyl.acetate")]
row.names(gm.pca_df) <- must_df$Sample_ID

gm_pca <- prcomp(~ ., data = gm.pca_df, scale. = TRUE)

# Points
gm_pca.plot <- as.data.frame(scores(gm_pca))
gm_pca.plot <- merge(must_df[,c(1:4)], gm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(gm_pca.plot) <- gm_pca.plot$Sample_ID

## PCA Plot (all points)
gg.gm_pca <- ggplot(gm_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.gm_pca

#

gm_dist <- as.matrix(dist(scale(gm.pca_df)))
gm_dist <- gm_dist[-34,-34]

set.seed(1)
adonis2(gm_dist ~ Origin*Farming + Condition, data = must_df[row.names(gm_dist),])

#
#### EXPORT FIGURE S4 ####

gg.figureS4 <- plot_grid(gg.gen, gg.gm_pca, ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figureS4

ggsave("Figures/Figure_S4.png", gg.figureS4, bg = "white", width = 10, height = 10)

#