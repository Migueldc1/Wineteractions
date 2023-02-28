# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metabolites

library(vegan)
library(ggplot2)
library(cowplot)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

#
#### LOAD DATA ####

## SAMPLE DATA
# Must
sample_GM <- read.table("Data/Metadata/sample_GM.txt", sep = "\t", header = TRUE)
sample_GM$Condition <- factor(sample_GM$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_GM <- subset(sample_GM, Stage == "0_initial")
sample_GM <- cbind.data.frame(Sample_ID = paste(sample_GM$Origin, sample_GM$Farming, sample_GM$Condition, sep = "-"), 
                              sample_GM)
sample_GM$Origin <- factor(sample_GM$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

# Wine
sample_GM.end <- read.table("Data/Metadata/sample_GM.end.txt", sep = "\t", header = TRUE)
sample_GM.end$Condition <- factor(sample_GM.end$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_GM.end$Origin <- factor(sample_GM.end$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))


## Ecosystem Services
sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(1, rep(3, 5), 2, 1, rep(4, 6)))


## COLORS
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")

#
#### PCA - GRAPE MUST ####

gm.pca_df <- sample_GM[,c(12,13,16,19)]
row.names(gm.pca_df) <- sample_GM[,1]

gm_pca <- prcomp(~ ., data = gm.pca_df, scale. = TRUE)

# Points
gm_pca.plot <- as.data.frame(scores(gm_pca))
gm_pca.plot <- merge(sample_GM[,c(1:5)], gm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

## PCA Plot (all points)
gg.gm_pca <- ggplot() + 
  geom_point(data = gm_pca.plot, aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  labs(caption = "Grape Must") + 
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        axis.text.x = element_text(size = 15, color = "black"))

gg.gm_pca

#
#### PCA - GRAPE WINE ####

gw.pca_df <- sample_GM.end[,7:24]
row.names(gw.pca_df) <- sample_GM.end[,1]

#
## MICE for filling missing data
mice_pca <- mice::mice(gw.pca_df, maxit = 999, method = "pmm", seed = 1)
gw.pca_df2 <- mice::complete(mice_pca, 1)

## Calculate PCA with filled missing data
gw_pca <- prcomp(~ ., data = gw.pca_df2[,c(1,3:6,8:11,13:18)], scale. = TRUE)

# Points
gw_pca.plot <- as.data.frame(scores(gw_pca))
gw_pca.plot <- merge(sample_GM.end[,c(1:5)], gw_pca.plot, by.x = "Sample_ID", by.y = "row.names")

## PCA Plot (all points)
gg.gw_pca <- ggplot() + 
  geom_point(data = gw_pca.plot, aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  labs(caption = "Final Wine") + 
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        axis.text.x = element_text(size = 15, color = "black"))

gg.gw_pca

#
#### EXPORT FIGURE 2 ####

gg.figure2 <- plot_grid(gg.gm_pca, gg.gw_pca, labels = c("A", "B"), label_size = 18)
gg.figure2

ggsave("Figures/Figure_2.png", gg.figure2, bg = "white", width = 18, height = 8)

#