# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(vegan)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(igraph)
library(grid)
#library(cowplot)
#library(ggforce)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

## Load Environment
#load("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/RNAseq-meta.RData")
#save.image("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/RNAseq-meta.RData")

#
#### CUSTOM FUNCTION ####

#
#### LOAD DATA ####

## KEGG ORTHOLOGY TABLE
ko_df <- readRDS("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/ko_df.rds")


## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Genus <- factor(sample_df$Genus, levels = c("Citeromyces", "Hanseniaspora", "Kluyveromyces", "Lachancea", 
                                                      "Saccharomyces", "Other"))

sample_df$Specie <- gsub("Hanseniaspora ", "H. ", sample_df$Specie)


## COLORS
col.cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col.gen <- c("#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77", "gray70")


## Ecosystem Services
sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(6, rep(8, 5), 7, 6, rep(9, 6)))
#
#### PCA ####

sgm.pca_df <- sample_df[,c(5:16,19:24)]
row.names(sgm.pca_df) <- sample_df[,1]

#
## MICE for filling missing data
mice_pca <- mice::mice(sgm.pca_df, maxit = 999, method = "pmm", seed = 1)
sgm.pca_df2 <- mice::complete(mice_pca, 1)

## Calculate PCA with filled missing data
sgm_pca <- prcomp(~ ., data = sgm.pca_df2[,c(1,4:8,11:18)], scale. = TRUE)

# Points
sgm_pca.plot <- as.data.frame(scores(sgm_pca))
sgm_pca.plot <- merge(sample_df[,c(1:4,17)], sgm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

# Vectors
sgm_pca.vplot <- scores(sgm_pca, display = "species")
sgm_pca.vplot <- cbind.data.frame(variable = row.names(sgm_pca.vplot), sgm_pca.vplot)
sgm_pca.vplot <- merge(sgm_group, sgm_pca.vplot[,1:3], by = "variable")
sgm_pca.vplot$Group <- factor(sgm_pca.vplot$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))

## PCA Plot (all points)
gg.sgm_pca <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Genus), shape = 21, size = 3) +
  scale_fill_manual(values = col.gen) +
  geom_segment(data = sgm_pca.vplot, aes(x = 0, y = 0, xend = PC1*11, yend = PC2*11, color = Group), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) + 
  geom_text(data = sgm_pca.vplot, aes(x = PC1*10.5, y = PC2*11.5, label = variable, color = Group), 
            size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("#bbcc06", "#81c236", "#30b3bf", "#82107a")) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21), order = 1), 
         shape = guide_legend(override.aes = list(fill = "black"), order = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        axis.text.x = element_text(size = 15, color = "black"))

gg.sgm_pca

ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/sgm_pca.png", gg.sgm_pca, 
       width = 11, height = 10)


#
#### METABOLITE - ORTHOLOG NETWORK ####

ko_df.t <- as.data.frame(t(ko_df[,-1]))
colnames(ko_df.t) <- ko_df[,1]

ko.sgm_df <- merge(sgm.pca_df2[,c(1,4:8,11:18)], ko_df.t, by = "row.names")
row.names(ko.sgm_df) <- ko.sgm_df[,1]
ko.sgm_df <- ko.sgm_df[,-1]

ko.sgm_corr <- rcorr(as.matrix(ko.sgm_df), type = "spearman")$r
ko.sgm_corr[upper.tri(ko.sgm_corr, diag = TRUE)] <- NA
ko.sgm_corp <- rcorr(as.matrix(ko.sgm_df), type = "spearman")$P
ko.sgm_corp[upper.tri(ko.sgm_corp, diag = TRUE)] <- NA

ko.sgm_df <- cbind.data.frame(melt(ko.sgm_corr[15:2912,1:14]),
                              melt(ko.sgm_corp[15:2912,1:14]))[,c(1:3,6)]

colnames(ko.sgm_df) <- c("Source", "Target", "value", "p.value")
ko.sgm_df <- ko.sgm_df[complete.cases(ko.sgm_df),]
ko.sgm_df$p.adjust <- p.adjust(ko.sgm_df$p.value, method = "fdr")

ko.sgm_df <- subset(ko.sgm_df, p.adjust <= 0.05)

##Draw Network
ko.sgm_net <- graph_from_data_frame(as.matrix(ko.sgm_df[,1:2]), directed = FALSE)

set.seed(1)
l <- layout_with_fr(ko.sgm_net, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

ko.sgm_cwt <- cluster_walktrap(ko.sgm_net)
max(ko.sgm_cwt$membership)

V(ko.sgm_net)$community <- ko.sgm_cwt$membership
V(ko.sgm_net)$community <- ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, sgm_group$cols, V(ko.sgm_net)$community)

colrs <- adjustcolor(c("darkorchid", "deepskyblue", "sienna1", "yellowgreen", "#a52a2a", 
                                   "yellow", "#b2ed6d", "#375bb0", "#d128c5"), alpha = 0.8)

png(file = "Workflow/4_Metabolomic_analysis/Outputs/Figures/ko.sgm_net.png", width = 1920, height = 1080)

plot(ko.sgm_net, vertex.size = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, 8, 5), 
     vertex.label = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, names(V(ko.sgm_net)), NA), 
     vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
     vertex.shape = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, "square", "circle"),
     vertex.color = colrs[V(ko.sgm_net)$community], 
     rescale = FALSE, layout = l*1.3, edge.color = "gray70", 
     edge.curved = 0.5)

dev.off()

library("ggplotify")

p1 <- as.grob(~plot(ko.sgm_net, vertex.size = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, 8, 5), 
                    vertex.label = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, names(V(ko.sgm_net)), NA), 
                    vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
                    vertex.shape = ifelse(names(V(ko.sgm_net)) %in% sgm_group$variable, "square", "circle"),
                    vertex.color = colrs[V(ko.sgm_net)$community], 
                    rescale = FALSE, layout = l*1.3, edge.color = "gray70", 
                    edge.curved = 0.5))


gg.ko.sgm_net <- ggplotify::as.ggplot(p.ko.sgm_net)

ggsave("Workflow/soil.bs.rr_corr.png", p1, bg = "white",
       width = 11, height = 7)




