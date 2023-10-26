# Author: M, de Celis Rodriguez
# Date: 03/08/2023
# Project: Wineteractions - Metabolite production by different fermentative yeasts

library(vegan)
library(ggplot2)
library(cowplot)
library(reshape2)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - Wineteractions/GitHub/Wineteractions/")

#load("Scripts/Figure2.RData")

#
#### LOAD DATA ####

## SAMPLE DATA
# Synthetic must
sample_sgm <- readRDS("Data/Meta-transcriptomics/sample_sgm.rds")
sample_sgm$Condition <- factor(sample_sgm$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_sgm$Origin <- factor(sample_sgm$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
row.names(sample_sgm) <- sample_sgm$Sample_ID

sample_sgm$Farming <- ifelse(sample_sgm$Farming == "ECO", "Organic", "Conventional")


## METABOLITE FAMILY
sgm_fam <- data.frame(Family = c(rep("Sugars", 3), rep("Alcohols", 2), rep("Acidity", 7), rep("Volatiles", 6)),
                      variable = c("Sugars", "Glucose", "Fructose", "Ethanol", "Glycerol", "pH", "Total_acidity",
                                   "Acetic_acid", "Lactic_acid", "Tartaric_acid", "Citric_acid", "Succinic_acid",
                                   "EEFA", "SCFA", "MCFA", "Fusel.alcohols", "Fusel.alcohol.acetates", "Ethyl.acetate"))


## COMMUNITY DATA
# Synthetic Grape Must (SGM) samples

#RNAseq derived data
rna_sgm <- read.table("Data/Meta-transcriptomics/bracken.report.txt", check.names = FALSE)

#ITS synthetic Grape Must (SGM)
asv_sgm <- readRDS("Data/Sequencing/Outputs/ASV_sgm.rds")
asv_sgm <- apply(asv_sgm, 1, function(x) x/sum(x))

tax_sgm <- as.data.frame(readRDS("Data/Sequencing/Outputs/tax_sgm.rds"))
tax_sgm[is.na(tax_sgm)] <- "Unidentified"
tax_sgm <- gsub("^[a-z]__", "", as.matrix(tax_sgm))


## COLORS
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")
names(col_origin) <- c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C")

col_condition <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
names(col_condition) <- c("Control", "18C", "NH4", "SO2")

col_family <- c("#bbcc06", "#81c236", "#30b3bf", "#82107a")
names(col_family) <- c("Sugars", "Alcohols", "Acidity","Volatiles")


#
################################################################################ FIGURE 2 ####
#### TAXONOMIC PROFILE SYNTHETIC MUST - RNA ####

rna_sgm.genus <- cbind.data.frame(Genus = colsplit(row.names(rna_sgm), pattern = " ", names = c("Genus", ""))[1], rna_sgm)
rna_sgm.genus <- aggregate(. ~ Genus, rna_sgm.genus, sum)

rna_sgm.genus_plot <- melt(rna_sgm.genus)
rna_sgm.genus_plot <- merge(sample_sgm[,1:4], rna_sgm.genus_plot, by.x = "Sample_ID", by.y = "variable")

rna_sgm.genus_plot[rna_sgm.genus_plot$value < 0.025, "Genus"] <- "Other"
rna_sgm.genus_plot <- aggregate(value ~ Genus + Sample_ID + Origin + Farming + Condition, rna_sgm.genus_plot, sum)

orderG_SGM <- levels(factor(rna_sgm.genus_plot$Genus))
orderG_SGM <- orderG_SGM[! orderG_SGM %in% c("Other", "Unidentified")]
orderG_SGM <- append(orderG_SGM, c("Other", "Unidentified"))

rna_sgm.genus_plot$Origin <- factor(rna_sgm.genus_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
rna_sgm.genus_plot$Farming <- factor(rna_sgm.genus_plot$Farming, levels = c("Conventional", "Organic"))
rna_sgm.genus_plot$Condition <- factor(rna_sgm.genus_plot$Condition, levels = c("Control", "18C", "NH4", "SO2"))
colnames(rna_sgm.genus_plot)[6] <- "Abundance"

gg.tax.rna_SGM <- ggplot(rna_sgm.genus_plot, 
                     aes(x = Condition, y = Abundance, fill = factor(Genus, levels = orderG_SGM))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", 
                    values = c("#caf55d", "#5df5cc", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77", "#f78e4d", "#6a3d9a", "#c05b17")) +
  facet_grid(Farming ~ Origin) +
  guides(fill = guide_legend(nrow = 2)) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.margin = margin(t = -0.25))

gg.tax.rna_SGM

#
#### METABOLITE PROFILE BOXPLOT - GENUS ####

dom.gen_sgm <- aggregate(Abundance ~ Sample_ID, rna_sgm.genus_plot, max)
dom.gen_sgm <- merge(rna_sgm.genus_plot, dom.gen_sgm, by = c("Sample_ID", "Abundance"))

metabolite_sgm <- merge(sample_sgm, dom.gen_sgm[,c(1,3)], by = "Sample_ID")
metabolite_sgm$Genus <- ifelse(metabolite_sgm$Genus %in% c("Hanseniaspora", "Lachancea", "Saccharomyces"),
                               metabolite_sgm$Genus, "Other")

metabolite.sgm_plot <- melt(metabolite_sgm[,-c(7,15)])

metabolite.sgm_plot$variable <- factor(metabolite.sgm_plot$variable,
                                       levels = c("Glucose", "Fructose", "Ethanol", "Glycerol",  "pH", 
                                                  "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                                  "Citric_acid", "Succinic_acid", "EEFA", "SCFA", "MCFA", "Fusel.alcohols", 
                                                  "Fusel.alcohol.acetates", "Ethyl.acetate"))
metabolite.sgm_plot$Genus <- factor(metabolite.sgm_plot$Genus, 
                                    levels = c("Saccharomyces", "Lachancea", "Hanseniaspora", "Other"))

## ANOVA
anova_df <- NULL
for (var in levels(metabolite.sgm_plot$variable)) {
  
  groups <- agricolae::LSD.test(aov(value ~ Genus, data = subset(metabolite.sgm_plot, variable == var)), "Genus")$groups
  groups$Genus <- row.names(groups)
  groups$nosig <- ifelse(sum(grepl("a", groups$groups)) == 4, NA, groups$groups)
  groups$groups <- ifelse(is.na(groups$nosig), NA, groups$groups)
  anova_df <- rbind(anova_df, cbind.data.frame(variable = var, groups))
  
}

anova_df$variable <- factor(anova_df$variable, levels = levels(metabolite.sgm_plot$variable))

gg.meta.sgm_sgm <- ggplot(metabolite.sgm_plot) +
  geom_boxplot(aes(x = Genus, y = value, color = Genus), size = 1, show.legend = FALSE) +
  geom_label(data = anova_df, aes(x = Genus, y = Inf, label = groups), vjust = 1,
             fill = "white", alpha = 0.5, label.size = NA, size = 6.5) +
  scale_color_manual(values = c("#1b9e77", "#8da0cb", "#cc3939", "gray70")) +
  facet_wrap(~ variable, ncol = 4, scales = "free_y") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 16, color = "black"))


## Set Strip background colors
gg.meta.sgm_sgm <- ggplot_gtable(ggplot_build(gg.meta.sgm_sgm))

strip_both <- which(grepl("strip-", gg.meta.sgm_sgm$layout$name))
fills <- alpha(c("#82107a", "#82107a", "#82107a", "#82107a", "#30b3bf", "#30b3bf", "#82107a", "#82107a", "#30b3bf", 
                 "#30b3bf", "#30b3bf", "#30b3bf", "#81c236", "#81c236", "#81c236", "#bbcc06", "#bbcc06", "#30b3bf"), 
               alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.meta.sgm_sgm$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid::grid.draw(gg.meta.sgm_sgm)

#
#### EXPORT FIGURE 2 ####

gg.figure2 <- plot_grid(gg.tax.rna_SGM, gg.meta.sgm_sgm, ncol = 1, labels = c("A", "B"), rel_heights = c(0.65,1), label_size = 18)
gg.figure2

ggsave("Figures/Figure_2.png", gg.figure2, bg = "white", width = 12.6, height = 14)

#
################################################################################ SUPPLEMENTARY FIGURE S3 ####
#### PCA SYNTHETIC MUST - CONDITION ####

sgm.pca_df <- sample_sgm[,c(9:13,16,17,6,18:23)]

## Calculate PCA with filled missing data
sgm_pca <- prcomp(~ ., data = sgm.pca_df[complete.cases(sgm.pca_df),], scale. = TRUE)

# Points
sgm_pca.plot <- as.data.frame(scores(sgm_pca))
sgm_pca.plot <- merge(sample_sgm[,1:5], sgm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(sgm_pca.plot) <- sgm_pca.plot$Sample_ID

## PCA Plot (all points)
gg.sgm_pca <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Condition), shape = 21, size = 4, stroke = 1) +
  scale_fill_manual(values = col_condition) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21), order = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        axis.text.x = element_text(size = 15, color = "black"))

gg.sgm_pca

sgm_dist <- as.matrix(dist(scale(sgm.pca_df[complete.cases(sgm.pca_df),])))

set.seed(1)
adonis2(sgm_dist ~ Condition, data = sgm_pca.plot[row.names(sgm_dist),])


#
#### PCA SYNTHETIC MUST - GENUS ####

sgm_pca.plot$Genus <- factor(sgm_pca.plot$Genus, levels = c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"))

## PCA Plot (all points)
gg.sgm_pca.genus <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Genus), shape = 21, size = 4, stroke = 1) +
  scale_fill_manual(values = c("#cc3939", "#8da0cb", "#1b9e77", "gray70")) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21), order = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"))

gg.sgm_pca.genus

set.seed(1)
adonis2(sgm_dist ~ Genus, data = sgm_pca.plot[row.names(sgm_dist),])


#
#### EXPORT FIGURE S3 ####

gg.figureS3 <- plot_grid(gg.sgm_pca, gg.sgm_pca.genus, ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figureS3

ggsave("Figures/Figure_S3.png", gg.figureS3, bg = "white", width = 10, height = 12)

#
################################################################################ SUPPLEMENTARY FIGURE S4 ####
#### METABOLITE PROFILE BOXPLOT - CONDITION ####

metabolite_plot <- melt(sample_sgm[,-c(7,15)])

metabolite_plot$variable <- factor(metabolite_plot$variable, 
                                   levels = c("Glucose", "Fructose", "Ethanol", "Glycerol",  "pH", 
                                              "Acetic_acid", "Lactic_acid", "Tartaric_acid", "Citric_acid", "Succinic_acid", 
                                              "EEFA", "SCFA", "MCFA", "Fusel.alcohols", "Fusel.alcohol.acetates", 
                                              "Ethyl.acetate"))

## ANOVA
anova_df <- NULL
for (var in levels(metabolite_plot$variable)) {
  
  groups <- agricolae::LSD.test(aov(value ~ Condition, data = subset(metabolite_plot, variable == var)), "Condition")$groups
  groups$Genus <- row.names(groups)
  groups$nosig <- ifelse(sum(grepl("a", groups$groups)) == 4, NA, groups$groups)
  groups$groups <- ifelse(is.na(groups$nosig), NA, groups$groups)
  anova_df <- rbind(anova_df, cbind.data.frame(variable = var, groups))
  
}

anova_df$variable <- factor(anova_df$variable, levels = levels(metabolite_plot$variable))

gg.meta_sgm <- ggplot(metabolite_plot) +
  geom_boxplot(aes(x = Condition, y = value, color = Condition), size = 1, show.legend = FALSE) +
  geom_label(data = anova_df, aes(x = Genus, y = Inf, label = groups), vjust = 1,
             fill = "white", alpha = 0.5, label.size = NA, size = 6.5) +
  scale_color_manual(values = col_condition) +
  facet_wrap(~ variable, nrow = 6, scales = "free_y") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 16, color = "black"))


## Set Strip background colors
gg.meta_sgm <- ggplot_gtable(ggplot_build(gg.meta_sgm))

strip_both <- which(grepl("strip-", gg.meta_sgm$layout$name))
fills <- alpha(c("#82107a", NA, NA, "#82107a", "#82107a", "#82107a", "#30b3bf", "#82107a", "#82107a", "#30b3bf", "#30b3bf", "#30b3bf", 
                 "#bbcc06", "#30b3bf", "#30b3bf", "#81c236", "#81c236", "#bbcc06"), 
               alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.meta_sgm$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid::grid.draw(gg.meta_sgm)

#
#### EXPORT FIGURE S4 ####
gg.figureS4 <- gg.meta_sgm
grid::grid.draw(gg.figureS4)

ggsave("Figures/Figure_S4.png", gg.figureS3, bg = "white", width = 12.6, height = 14)

#
################################################################################ SUPPLEMENTARY FIGURE S5 ####
#### TAXONOMIC PROFILE SYNTHETIC MUST - RNA vs ITS ####

## ITS

asv_sgm.plot <- melt(asv_sgm)
colnames(asv_sgm.plot) <- c("Id", "Sample_ID", "Abundance")
asv_sgm.plot <- merge(asv_sgm.plot, tax_sgm[,6, drop = FALSE], by.x = "Id", by.y = "row.names")

asv_sgm.plot_gen <- aggregate(Abundance ~ Genus + Sample_ID, asv_sgm.plot, sum)
asv_sgm.plot_gen[asv_sgm.plot_gen$Abundance < 0.025, "Genus"] <- "Other" 
asv_sgm.plot_gen <- aggregate(Abundance ~ Genus + Sample_ID, asv_sgm.plot_gen, sum)

asv_sgm.plot_gen <- merge(sample_sgm[,1:4], asv_sgm.plot_gen, by = "Sample_ID")

## COMBINE WITH RNA

tax_sgm.comp <- rbind(cbind.data.frame(rna_sgm.genus_plot[,c(2:5,1,6)], Type = "RNA reads"),
                      cbind.data.frame(asv_sgm.plot_gen, Type = "ITS reads"))

orderG_sgm <- levels(factor(tax_sgm.comp$Genus))
orderG_sgm <- orderG_sgm[! orderG_sgm %in% c("Other", "Unidentified")]
orderG_sgm <- append(orderG_sgm, c("Other", "Unidentified"))

tax_sgm.comp$Farming <- ifelse(tax_sgm.comp$Farming == "Conventional", "CONV", "ORG")

gg.genus_sgm <- ggplot(tax_sgm.comp, 
                       aes(x = Condition, y = Abundance, fill = factor(Genus, levels = orderG_sgm))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", 
                    values = c("#caf55d", "#ed3460", "#5df5cc", "#cc3939", "#cab2d6", "#8da0cb", "#d5eb26", "#1b9e77", 
                               "#f78e4d", "#e6d8bd", "#6a3d9a", "#c05b17")) +
  facet_grid(Type ~ Origin + Farming) +
  guides(fill = guide_legend(nrow = 2)) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.genus_sgm

#
#### EXPORT FIGURE S5 ####

gg.figureS5 <- gg.genus_sgm

ggsave("Figures/Figure_S5.png", gg.figureS3, bg = "white", width = 12.6, height = 6)

#