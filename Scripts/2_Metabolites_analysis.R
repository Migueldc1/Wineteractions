# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metabolites

library(vegan)
library(ggplot2)
library(cowplot)
library(reshape2)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

#load("Scripts/Figure2.RData")

#
#### LOAD DATA ####

## SAMPLE DATA
# Must
sample_GM <- read.table("Data/Metadata/sample_GM.txt", sep = "\t", header = TRUE)
sample_GM$Condition <- factor(sample_GM$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_GM <- cbind.data.frame(Sample_ID = paste(sample_GM$Origin, sample_GM$Farming, sample_GM$Condition, sep = "-"), 
                              sample_GM)
sample_GM$Origin <- factor(sample_GM$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
sample_GM$Farming <- ifelse(sample_GM$Farming == "ECO", "Organic", "Conventional")


# Wine
sample_GM.end <- read.table("Data/Metadata/sample_GM.end.txt", sep = "\t", header = TRUE)
sample_GM.end$Condition <- factor(sample_GM.end$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_GM.end$Origin <- factor(sample_GM.end$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
sample_GM.end$Farming <- ifelse(sample_GM.end$Farming == "ECO", "Organic", "Conventional")


# Synthetic must
sample_SGM <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_SGM$Condition <- factor(sample_SGM$Condition, levels = c("Control", "18C", "NH4", "SO2"))
row.names(sample_SGM) <- sample_SGM$Sample_ID

sample_SGM$Farming <- ifelse(sample_SGM$Farming == "ECO", "Organic", "Conventional")


## Taxonomic Information

# Grape Must
asv_GM <- readRDS("Data/Sequencing/Outputs/ASV_GM.rds")
asv_GM <- apply(asv_GM, 1, function(x) x/sum(x))
asv_GM <- asv_GM[, colnames(asv_GM) %in% subset(sample_GM, Stage == "1_middle")$Seq_ID]

tax_GM <- as.data.frame(readRDS("Data/Sequencing/Outputs/tax_GM.rds"))
tax_GM[is.na(tax_GM)] <- "Unidentified"


# Synthetic Must
tax_SGM <- read.table("Data/Meta-transcriptomics/bracken.report.txt", check.names = FALSE)


## Ecosystem Services
sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(1, rep(3, 5), 2, 1, rep(4, 6)))


## COLORS
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")

#
################################### FIGURE 2 ####
#### PCA - GRAPE MUST ####

sample_GM0 <- subset(sample_GM, Stage == "0_initial")

gm.pca_df <- sample_GM0[,c(13,15,16,19,20)]
row.names(gm.pca_df) <- sample_GM0[,1]

gm.pca_df[grepl("NH4", row.names(gm.pca_df)), c("Ammonia", "PAN")] <- NA

mice.gm_pca <- mice::mice(gm.pca_df, maxit = 999, method = "pmm", seed = 1)
gm.pca_df2 <- mice::complete(mice.gm_pca, 1)

gm_pca <- prcomp(~ ., data = gm.pca_df2, scale. = TRUE)

# Points
gm_pca.plot <- as.data.frame(scores(gm_pca))
gm_pca.plot <- merge(sample_GM0[,c(1:5)], gm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(gm_pca.plot) <- gm_pca.plot$Sample_ID

## PCA Plot (all points)
gg.gm_pca <- ggplot(gm_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  labs(caption = "Grape Must") + 
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        plot.margin = margin(t = 0.2, r = 0.2, b = -0.2, l = 0.2, unit = "cm"))

gg.gm_pca

gm_dist <- as.matrix(dist(scale(gm.pca_df2[-c(57:60),])))

set.seed(1)
adonis2(gm_dist ~ Origin*Farming, data = gm_pca.plot[row.names(gm_dist),])

#
#### PCA - GRAPE WINE ####

gw.pca_df <- sample_GM.end[,7:24]
row.names(gw.pca_df) <- sample_GM.end[,1]

#
## MICE for filling missing data
#mice.gw_pca <- mice::mice(gw.pca_df, maxit = 999, method = "pmm", seed = 1)
#gw.pca_df2 <- mice::complete(mice.gw_pca, 1)

## Calculate PCA with filled missing data
gw_pca <- prcomp(~ ., data = gw.pca_df[complete.cases(gw.pca_df),c(1,3:6,8:11,13:18)], scale. = TRUE)

# Points
gw_pca.plot <- as.data.frame(scores(gw_pca))
gw_pca.plot <- merge(sample_GM.end[,1:5], gw_pca.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(gw_pca.plot) <- gw_pca.plot$Sample_ID

## PCA Plot (all points)
gg.gw_pca <- ggplot(gw_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  labs(caption = "Fermented Must") + 
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        plot.margin = margin(t = 0.2, r = 0.2, b = -0.2, l = 0.2, unit = "cm"))

gg.gw_pca

gw_dist <- as.matrix(dist(scale(gw.pca_df[complete.cases(gw.pca_df),c(1,3:6,8:11,13:18)])))

set.seed(1)
adonis2(gw_dist ~ Origin*Farming, data = gw_pca.plot[row.names(gw_dist),])


#
#### VARIANCE PARTITIONING - GRAPE WINE ####
## Need to run first the Supplementary Figures Script

gw_varp <- varpart(gw_dist,
                   ~ Origin+Farming,
                   ~ Condition,
                   ~ Genus,
                   data = gw_pca.plot[row.names(gw_dist),])

plot(gw_varp, digits = 4, Xnames = c("Vineyard", "Condition", "Genus"),
     bg = c("#e3b352", "#9152e3", "#26962f"))

gw.varp_df <- gw_varp$part$indfract
gw.varp_df$variable <- c("Origin", "Condition", "Genus", "None", "None", "Origin|Genus", "All", "Residuals")
gw.varp_df$variable <- factor(gw.varp_df$variable, 
                              levels = c("Origin", "Origin|Genus", "Genus", "Condition", "All", "Residuals"))

gw.varp_df <- subset(gw.varp_df, Adj.R.square > 0)
gw.varp_df$Adj.R.square <- gw.varp_df$Adj.R.square / sum(gw.varp_df$Adj.R.square)

gw.varp <- ggplot(gw.varp_df) +
  geom_bar(aes(x = 1, y = Adj.R.square, fill = variable), stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("#f0d467", "#acd164", "#58bf6e", "#bf5037", "#b00707", "#9e9e9d")) +
  guides(fill = guide_legend(nrow = 1)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_line(linewidth = 1, color = "gray80"),
        panel.grid.minor.x = element_line(linewidth = 0.75, color = "gray90"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(t = 0.5, b = 0.2, l = 0.25, r = 0.25, unit = "cm"))

gw.varp

#
#### EXPORT FIGURE 2 ####
gg.dummy_pca <- ggplot(gw_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3), nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        legend.margin = margin(b = -0.25, unit = "cm"),
        plot.margin = margin(b = 0.5, unit = "cm"))

gg.dummy_pca

gg.pca_legend <- get_legend(gg.dummy_pca)

gg.figure2 <- plot_grid(plot_grid(gg.gm_pca, gg.gw_pca, labels = c("A", "B"), label_size = 18),
                        gg.pca_legend, ncol = 1, rel_heights = c(1.25, 0.2))
gg.figure2

gg.figure2_new <- plot_grid(gg.figure2, gw.varp, ncol = 1, rel_heights = c(1, 0.4), labels = c("", "C"), label_size = 18, label_y = 1.05)
gg.figure2_new

ggsave("Figures/Figure_2.png", gg.figure2_new, bg = "white", width = 12.6, height = 9)

#
################################### RELATED SUPPLEMENTARY INFORMATION ####
#### TAXONOMIC PROFILE - GRAPE MUST ####

tax_GM.mid <- aggregate(asv_GM, list(tax_GM$Genus), sum)
colnames(tax_GM.mid)[1] <- "Genus"
tax_GM.mid$Genus <- gsub("g__", "", tax_GM.mid$Genus)
tax_GM.mid <- melt(tax_GM.mid)

tax_GM.mid[tax_GM.mid$value < 0.025, "Genus"] <- "Other"
tax_GM.mid <- aggregate(value ~ Genus + variable, tax_GM.mid, sum)

orderG_GM <- levels(factor(tax_GM.mid$Genus))
orderG_GM <- orderG_GM[! orderG_GM %in% c("Other", "Unidentified")]
orderG_GM <- append(orderG_GM, c("Other", "Unidentified"))

tax.GM_plot <- cbind.data.frame(tax_GM.mid, colsplit(tax_GM.mid$variable, "-", c("Origin", "Farming", "Condition")))
tax.GM_plot <- merge(tax_GM.mid, sample_GM[,2:5], by.x = "variable", by.y = "Seq_ID")

tax.GM_plot$Origin <- factor(tax.GM_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
tax.GM_plot$Condition <- factor(tax.GM_plot$Condition, levels = c("Control", "18C", "NH4", "SO2"))

gg.tax_GM <- ggplot(tax.GM_plot, 
                 aes(x = Condition, y = value, fill = factor(Genus, levels = orderG_GM))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ffff33", "#80b1d3", "#ed3460", "#7fcc2d", "#cab2d6", "#8da0cb",
                                               "#1b9e77", "#e6d8bd", "#edc31c", "#6a3d9a", "#bf5b17")) +
  facet_grid(Farming ~ Origin) +
  guides(fill = guide_legend(nrow = 2)) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.margin = margin(t = -0.25))

gg.tax_GM

#
#### PCA - GRAPE WINE - GENUS ####

gw_pca.genus <- aggregate(value ~ variable, tax.GM_plot, max)
gw_pca.genus <- merge(gw_pca.genus, tax.GM_plot[,1:3], by = c("variable", "value"))
gw_pca.genus <- merge(gw_pca.genus, sample_GM[,1:2], by.x = "variable", by.y = "Seq_ID")

gw_pca.plot <- merge(gw_pca.plot, gw_pca.genus[,2:4], by = "Sample_ID", all.x = TRUE)
gw_pca.plot$Genus[is.na(gw_pca.plot$Genus)] <- "Unidentified"

row.names(gw_pca.plot) <- gw_pca.plot$Sample_ID

gg.gw_pca2 <- ggplot(gw_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Genus, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = c("#ffff33", "#ed3460", "#cab2d6", "#8da0cb", "#1b9e77", "#bf5b17")) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gw_pca$sdev)^2 / sum((gw_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  labs(caption = "Final Wine") + 
  theme_bw() +
  theme(aspect.ratio = 0.75,
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"))

gg.gw_pca2

set.seed(1)
adonis2(gw_dist ~ Origin*Farming + Condition + Genus, data = gw_pca.plot[row.names(gw_dist),])

#
#### EXPORT SUPPLEMENTARY FIGURE 2 ####

gg.figureS2 <- plot_grid(gg.tax_GM, gg.gw_pca2, ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figureS2

ggsave("Figures/Figure_S2.png", gg.figureS2, bg = "white", width = 11, height = 11)

#
################################### FIGURE 3 ####
#### TAXONOMIC PROFILE - SYNTHETIC MUST ####

tax_SGM.genus <- cbind.data.frame(Genus = colsplit(row.names(tax_SGM), pattern = " ", names = c("Genus", "to.rm"))[1], tax_SGM)
tax_SGM.g <- melt(tax_SGM.genus)

tax_SGM.g[tax_SGM.g$value < 0.025, "Genus"] <- "Other"
tax_SGM.g <- aggregate(value ~ Genus + variable, tax_SGM.g, sum)

orderG_SGM <- levels(factor(tax_SGM.g$Genus))
orderG_SGM <- orderG_SGM[! orderG_SGM %in% c("Other", "Unidentified")]
orderG_SGM <- append(orderG_SGM, c("Other", "Unidentified"))

tax.SGM_plot <- cbind.data.frame(tax_SGM.g, colsplit(tax_SGM.g$variable, "-", c("Origin", "Farming", "Condition")))
tax.SGM_plot$Farming <- ifelse(tax.SGM_plot$Farming == "ECO", "Organic", "Conventional")

tax.SGM_plot$Origin <- factor(tax.SGM_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
tax.SGM_plot$Condition <- factor(tax.SGM_plot$Condition, levels = c("Control", "18C", "NH4", "SO2"))

gg.tax_SGM <- ggplot(tax.SGM_plot, 
                     aes(x = Condition, y = value, fill = factor(Genus, levels = orderG_SGM))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", 
                    values = c("#caf55d", "#5df5cc", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77", "#f78e4d", "#6a3d9a", "#c05b17")) +
  facet_grid(Farming ~ Origin) +
  guides(fill = guide_legend(nrow = 2)) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.margin = margin(t = -0.25))

gg.tax_SGM

#
#### PCA - SYNTHETIC MUST ####

sgm.pca_df <- sample_SGM[,c(8:12,15,16,5,17:22)]

## Calculate PCA with filled missing data
sgm_pca <- prcomp(~ ., data = sgm.pca_df[complete.cases(sgm.pca_df),], scale. = TRUE)

# Points
sgm_pca.plot <- as.data.frame(scores(sgm_pca))
sgm_pca.plot <- merge(sample_SGM[,1:4], sgm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(sgm_pca.plot) <- sgm_pca.plot$Sample_ID

# Vectors
sgm_pca.vplot <- scores(sgm_pca, display = "species")
sgm_pca.vplot <- cbind.data.frame(variable = row.names(sgm_pca.vplot), sgm_pca.vplot)
sgm_pca.vplot <- merge(sgm_group, sgm_pca.vplot[,1:3], by = "variable")
sgm_pca.vplot$Group <- factor(sgm_pca.vplot$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))

sgm_pca.genus <- aggregate(value ~ variable, tax.SGM_plot, max)
sgm_pca.genus <- merge(sgm_pca.genus, tax.SGM_plot[,1:3], by = c("variable", "value"))

sgm_pca.plot <- merge(sgm_pca.plot, sgm_pca.genus[,1:3], by.x = "Sample_ID", by.y = "variable", all.x = TRUE)
row.names(sgm_pca.plot) <- sgm_pca.plot$Sample_ID

## PCA Plot (all points)
gg.sgm_pca <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Genus), shape = 21, size = 3) +
  scale_fill_manual(values = c("#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77")) +
  geom_segment(data = sgm_pca.vplot, aes(x = 0, y = 0, xend = PC1*11, yend = PC2*11, color = Group), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) + 
  geom_text(data = sgm_pca.vplot, aes(x = PC1*10.5, y = PC2*11.5, label = variable, color = Group), 
            size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("#bbcc06", "#81c236", "#30b3bf", "#82107a")) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(shape = 21), order = 1), 
         shape = guide_legend(override.aes = list(fill = "black"), order = 1)) +
  xlim(-4.5,5.5) +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        aspect.ratio = 0.75,
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black", hjust = 0),
        legend.title = element_text(size = 16, color = "black", hjust = 0),
        axis.text.x = element_text(size = 15, color = "black"))

gg.sgm_pca

sgm_dist <- as.matrix(dist(scale(sgm.pca_df[complete.cases(sgm.pca_df),])))

set.seed(1)
adonis2(sgm_dist ~ Origin*(Genus+Farming), data = sgm_pca.plot[row.names(sgm_dist),])

#
#### VARIANCE PARTITIONING - SYNTHETIC MUST ####
## Need to run first the Supplementary Figures Script

sgm_varp <- varpart(sgm_dist,
                   ~ Origin+Farming,
                   ~ Condition,
                   ~ Genus,
                   data = sgm_pca.plot[row.names(sgm_dist),])

plot(sgm_varp, digits = 4, Xnames = c("Vineyard", "Condition", "Genus"),
     bg = c("#e3b352", "#9152e3", "#26962f"))

sgm.varp_df <- sgm_varp$part$indfract
sgm.varp_df$variable <- c("Origin", "Condition", "Genus", "None", "Genus|Condition", "Origin|Genus", "None", "Residuals")
sgm.varp_df$variable <- factor(sgm.varp_df$variable, 
                              levels = c("Origin", "Origin|Genus", "Genus", "Genus|Condition", "Condition","Residuals"))

sgm.varp_df <- subset(sgm.varp_df, Adj.R.square > 0)
sgm.varp_df$Adj.R.square <- sgm.varp_df$Adj.R.square / sum(sgm.varp_df$Adj.R.square)

sgm.varp <- ggplot(sgm.varp_df) +
  geom_bar(aes(x = 1, y = Adj.R.square, fill = variable), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#f0d467", "#acd164", "#58bf6e", "#8c8853", "#bf5037", "#9e9e9d")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, color = "black"),
        panel.grid.major.y = element_line(linewidth = 1, color = "gray80"),
        panel.grid.minor.y = element_line(linewidth = 0.75, color = "gray90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(t = 0.5, b = 0.2, l = 0.25, r = 0.25, unit = "cm"))

sgm.varp

#
#### EXPORT FIGURE 3 ####

gg.figure3 <- plot_grid(gg.tax_SGM, plot_grid(sgm.varp, gg.sgm_pca, nrow = 1, rel_widths = c(0.5, 1), labels = c("B", "C"), label_size = 18), 
                        ncol = 1, rel_heights = c(0.8, 1), labels = "A", label_size = 18)
gg.figure3

ggsave("Figures/Figure_3.png", gg.figure3, bg = "white", width = 12.6, height = 10)

#