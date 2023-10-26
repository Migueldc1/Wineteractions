# Author: M, de Celis Rodriguez
# Date: 28/07/2023
# Project: Wineteractions - Must metabolite and fungal diversity

library(reshape2)
library(ggplot2)
library(hillR)
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

sample_df <- subset(sample_df, Stage == "0_initial")
row.names(sample_df) <- sample_df$Seq_ID

to.rm <- c("ITS-56", "ITS-74", "ITS-92", "ITS-110")
across_wa <- c("RdG", "VLP", "LM", "M", "R1")
within_wa <- c("R1", "R2", "R3A", "R3B", "R3C")


## Community data
asv_GM <- readRDS("Data/Sequencing/Outputs/ASV_GM.rds")
asv_GM <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM <- asv_GM[, colSums(asv_GM != 0) > 0]

asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Data/Sequencing/Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]
tax_GM[is.na(tax_GM)] <- "Unidentified"


## Colors
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")
names(col_origin) <- c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C")

#
################################################################################ FIGURE 1 ####
#### PCA - GRAPE MUST ####

gm.pca_df <- sample_df[, c(13,15,16,19,20)]
row.names(gm.pca_df) <- sample_df$Sample_ID

gm.pca_df[grepl("NH4", row.names(gm.pca_df)), c("Ammonia", "PAN")] <- NA

mice.gm_pca <- mice::mice(gm.pca_df, maxit = 999, method = "pmm", seed = 1)
gm.pca_df2 <- mice::complete(mice.gm_pca, 1)


## ACROSS WINE APPELLATIONS
gm.pca_df <- gm.pca_df2[-c(57:60),]
gm_pca.ac <- prcomp(~ ., data = gm.pca_df[row.names(gm.pca_df) %in% subset(sample_df, Origin %in% across_wa)$Sample_ID, ], 
                    scale. = TRUE)

# Points
gm_pca.ac.plot <- as.data.frame(scores(gm_pca.ac))
gm_pca.ac.plot <- merge(sample_df[,c(1:5)], gm_pca.ac.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(gm_pca.ac.plot) <- gm_pca.ac.plot$Sample_ID

## PCA Plot (all points)
gg.gm.ac_pca <- ggplot(gm_pca.ac.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca.ac$sdev)^2 / sum((gm_pca.ac$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca.ac$sdev)^2 / sum((gm_pca.ac$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.gm.ac_pca

gm_dist.ac <- as.matrix(dist(scale(gm.pca_df2[row.names(gm_pca.ac.plot),])))


#
## WITHIN WINE APPELLATIONS
gm_pca.wtn <- prcomp(~ ., data = gm.pca_df2[row.names(gm.pca_df2) %in% subset(sample_df, Origin %in% within_wa)$Sample_ID, ], 
                     scale. = TRUE)

# Points
gm_pca.wtn.plot <- as.data.frame(scores(gm_pca.wtn))
gm_pca.wtn.plot <- merge(sample_df[,c(1:5)], gm_pca.wtn.plot, by.x = "Sample_ID", by.y = "row.names")

row.names(gm_pca.wtn.plot) <- gm_pca.wtn.plot$Sample_ID

## PCA Plot (all points)
gg.gm.wtn_pca <- ggplot(gm_pca.wtn.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca.wtn$sdev)^2 / sum((gm_pca.wtn$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca.wtn$sdev)^2 / sum((gm_pca.wtn$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.gm.wtn_pca

gm_dist.wtn <- as.matrix(dist(scale(gm.pca_df2[row.names(gm_pca.wtn.plot),])))


#
#### BETA DIVERSITY - GENUS LEVEL ####

asv.t_gen <- merge(tax_GM[,c(6,8)], asv.t_GM, by.x = "Id", by.y = "row.names")
asv.t_gen <- aggregate(. ~ Genus, asv.t_gen[,-1], sum)


## ACROSS WINE APPELLATIONS
bray_gen.ac <- as.matrix(vegdist(t(asv.t_gen[,colnames(asv.t_gen) %in% subset(sample_df, Origin %in% across_wa)$Seq_ID]), 
                              method = "bray", binary = TRUE))

set.seed(1)
nMDS_gen.ac <- metaMDS(bray_gen.ac)
nMDS_gen.ac$stress

nMDS_gen.ac.plot <- as.data.frame(nMDS_gen.ac[["points"]])
nMDS_gen.ac.plot$Seq_ID <- rownames(nMDS_gen.ac.plot)
nMDS_gen.ac.plot <- merge(nMDS_gen.ac.plot, sample_df, by = "Seq_ID")

gg.nmds.gen.ac_bray <- ggplot(nMDS_gen.ac.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  scale_color_manual(values = col_origin) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_gen.ac$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.nmds.gen.ac_bray


#
## WITHIN WINE APPELLATIONS
bray_gen.wtn <- as.matrix(vegdist(t(asv.t_gen[,colnames(asv.t_gen) %in% subset(sample_df, Origin %in% within_wa)$Seq_ID]), 
                                  method = "bray", binary = TRUE))

set.seed(1)
nMDS_gen.wtn <- metaMDS(bray_gen.wtn)
nMDS_gen.wtn$stress

nMDS_gen.wtn.plot <- as.data.frame(nMDS_gen.wtn[["points"]])
nMDS_gen.wtn.plot$Seq_ID <- rownames(nMDS_gen.wtn.plot)
nMDS_gen.wtn.plot <- merge(nMDS_gen.wtn.plot, sample_df, by = "Seq_ID")

gg.nmds.gen.wtn_bray <- ggplot(nMDS_gen.wtn.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Farming, size = Farming), stroke = 3, show.legend = FALSE) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  scale_color_manual(values = col_origin) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_gen.wtn$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.nmds.gen.wtn_bray


#
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
  scale_fill_manual(name = "Genus", values = c("#b15928", "#ccebc5", "#ffff33", "#e6ab02", "#80b1d3", "#bf5b17", "#e6f5c9", 
                                               "#ff7f00", "#fff2ae", "#8da0cb", "#bebada", "#1b9e77", "#fbb4ae", "#6a3d9a")) +
  guides(fill = guide_legend(nrow = 3)) + 
  facet_grid(Farming ~ Origin) +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.gen


### STATISTICs

asv.t_GM.g <- melt(asv.t_GM)
colnames(asv.t_GM.g) <- c("Id", "Seq_ID", "value")
asv.t_GM.g <- merge(asv.t_GM.g, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.g <- merge(asv.t_GM.g, sample_df[c(1:5)], by = "Seq_ID")

asv.t_GM.g <- aggregate(value ~ Seq_ID + Genus + Farming, asv.t_GM.g, sum)

wilcox.test(subset(asv.t_GM.g, Genus == "Lachancea" & Farming == "Conventional")$value,
            subset(asv.t_GM.g, Genus == "Lachancea" & Farming == "Organic")$value)


#
#### ALPHA DIVERSITY ANALYSIS ####

alpha_GM <- cbind.data.frame(# Hill based taxonomic alpha diversity
                             Richness = hill_taxa(t(asv.t_GM), q = 0),
                             `Hill–Shannon` = hill_taxa(t(asv.t_GM), q = 1),
                             `Hill–Simpson` = hill_taxa(t(asv.t_GM), q = 2))

alpha_GM$Seq_ID <- row.names(alpha_GM)
alpha_GM <- merge(alpha_GM, sample_df, by = "Seq_ID")

alpha_GM.plot <- melt(alpha_GM[,1:7])

wilcox_df <- NULL
for (var in unique(alpha_GM.plot$variable)) {
  
  alpha_sub <- subset(alpha_GM.plot, variable == var)
  wilcox_df <- rbind(wilcox_df,
                     cbind.data.frame(variable = var, 
                                      p.value = wilcox.test(alpha_sub[alpha_sub$Farming == "Conventional", "value"],
                                                            alpha_sub[alpha_sub$Farming == "Organic", "value"])$p.value))
  
}

alpha_GM.plot <- merge(alpha_GM.plot, wilcox_df, by = "variable")
alpha_GM.plot$p.value <- ifelse(alpha_GM.plot$p.value <= 0.001, "***",
                                ifelse(alpha_GM.plot$p.value <= 0.01, "**",
                                       ifelse(alpha_GM.plot$p.value <= 0.05, "*", "")))

## Management
gg.alpha_far <- ggplot(data = alpha_GM.plot, aes(x = Farming, y = value)) +
  geom_jitter(size = 2) +
  geom_boxplot(size = 1, alpha = 0.75) +
  geom_text(aes(x = 1.5, y = Inf, label = p.value), hjust = 0.5, vjust = 1.25, size = 8) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.alpha_far

#
#### EXPORT FIGURE 1 ####

## LEGEND FARMING PRACTICES
gg.legend.fp <- ggplot(nMDS_gen.ac.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, shape = Farming), stroke = 3) +
  scale_shape_manual(values = c(16, 5)) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1)))) +
  theme_bw() +
  theme(legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.box.just = "center")

gg.legend.fp <- get_legend(gg.legend.fp)


## LEGEND ACROSS WINE APPELLATIONS
gg.legend.ac <- ggplot(nMDS_gen.ac.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin), stroke = 3) +
  scale_color_manual(values = col_origin) +
  guides(color = guide_legend(title = "Across WA", override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.box.just = "center")

gg.legend.ac <- get_legend(gg.legend.ac)


## LEGEND WITHIN WINE APPELLATIONS
gg.legend.wtn <- ggplot(nMDS_gen.wtn.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin), stroke = 3) +
  scale_color_manual(values = col_origin) +
  guides(color = guide_legend(title = "Within WA", override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.box.just = "center")

gg.legend.wtn <- get_legend(gg.legend.wtn)


## CREATE THE FIGURE
gg.figure1 <- plot_grid(plot_grid(plot_grid(gg.gm.ac_pca, gg.gm.wtn_pca, ncol = 1),
                                  plot_grid(gg.legend.ac, gg.legend.fp, gg.legend.wtn, ncol = 1, rel_heights = c(1,0.25,1)),
                                  plot_grid(gg.nmds.gen.ac_bray, gg.nmds.gen.wtn_bray, ncol = 1),
                                  nrow = 1, rel_widths = c(1, 0.3, 1), labels = c("A", "", "B"), label_size = 18), 
                        plot_grid(gg.gen, gg.alpha_far, rel_widths = c(1, 0.4), labels = c("C", "D"), label_size = 18),
                        ncol = 1)
gg.figure1

ggsave("Figures/Figure_1.png", gg.figure1, bg = "white", width = 12.6, height = 11)

#
################################################################################ SUPPLEMENTARY FIGURE S2 ####
#### PCA - GRAPE MUST ####
## ACROSS WINE APPELLATIONS
gm_pca.ac2 <- prcomp(~ ., data = gm.pca_df2[row.names(gm.pca_df2) %in% subset(sample_df, Origin %in% across_wa)$Sample_ID, ], 
                     scale. = TRUE)

# Points
gm_pca.ac.plot2 <- as.data.frame(scores(gm_pca.ac2))
gm_pca.ac.plot2 <- merge(sample_df[,c(1:5)], gm_pca.ac.plot2, by.x = "Sample_ID", by.y = "row.names")

row.names(gm_pca.ac.plot2) <- gm_pca.ac.plot2$Sample_ID

## PCA Plot (all points)
gg.gm.ac_pca2 <- ggplot(gm_pca.ac.plot2) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_color_manual(values = col_origin) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  xlab(paste("PC1: ", round(((gm_pca.ac2$sdev)^2 / sum((gm_pca.ac2$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca.ac2$sdev)^2 / sum((gm_pca.ac2$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.gm.ac_pca2

#
#### EXPORT SUPPLEMENTARY FIGURE S2 ####

gg.figureS2 <- gg.gm.ac_pca2

ggsave("Figures/Figure_S2.png", gg.figureS2, bg = "white", width = 7, height = 7)

#
################################################################################ SUPPLEMENTARY FIGURE S3 ####
#### BOXPLOT - GRAPE MUST ####

gm.df_plot <- merge(sample_df[,c(1,3:5)], gm.pca_df2, by.x = "Sample_ID", by.y = "row.names")
gm.df_plot <- melt(gm.df_plot[gm.df_plot$Origin != "RdG", -5], id.vars = c("Sample_ID", "Origin", "Farming", "Condition"))
gm.df_plot$Scale <- ifelse(gm.df_plot$Origin %in% across_wa, "Across WA", "Within WA")

gm.df_plot$variable <- factor(gm.df_plot$variable, levels = c("Sugars", "Ammonia", "PAN", "Malic_acid"))

## WILCOXON TEST
wilcox_df <- NULL
for (var in levels(gm.df_plot$variable)) {
  for (sc in unique(gm.df_plot$Scale)) {
    
    metab_sub.conv <- subset(gm.df_plot, variable == var & Scale == sc & Farming == "Conventional")
    metab_sub.org <- subset(gm.df_plot, variable == var & Scale == sc & Farming == "Organic")
    
    wilcox_df <- rbind(wilcox_df,
                       cbind.data.frame(variable = var, Scale = sc,
                                        p.value = wilcox.test(metab_sub.conv$value, metab_sub.org$value, paired = TRUE)$p.value))
  }
}

wilcox_df$p.value <- ifelse(wilcox_df$p.value <= 0.001, "***",
                            ifelse(wilcox_df$p.value <= 0.01, "**",
                                   ifelse(wilcox_df$p.value <= 0.05, "*", "")))

## PLOT
gg.gm_box <- ggplot() + 
  geom_jitter(data = gm.df_plot[gm.df_plot$variable != "pH_must",], aes(x = Farming, y = value), size = 2) +
  geom_boxplot(data = gm.df_plot[gm.df_plot$variable != "pH_must",], aes(x = Farming, y = value), size = 1, alpha = 0.75) +
  facet_wrap(~ Scale + variable, scales = "free_y", nrow = 2) +
  geom_text(data = wilcox_df, aes(x = 1.5, y = Inf, label = p.value), hjust = 0.5, vjust = 1.25, size = 8) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"))

gg.gm_box

#
#### EXPORT SUPPLEMENTARY FIGURE S2 ####

gg.figureS3 <- gg.gm_box

ggsave("Figures/Figure_S3.png", gg.figureS3, bg = "white", width = 7, height = 7)

#
################################################################################ SUPPLEMENTARY TABLE S1 ####
#### PERMANOVA - GRAPE MUST ####

## ACROSS WINE APPLETATIONS
set.seed(1)
adonis2(gm_dist.ac ~ Origin*Farming, data = gm_pca.ac.plot[row.names(gm_dist.ac),])


## WITHIN WINE APPLETATIONS
set.seed(1)
adonis2(gm_dist.wtn ~ Origin*Farming, data = gm_pca.wtn.plot[row.names(gm_dist.wtn),])


#
################################################################################ SUPPLEMENTARY TABLE S2 ####
#### PERMANOVA - FUNGAL COMMUNITY ####

## ACROSS WINE APPLETATIONS
set.seed(1)
adonis2(bray_gen.ac ~ Origin*Farming, data = sample_df[row.names(bray_gen.ac),])


## WITHIN WINE APPLETATIONS
set.seed(1)
adonis2(bray_gen.wtn ~ Origin*Farming, data = sample_df[row.names(bray_gen.wtn),])


#