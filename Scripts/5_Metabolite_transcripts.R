# Author: M, de Celis Rodriguez
# Date: 04/07/2023
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(cowplot)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

#
#### LOAD DATA ####

## ORTHOLOG TABLES
#Normalized count table
ko.n_df <- readRDS("Data/Meta-transcriptomics/ko.n_df.rds")

#DEO genus
#DEO.gen_df <- readRDS("Data/Meta-transcriptomics/DEO.gen_df.rds")


## SYNTHETIC MUST
sample_df <- readRDS("Data/Metadata/sample_SGM.rds")
row.names(sample_df) <- sample_df$Sample_ID
sample_df$Farming <- ifelse(sample_df$Farming == "ECO", "Organic", "Conventional")

sgm.pca_df <- sample_df[,c(8:12,15,16,5,17:22)]


## METABOLITE GROUP
sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(1, rep(3, 5), 2, 1, rep(4, 6)))

#
#### METABOLITE - ORTHOLOG NETWORK - GLOBAL ####

ko.t_df <- as.data.frame(t(ko.n_df))
ko.sgm_df <- merge(sgm.pca_df, ko.t_df, by = "row.names")
row.names(ko.sgm_df) <- ko.sgm_df[,1]
ko.sgm_df <- ko.sgm_df[,-1]

## SPEARMAN'S CORRELATIONS
net.cor_df <- rcorr(as.matrix(ko.sgm_df), type = "spearman")

net.cor.r_df <- net.cor_df$r
net.cor.r_df[upper.tri(net.cor.r_df, diag = TRUE)] <- NA
net.cor.p_df <- net.cor_df$P
net.cor.p_df[upper.tri(net.cor.p_df, diag = TRUE)] <- NA

net.cor_df <- cbind.data.frame(melt(net.cor.r_df[15:2912,1:14]),
                               melt(net.cor.p_df[15:2912,1:14]))[,c(1:3,6)]

colnames(net.cor_df) <- c("Ortholog", "Metabolite", "cor.value", "p.value")
net.cor_df$p.adjust <- p.adjust(net.cor_df$p.value, method = "fdr")

net.cor_df <- subset(net.cor_df, p.adjust <= 0.05)

#
#### DRAW BIPARTITE NETWORK - POSITIVE ####
ko.sgm_net <- graph_from_data_frame(as.matrix(subset(net.cor_df, cor.value > 0)[,c(1:2)]), directed = FALSE)

set.seed(1)
l <- layout_with_fr(ko.sgm_net, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
l_df <- as.data.frame(l)  ## convert the layout to a data.frame
l_df$names <- names(V(ko.sgm_net))  ## add in the species codes

ko.sgm_cwt <- cluster_walktrap(ko.sgm_net)
max(ko.sgm_cwt$membership)

## NODE TABLE
ko.sgm.node_df <- data.frame(names = ko.sgm_cwt$names, membership = ko.sgm_cwt$membership)
ko.sgm.node_df <- merge(ko.sgm.node_df, sgm_group, by.x = "names", by.y = "variable", all.x = TRUE)
ko.sgm.node_df$membership <- ifelse(is.na(ko.sgm.node_df$cols), ko.sgm.node_df$membership, 
                                    max(ko.sgm.node_df$membership) + ko.sgm.node_df$cols)
ko.sgm.node_df$Group <- factor(ko.sgm.node_df$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))
row.names(ko.sgm.node_df) <- ko.sgm.node_df$names

ko.sgm.node_df <- merge(ko.sgm.node_df, l_df, by = "names")
ko.sgm.node_df$Type <- ifelse(is.na(ko.sgm.node_df$Group), "Ortholog", "Metabolite")
ko.sgm.node_df <- ko.sgm.node_df[order(ko.sgm.node_df$membership, decreasing = FALSE),]

## EDGE TABLE
ko.sgm.edge_df <- get.data.frame(ko.sgm_net)
ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,c(1,5,6)], by.x = "from", by.y = "names")
ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,c(1,3,5,6)], by.x = "to", by.y = "names")
colnames(ko.sgm.edge_df) <- c("Metabolite", "Ortholog", "Metabolite.x", "Metabolite.y", "Group", "Ortholog.x", "Ortholog.y")

ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,1:2], by.x = "Metabolite", by.y = "names")
ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,1:2], by.x = "Ortholog", by.y = "names")

## DRAW NETWORK
gg.bipnet <- ggplot() +
  geom_curve(data = ko.sgm.edge_df, aes(x = Metabolite.x, xend = Ortholog.x, y = Metabolite.y, yend = Ortholog.y), 
             color = "gray70", linewidth = 1.25, alpha = 0.8) +
  geom_point(data = ko.sgm.node_df, aes(x = V1, y = V2, shape = Type, size = Type, alpha = Type, fill = factor(membership)), 
             stroke = 1.25, color = "black", show.legend = FALSE) +
  geom_text(data = subset(ko.sgm.node_df, membership >= 6), aes(x = V1, y = V2, label = names), size = 10) +
  scale_alpha_manual(values = c(1, 0.6)) +
  scale_size_manual(values = c(20, 6)) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c("darkorchid", "deepskyblue", "sienna1", "yellowgreen", "#910a25", 
                               "#bbcc06", "#81c236", "#30b3bf", "#82107a")) +
  theme_void() 

gg.bipnet

#
#### MODULES ABUNDANCE - GLOBAL (genus) ####
## DEEPSKYBLUE MODULE - Sugars, Fusel.alcohol.acetates

ko.sgm.mod1_df <- ko.t_df[, colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 2)$names]
ko.sgm.mod1_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod1_df), Module.Ab = rowSums(ko.sgm.mod1_df),
                                   sample_df[row.names(ko.sgm.mod1_df), 
                                              colnames(sample_df) %in% c("Sugars", "Fusel.alcohol.acetates")])

cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Sugars)
cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Fusel.alcohol.acetates)

#
## YELLOWGREEN MODULE - Tartaric_acid, Acetic_acid, EEFA

ko.sgm.mod2_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 4)$names]
ko.sgm.mod2_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod2_df), Module.Ab = rowSums(ko.sgm.mod2_df),
                                   sample_df[row.names(ko.sgm.mod2_df),
                                             colnames(sample_df) %in% c("Tartaric_acid", "Acetic_acid", "EEFA")])

cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Tartaric_acid)
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Acetic_acid)
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$EEFA)

#
## REDBROWN (#910a25) MODULE - Ethanol, Glycerol

ko.sgm.mod3_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 5)$names]
ko.sgm.mod3_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod3_df), Module.Ab = rowSums(ko.sgm.mod3_df),
                                   sample_df[row.names(ko.sgm.mod3_df),
                                             colnames(sample_df) %in% c("Ethanol", "Glycerol")])

cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Ethanol)
cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Glycerol)

#
## YELLOWGREEN MODULE - Lactic_acid, Succinic_acid, Fusel.alcohols

ko.sgm.mod4_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 3)$names]
ko.sgm.mod4_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod4_df), Module.Ab = rowSums(ko.sgm.mod4_df),
                                   sample_df[row.names(ko.sgm.mod4_df),
                                              colnames(sample_df) %in% c("Lactic_acid", "Succinic_acid", "Fusel.alcohols")])

cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Lactic_acid)
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Succinic_acid)
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Fusel.alcohols)

#
## Draw boxplot

ko.sgm.mod_df <- rbind(cbind.data.frame(melt(ko.sgm.mod1_df, id.vars = c("Sample_ID", "Module.Ab")), Module = "Module 1"),
                       cbind.data.frame(melt(ko.sgm.mod2_df, id.vars = c("Sample_ID", "Module.Ab")), Module = "Module 2"),
                       cbind.data.frame(melt(ko.sgm.mod3_df, id.vars = c("Sample_ID", "Module.Ab")), Module = "Module 3"),
                       cbind.data.frame(melt(ko.sgm.mod4_df, id.vars = c("Sample_ID", "Module.Ab")), Module = "Module 4"))

ko.sgm.mod_df <- merge(ko.sgm.mod_df, sample_df[,c(1,23)], by = "Sample_ID")
ko.sgm.mod_df$Genus <- factor(ko.sgm.mod_df$Genus, levels = c("Saccharomyces", "Lachancea", "Hanseniaspora", "Other"))

gg.mod_gen <- ggplot(ko.sgm.mod_df) +
  geom_point(aes(x = Module.Ab, y = value, color = Genus), size = 3) +
  scale_color_manual(values = c("#1b9e77", "#8da0cb", "#cc3939", "gray80"), name = "Dominant Genus") +
  ylab("Metabolite Concentration") + xlab("Normalized Module Expression") + 
  facet_wrap(~ Module + variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.mod_gen

## Set Strip background colors
gg.mod_gen <- ggplot_gtable(ggplot_build(gg.mod_gen))

strip_both <- which(grepl("strip-", gg.mod_gen$layout$name))
fills <- alpha(c("#910a25", "#910a25", "sienna1", "sienna1", "sienna1",
                 "deepskyblue", "deepskyblue", "yellowgreen", "yellowgreen", "yellowgreen"), alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.mod_gen$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid::grid.draw(gg.mod_gen)

#
#### EXPORT FIGURE 6 ####

gg.figure6 <- plot_grid(gg.bipnet, gg.mod_gen, rel_heights = c(1.35, 1), ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figure6

ggsave("Figures/Figure_6.png", gg.figure6, bg = "white", width = 12.6, height = 14)


#