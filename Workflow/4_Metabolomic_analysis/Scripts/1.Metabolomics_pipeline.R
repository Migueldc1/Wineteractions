# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(vegan)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(igraph)
library(grid)
library(ggplotify)
library(goseq)
library(cowplot)

#library(ggforce)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

## Load Environment
#load("Workflow/4_Metabolomic_analysis/Outputs/RData/Wine_net.RData")
#save.image("Workflow/4_Metabolomic_analysis/Outputs/RData/Wine_net.RData")

#
#### CUSTOM FUNCTION ####

#
#### LOAD DATA ####

## KEGG ORTHOLOGY TABLE
ko.n_df <- readRDS("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/ko.n_df.rds")
ko.deg_df <- readRDS("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/ko.deg_df.rds")

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

sample_df$Dom.Genus <- ifelse(sample_df$Genus %in% c("Citeromyces", "Kluyveromyces"), "Other", sample_df$Genus)
sample_df$Dom.Genus <- factor(sample_df$Dom.Genus, levels = c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"))

## COLORS
col.cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col.gen <- c("#cc3939", "#8da0cb", "#1b9e77", "gray70")


## Ecosystem Services
sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(1, rep(3, 5), 2, 1, rep(4, 6)))
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
sgm_pca.plot <- merge(sample_df[,c(1:4,25)], sgm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

# Vectors
sgm_pca.vplot <- scores(sgm_pca, display = "species")
sgm_pca.vplot <- cbind.data.frame(variable = row.names(sgm_pca.vplot), sgm_pca.vplot)
sgm_pca.vplot <- merge(sgm_group, sgm_pca.vplot[,1:3], by = "variable")
sgm_pca.vplot$Group <- factor(sgm_pca.vplot$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))

## PCA Plot (all points)
gg.sgm_pca <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Dom.Genus), shape = 21, size = 3) +
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

ko_df.t <- as.data.frame(t(ko.n_df))

ko.sgm_df <- merge(sgm.pca_df2[,c(1,4:8,11:18)], ko_df.t, by = "row.names")
row.names(ko.sgm_df) <- ko.sgm_df[,1]
ko.sgm_df <- ko.sgm_df[,-1]

ko.sgm_cor <- rcorr(as.matrix(ko.sgm_df), type = "spearman")

ko.sgm_corr <- ko.sgm_cor$r
ko.sgm_corr[upper.tri(ko.sgm_corr, diag = TRUE)] <- NA
ko.sgm_corp <- ko.sgm_cor$P
ko.sgm_corp[upper.tri(ko.sgm_corp, diag = TRUE)] <- NA

ko.sgm_df <- cbind.data.frame(melt(ko.sgm_corr[15:2912,1:14]),
                              melt(ko.sgm_corp[15:2912,1:14]))[,c(1:3,6)]

colnames(ko.sgm_df) <- c("Source", "Target", "value", "p.value")
ko.sgm_df <- ko.sgm_df[complete.cases(ko.sgm_df),]
ko.sgm_df$p.adjust <- p.adjust(ko.sgm_df$p.value, method = "fdr")

ko.sgm_df <- subset(ko.sgm_df, p.adjust <= 0.05)
ko.sgm_df$Interaction <- ifelse(ko.sgm_df$value > 0, "Positive", "Negative")

#### POSITIVE NETWORK ####
ko.sgm.pos_net <- graph_from_data_frame(as.matrix(subset(ko.sgm_df, Interaction == "Positive")[,c(1:2)]), directed = FALSE)

set.seed(1)
l <- layout_with_fr(ko.sgm.pos_net, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

ko.sgm.pos_cwt <- cluster_walktrap(ko.sgm.pos_net)
max(ko.sgm.pos_cwt$membership)

ko.sgm.node_df <- data.frame(names = ko.sgm.pos_cwt$names, membership = ko.sgm.pos_cwt$membership)
ko.sgm.node_df <- merge(ko.sgm.node_df, ko.deg_df[,c(1,4)], by.x = "names", by.y = "KEGG_ko", all.x = TRUE)
ko.sgm.node_df <- merge(ko.sgm.node_df, sgm_group, by.x = "names", by.y = "variable", all.x = TRUE)
ko.sgm.node_df$membership <- ifelse(is.na(ko.sgm.node_df$cols), ko.sgm.node_df$membership, 
                                    max(ko.sgm.node_df$membership) + ko.sgm.node_df$cols)

ko.sgm.node_df$Group <- ifelse(is.na(ko.sgm.node_df$cols), ko.sgm.node_df$OverExpr, ko.sgm.node_df$Group)
ko.sgm.node_df$Group[is.na(ko.sgm.node_df$Group)] <- "none"
ko.sgm.node_df$Group <- factor(ko.sgm.node_df$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles", "Saccharomyces",
                                                                "Lachancea", "Hanseniaspora", "Lt&Hs", "none"))

ko.sgm.node_df <- ko.sgm.node_df[,c(1,2,4)]
row.names(ko.sgm.node_df) <- ko.sgm.node_df$names

V(ko.sgm.pos_net)$community <- ko.sgm.node_df[V(ko.sgm.pos_net)$name,"membership"]
V(ko.sgm.pos_net)$Group <- as.numeric(ko.sgm.node_df[V(ko.sgm.pos_net)$name,"Group"])

colrs <- adjustcolor(c("darkorchid", "deepskyblue", "sienna1", "yellowgreen", "#a54a2a", "#910a25", 
                       "yellow", "#b2ed6d", "#375bb0", "#d128c5"), alpha = 0.8)
                                   
colrs1 <- adjustcolor(c("yellow", "#b2ed6d", "#375bb0", "#d128c5", "#1b9e77", "#8da0cb", "#cc3939", "#AD6D82", "gray70"), 
                      alpha = 0.8)
                                   
gg.ko.sgm_pos <- as.grob(~plot(ko.sgm.pos_net, vertex.size = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, 8, 5), 
                    vertex.label = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, names(V(ko.sgm.pos_net)), NA), 
                    vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
                    vertex.shape = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, "square", "circle"),
                    vertex.color = colrs[V(ko.sgm.pos_net)$community], 
                    rescale = FALSE, layout = l*1.3, edge.color = "gray70",
                    edge.curved = 0.5))

gg.ko.sgm_pos <- ggplotify::as.ggplot(gg.ko.sgm_pos)
gg.ko.sgm_pos

ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/ko.sgm_pos.png", gg.ko.sgm_pos, bg = "white",
       width = 11, height = 7)

gg.ko.sgm_pos1 <- as.grob(~plot(ko.sgm.pos_net, vertex.size = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, 8, 5), 
                    vertex.label = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, names(V(ko.sgm.pos_net)), NA), 
                    vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
                    vertex.shape = ifelse(names(V(ko.sgm.pos_net)) %in% sgm_group$variable, "square", "circle"),
                    vertex.color = colrs1[as.numeric(V(ko.sgm.pos_net)$Group)], 
                    rescale = FALSE, layout = l*1.3, edge.color = "gray70",
                    edge.curved = 0.5))

gg.ko.sgm_pos1 <- ggplotify::as.ggplot(gg.ko.sgm_pos1)
gg.ko.sgm_pos1
ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/ko.sgm_pos1.png", gg.ko.sgm_pos1, bg = "white",
       width = 11, height = 7)

#
#### NEGATIVE NETWORK ####
ko.sgm.neg_net <- graph_from_data_frame(as.matrix(subset(ko.sgm_df, Interaction == "Negative")[,c(1:2)]), directed = FALSE)

set.seed(1)
l <- layout_with_fr(ko.sgm.neg_net, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

ko.sgm.neg_cwt <- cluster_walktrap(ko.sgm.neg_net)
max(ko.sgm.neg_cwt$membership)

ko.sgm.node.n_df <- data.frame(names = ko.sgm.neg_cwt$names, membership = ko.sgm.neg_cwt$membership)
ko.sgm.node.n_df <- merge(ko.sgm.node.n_df, ko.deg_df[,c(1,4)], by.x = "names", by.y = "KEGG_ko", all.x = TRUE)
ko.sgm.node.n_df <- merge(ko.sgm.node.n_df, sgm_group, by.x = "names", by.y = "variable", all.x = TRUE)
ko.sgm.node.n_df$membership <- ifelse(is.na(ko.sgm.node.n_df$cols), ko.sgm.node.n_df$membership, 
                                    max(ko.sgm.node.n_df$membership) + ko.sgm.node.n_df$cols)

ko.sgm.node.n_df$Group <- ifelse(is.na(ko.sgm.node.n_df$cols), ko.sgm.node.n_df$OverExpr, ko.sgm.node.n_df$Group)
ko.sgm.node.n_df$Group[is.na(ko.sgm.node.n_df$Group)] <- "none"
ko.sgm.node.n_df$Group <- factor(ko.sgm.node.n_df$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles", "Saccharomyces",
                                                                "Lachancea", "Hanseniaspora", "Lt&Hs", "none"))

ko.sgm.node.n_df <- ko.sgm.node.n_df[,c(1,2,4)]
row.names(ko.sgm.node.n_df) <- ko.sgm.node.n_df$names

V(ko.sgm.neg_net)$community <- ko.sgm.node.n_df[V(ko.sgm.neg_net)$name,"membership"]
V(ko.sgm.neg_net)$Group <- as.numeric(ko.sgm.node.n_df[V(ko.sgm.neg_net)$name,"Group"])

colrs <- adjustcolor(c("darkorchid", "deepskyblue", "sienna1", "yellowgreen", "#a52a2a", 
                                   "yellow", "#b2ed6d", "#375bb0", "#d128c5"), alpha = 0.8)
                                   
colrs1 <- adjustcolor(c("yellow", "#b2ed6d", "#375bb0", "#d128c5", "#1b9e77", "#8da0cb", "#cc3939", "#AD6D82", "gray70"), 
                      alpha = 0.8)

gg.ko.sgm_neg <- as.grob(~plot(ko.sgm.neg_net, vertex.size = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, 8, 5), 
                               vertex.label = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, names(V(ko.sgm.neg_net)), NA), 
                               vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
                               vertex.shape = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, "square", "circle"),
                               vertex.color = colrs[V(ko.sgm.neg_net)$community], 
                               rescale = FALSE, layout = l*1.3, edge.color = "gray70",
                               edge.curved = 0.5))

gg.ko.sgm_neg <- ggplotify::as.ggplot(gg.ko.sgm_neg)
gg.ko.sgm_neg
ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/ko.sgm_neg.png", gg.ko.sgm_neg, bg = "white",
       width = 11, height = 7)

gg.ko.sgm_neg1 <- as.grob(~plot(ko.sgm.neg_net, vertex.size = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, 8, 5), 
                                vertex.label = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, names(V(ko.sgm.neg_net)), NA), 
                                vertex.label.cex = 1.5, vertex.label.color = "black", vertex.label.family = "sans",
                                vertex.shape = ifelse(names(V(ko.sgm.neg_net)) %in% sgm_group$variable, "square", "circle"),
                                vertex.color = colrs1[as.numeric(V(ko.sgm.neg_net)$Group)], 
                                rescale = FALSE, layout = l*1.3, edge.color = "gray70",
                                edge.curved = 0.5))

gg.ko.sgm_neg1 <- ggplotify::as.ggplot(gg.ko.sgm_neg1)
gg.ko.sgm_neg1
ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/ko.sgm_neg1.png", gg.ko.sgm_neg1, bg = "white",
       width = 11, height = 7)

#
############################################## METABOLITE BIOLOGICAL PROCESSES ####
ko.sgm.tot_df <- merge(ko.sgm.node_df, ko.sgm.node.n_df, by = "names", all = TRUE)
ko.sgm.tot_df$Group <- ifelse(is.na(ko.sgm.tot_df$Group.x), 
                              as.character(ko.sgm.tot_df$Group.y), as.character(ko.sgm.tot_df$Group.x))

ko.sgm.tot_df <- cbind.data.frame(ko.sgm.tot_df[,c(1,6)], 
                                  Positive = ko.sgm.tot_df$membership.x, Negative = ko.sgm.tot_df$membership.y)

go_df <- readRDS("Workflow/3_Meta-transcriptomics_analysis/Outputs/RData/go_df.rds")

#
#### POSITIVE MODULES ####

count.bias <- rowMeans(ko_df[,-1])
names(count.bias) <- ko_df[,1]

## DARKORCHID MODULE - acetic_acid, tartaric_acid, EEFA, fusel.alcohol.acetates

ko.sgm.mod1_df <- subset(ko.sgm.node_df, membership == 1)

DEG_mod1 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.mod1_df$names, 1, 0)
names(DEG_mod1) <- ko_df$KEGG_ko

pwf_mod1 <- nullp(DEgenes = DEG_mod1, bias.data = jitter(rep(1000, length(DEG_mod1))))
go_mod1 <- goseq(pwf_mod1, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO_mod1 <- subset(go_mod1, ontology == "BP" & numDEInCat > 5)

## DEEPSKYBLUE MODULE - ethanol, glycerol, ethyl.acetate

ko.sgm.mod2_df <- subset(ko.sgm.node_df, membership == 2)

DEG_mod2 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.mod2_df$names, 1, 0)
names(DEG_mod2) <- ko_df$KEGG_ko

pwf_mod2 <- nullp(DEgenes = DEG_mod2, bias.data = jitter(rep(1000, length(DEG_mod2))))
go_mod2 <- goseq(pwf_mod2, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO_mod2 <- subset(go_mod2, ontology == "BP" & numDEInCat > 5)

## SIENNA MODULE - succinic_acid, lactic_acid, fusel.alcohols

ko.sgm.mod3_df <- subset(ko.sgm.node_df, membership == 3)

DEG_mod3 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.mod3_df$names, 1, 0)
names(DEG_mod3) <- ko_df$KEGG_ko

pwf_mod3 <- nullp(DEgenes = DEG_mod3, bias.data = jitter(rep(1000, length(DEG_mod3))))
go_mod3 <- goseq(pwf_mod3, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO_mod3 <- subset(go_mod3, ontology == "BP" & numDEInCat > 5)

## GREEN (#349e87) MODULE - sugars, acetic_acid

ko.sgm.mod4_df <- subset(ko.sgm.node_df, membership == 7)

DEG_mod4 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.mod4_df$names, 1, 0)
names(DEG_mod4) <- ko_df$KEGG_ko

pwf_mod4 <- nullp(DEgenes = DEG_mod4, bias.data = jitter(rep(1000, length(DEG_mod4))))
go_mod4 <- goseq(pwf_mod4, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO_mod4 <- subset(go_mod4, ontology == "BP" & numDEInCat > 5)

# Plot
go_modpos <- rbind.data.frame(cbind(enrichGO_mod1[1:20,], Module = "Purple"),
                              cbind(enrichGO_mod2[1:20,], Module = "Blue"),
                              cbind(enrichGO_mod3[1:20,], Module = "Brown"),
                              cbind(enrichGO_mod4[1:20,], Module = "Green"))

go_modpos$Module <- factor(go_modpos$Module, levels = c("Purple", "Blue", "Brown", "Green"))

go_modpos$comp.cat <- ifelse(duplicated(go_modpos$term) | duplicated(go_modpos$term, fromLast = TRUE), 
                             "Both", go_modpos$Module)

go_modpos$term <- factor(go_modpos$term, 
                         levels = unique(go_modpos$term[order(go_modpos$numDEInCat, go_modpos$comp.cat, decreasing = TRUE)]))

gg.go_modpos <- ggplot(go_modpos) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Module), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("darkorchid", "deepskyblue", "sienna1", "#349e87")) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        aspect.ratio = 1.2,
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_modpos

ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/net.pos_enrichGO.png", gg.go_modpos, 
       width = 18, height = 11, dpi = 300, bg = "white")

#
#### NEGATIVE MODULES ####

## DEEPSKYBLUE MODULE - ethanol, glycerol, ethyl.acetate

ko.sgm.n.mod1_df <- subset(ko.sgm.node.n_df, membership == 2)

DEG.n_mod1 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.n.mod1_df$names, 1, 0)
names(DEG.n_mod1) <- ko_df$KEGG_ko

pwf.n_mod1 <- nullp(DEgenes = DEG.n_mod1, bias.data = jitter(rep(1000, length(DEG.n_mod1))))
go.n_mod1 <- goseq(pwf.n_mod1, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO.n_mod1 <- subset(go.n_mod1, ontology == "BP" & numDEInCat > 5)

## SIENNA MODULE - acetic_acid, tartaric_acid, EEFA, fusel.alcohols

ko.sgm.n.mod2_df <- subset(ko.sgm.node.n_df, membership == 3)

DEG.n_mod2 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.n.mod2_df$names, 1, 0)
names(DEG.n_mod2) <- ko_df$KEGG_ko

pwf.n_mod2 <- nullp(DEgenes = DEG.n_mod2, bias.data = jitter(rep(1000, length(DEG.n_mod2))))
go.n_mod2 <- goseq(pwf.n_mod2, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO.n_mod2 <- subset(go.n_mod2, ontology == "BP" & numDEInCat > 5)

## YELLOWGREEN MODULE - succinic_acid, lactic_acid, fusel.alcohols, MCFA

ko.sgm.n.mod3_df <- subset(ko.sgm.node.n_df, membership == 4)

DEG.n_mod3 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.n.mod3_df$names, 1, 0)
names(DEG.n_mod3) <- ko_df$KEGG_ko

pwf.n_mod3 <- nullp(DEgenes = DEG.n_mod3, bias.data = jitter(rep(1000, length(DEG.n_mod3))))
go.n_mod3 <- goseq(pwf.n_mod3, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO.n_mod3 <- subset(go.n_mod3, ontology == "BP" & numDEInCat > 5)

## REDBROWN (#a52a2a) MODULE - sugars

ko.sgm.n.mod4_df <- subset(ko.sgm.node.n_df, membership == 5)

DEG.n_mod4 <- ifelse(ko_df$KEGG_ko %in% ko.sgm.n.mod4_df$names, 1, 0)
names(DEG.n_mod4) <- ko_df$KEGG_ko

pwf.n_mod4 <- nullp(DEgenes = DEG.n_mod4, bias.data = jitter(rep(1000, length(DEG.n_mod4))))
go.n_mod4 <- goseq(pwf.n_mod4, gene2cat = go_df, test.cats = c("GO:BP"))

enrichGO.n_mod4 <- subset(go.n_mod4, ontology == "BP" & numDEInCat > 5)

# Plot
go.n_modpos <- rbind.data.frame(cbind(enrichGO.n_mod1[1:20,], Module = "Blue"),
                                cbind(enrichGO.n_mod2[1:20,], Module = "Brown"),
                                cbind(enrichGO.n_mod3[1:20,], Module = "Green"),
                                cbind(enrichGO.n_mod4[1:20,], Module = "RedBrown"))

go.n_modpos$Module <- factor(go.n_modpos$Module, levels = c("Blue", "Brown", "Green", "RedBrown"))

go.n_modpos$comp.cat <- ifelse(duplicated(go.n_modpos$term) | duplicated(go.n_modpos$term, fromLast = TRUE), 
                               "Both", go.n_modpos$Module)

go.n_modpos$term <- factor(go.n_modpos$term, 
                           levels = unique(go.n_modpos$term[order(go.n_modpos$numDEInCat, go.n_modpos$comp.cat, 
                                                                  decreasing = TRUE)]))

gg.go.n_modpos <- ggplot(go.n_modpos) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Module), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("deepskyblue", "sienna1", "#349e87", "yellowgreen", "#a52a2a")) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        aspect.ratio = 1.2,
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go.n_modpos

ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/net.neg_enrichGO.png", gg.go.n_modpos, 
       width = 18, height = 13, dpi = 300, bg = "white")

#
############################################## METABOLITE MODULES ####
#
#### POSITIVE MODULES ####
## DARKORCHID MODULE - Sugars, Fusel.alcohol.acetates

ko.sgm.mod1_df <- ko_df.t[, colnames(ko_df.t) %in% subset(ko.sgm.node_df, membership == 1)$names]
ko.sgm.mod1_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.mod1_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Dom.Genus", "Sugars", "Fusel.alcohol.acetates")])

cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Sugars)
cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Fusel.alcohol.acetates)

gg.mod1_gen <- ggplot(ko.sgm.mod1_df) +
  geom_boxplot(aes(x = Dom.Genus, y = Module.Ab, color = Dom.Genus), size = 1.5) +
  scale_color_manual(values = col.gen) +
  xlab("Dominant Genus") + ylab("Normalized Module Counts") + 
  annotate(geom = "text", x = 2.5, y = Inf, label = "Module 1", size = 6, vjust = 3, color = "darkorchid", fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1.75,
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"))

gg.mod1_gen

ggplot(ko.sgm.mod1_df) +
  geom_point(aes(x = Module.Ab, y = Sugars, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod1_df) +
  geom_point(aes(x = Module.Ab, y = Fusel.alcohol.acetates, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

#
## SIENNA MODULE - Lactic_acid, Succinic_acid, Fusel.alcohols

ko.sgm.mod2_df <- ko_df.t[,colnames(ko_df.t) %in% subset(ko.sgm.node_df, membership == 3)$names]
ko.sgm.mod2_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.mod2_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Dom.Genus", "Lactic_acid", "Succinic_acid", "Fusel.alcohols")])

cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Lactic_acid)
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Succinic_acid)
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Fusel.alcohols)

gg.mod2_gen <- ggplot(ko.sgm.mod2_df) +
  geom_boxplot(aes(x = Dom.Genus, y = Module.Ab, color = Dom.Genus), size = 1.5) +
  scale_color_manual(values = col.gen) +
  xlab("Dominant Genus") + ylab("Normalized Module Counts") + 
  annotate(geom = "text", x = 2.5, y = Inf, label = "Module 2", size = 6, vjust = 3, color = "sienna1", fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1.75,
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"))

gg.mod2_gen

ggplot(ko.sgm.mod2_df) +
  geom_point(aes(x = Module.Ab, y = Lactic_acid, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod2_df) +
  geom_point(aes(x = Module.Ab, y = Succinic_acid, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod2_df) +
  geom_point(aes(x = Module.Ab, y = Fusel.alcohols, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

#
## YELLOWGREEN MODULE - Ethanol, Glycerol

ko.sgm.mod3_df <- ko_df.t[,colnames(ko_df.t) %in% subset(ko.sgm.node_df, membership == 4)$names]
ko.sgm.mod3_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.mod3_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Dom.Genus", "Ethanol", "Glycerol")])

cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Ethanol)
cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Glycerol)

gg.mod3_gen <- ggplot(ko.sgm.mod3_df) +
  geom_boxplot(aes(x = Dom.Genus, y = Module.Ab, color = Dom.Genus), size = 1.5) +
  scale_color_manual(values = col.gen) +
  xlab("Dominant Genus") + ylab("Normalized Module Counts") + 
  annotate(geom = "text", x = 2.5, y = Inf, label = "Module 3", size = 6, vjust = 3, color = "yellowgreen", fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1.75,
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"))

gg.mod3_gen

ggplot(ko.sgm.mod3_df) +
  geom_point(aes(x = Module.Ab, y = Ethanol, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod3_df) +
  geom_point(aes(x = Module.Ab, y = Glycerol, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

#
## REDBROWN (#910a25) MODULE - Tartaric_acid, Acetic_acid, EEFA

ko.sgm.mod4_df <- ko_df.t[,colnames(ko_df.t) %in% subset(ko.sgm.node_df, membership == 6)$names]
ko.sgm.mod4_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.mod4_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Dom.Genus", "Tartaric_acid", "Acetic_acid", "EEFA")])

cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Tartaric_acid)
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Acetic_acid)
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$EEFA)

gg.mod4_gen <- ggplot(ko.sgm.mod4_df) +
  geom_boxplot(aes(x = Dom.Genus, y = Module.Ab, color = Dom.Genus), size = 1.5) +
  scale_color_manual(values = col.gen) +
  xlab("Dominant Genus") + ylab("Normalized Module Counts") + 
  annotate(geom = "text", x = 2.5, y = Inf, label = "Module 4", size = 6, vjust = 3, color = "#910a25", fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1.75,
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"))

gg.mod4_gen

ggplot(ko.sgm.mod4_df) +
  geom_point(aes(x = Lactic_acid, y = Module.Ab, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod4_df) +
  geom_point(aes(x = Succinic_acid, y = Module.Ab, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

ggplot(ko.sgm.mod4_df) +
  geom_point(aes(x = Fusel.alcohols, y = Module.Ab, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen)

#
#### GENERATE BOXPLOT ####
ko.sgm.mod_df <- rbind.data.frame(cbind(ko.sgm.mod1_df[,c(1,4)], module = "Module 1"),
                                  cbind(ko.sgm.mod2_df[,c(1,5)], module = "Module 2"), 
                                  cbind(ko.sgm.mod3_df[,c(1,4)], module = "Module 3"), 
                                  cbind(ko.sgm.mod4_df[,c(1,5)], module = "Module 4"))

gg.mod_gen <- ggplot(ko.sgm.mod_df) +
  geom_boxplot(aes(x = Dom.Genus, y = Module.Ab, color = Dom.Genus), size = 1.5) +
  scale_color_manual(values = col.gen, name = "Dominant Genus") +
  xlab("Dominant Genus") + ylab("Normalized Module Counts") + 
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(aspect.ratio = 1.75,
        legend.position = "bottom",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 14, color = "black"))

gg.mod_gen

## Set Strip background colors
gg.mod_gen <- ggplot_gtable(ggplot_build(gg.mod_gen))

strip_both <- which(grepl("strip-", gg.mod_gen$layout$name))
fills <- alpha(c("darkorchid", "sienna1", "yellowgreen", "#910a25"), alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.mod_gen$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid.draw(gg.mod_gen)
ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/mod.ab_gen.png", gg.mod_gen, 
       width = 16, height = 7)

#
#### GENERATE CORRELATIONS ####

ko.sgm.cor_df <- rbind(cbind(ko.sgm.mod1_df[,c(1,4)], value = ko.sgm.mod1_df[,2], metabolite = "Sugars", module = "Module 1"),
                       cbind(ko.sgm.mod1_df[,c(1,4)], value = ko.sgm.mod1_df[,3], metabolite = "Fusel.alcohol.acetates", module = "Module 1"), 
                       cbind(ko.sgm.mod2_df[,c(1,5)], value = ko.sgm.mod2_df[,2], metabolite = "Lactic_acid", module = "Module 2"), 
                       cbind(ko.sgm.mod2_df[,c(1,5)], value = ko.sgm.mod2_df[,3], metabolite = "Succinic_acid", module = "Module 2"),
                       cbind(ko.sgm.mod2_df[,c(1,5)], value = ko.sgm.mod2_df[,4], metabolite = "Fusel.alcohols", module = "Module 2"),
                       cbind(ko.sgm.mod3_df[,c(1,4)], value = ko.sgm.mod3_df[,2], metabolite = "Ethanol", module = "Module 3"),
                       cbind(ko.sgm.mod3_df[,c(1,4)], value = ko.sgm.mod3_df[,3], metabolite = "Glycerol", module = "Module 3"),
                       cbind(ko.sgm.mod4_df[,c(1,5)], value = ko.sgm.mod4_df[,2], metabolite = "Acetic_acid", module = "Module 4"), 
                       cbind(ko.sgm.mod4_df[,c(1,5)], value = ko.sgm.mod4_df[,3], metabolite = "Tartaric_acid", module = "Module 4"),
                       cbind(ko.sgm.mod4_df[,c(1,5)], value = ko.sgm.mod4_df[,4], metabolite = "EEFA", module = "Module 4"))

gg.cor_gen <- ggplot(ko.sgm.cor_df) +
  geom_point(aes(x = Module.Ab, y = value, color = Dom.Genus), size = 3) +
  scale_color_manual(values = col.gen, name = "Dominant Genus") +
  xlab("Normalized Module Counts") + ylab("Metabolite Concentration") + 
  facet_wrap(~ module + metabolite, scales = "free", nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 14, color = "black"))

gg.cor_gen

## Set Strip background colors
gg.cor_gen <- ggplot_gtable(ggplot_build(gg.cor_gen))

strip_both <- which(grepl("strip-", gg.cor_gen$layout$name))
fills <- alpha(c(rep("yellowgreen", 2), rep("#910a25", 3), rep("darkorchid", 2), rep("sienna1", 3)), alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.cor_gen$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid.draw(gg.cor_gen)
ggsave("Workflow/4_Metabolomic_analysis/Outputs/Figures/mod.cor_gen.png", gg.cor_gen, 
       width = 16, height = 8)

#





#### BIOLOGICAL ENRICHMENT - MODULE 1 ####


##KEGG
kk.AHLA <- enrichKEGG(gene = names(genes.AHLA[genes.AHLA == 1]), organism = 'pae')
head(kk.AHLA, n=10)

browseKEGG(kk.AHLA, 'pae00405')
browseKEGG(kk.AHLA, 'pae02024')
browseKEGG(kk.AHLA, 'pae02025')
browseKEGG(kk.AHLA, 'pae00460')

logFC <- res.AHLA$log2FoldChange
names(logFC) <- row.names(res.AHLA)

logFC <- logFC[names(logFC) %in% names(genes.AHLA[genes.AHLA == 1])]
pathview(gene.data = logFC, 
         pathway.id = "pae02024", 
         species = "pae",  
         gene.idtype="kegg")
#




#### NEGATIVE MODULES ####
## DEEPSKYBLUE MODULE - Ethanol, Glycerol, Ethyl.acetate

ko.sgm.n.mod1_df <- ko_df.t[,ko_df$KEGG_ko %in% subset(ko.sgm.node.n_df, membership == 2)$names]
ko.sgm.n.mod1_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.n.mod1_df),
                                     sample_df[,colnames(sample_df) %in%
                                                 c("Genus", "Ethanol", "Glycerol", "Ethyl.acetate")])

summary(lm(Ethyl.acetate ~ Module.Ab:Genus, ko.sgm.n.mod1_df))

ggplot(ko.sgm.n.mod1_df) +
  geom_point(aes(x = Module.Ab, y = Glycerol, color = Genus)) +
  geom_smooth(aes(x = Module.Ab, y = Glycerol, color = Genus), method = "lm", se = FALSE)


## Con el ethanol, solo sale significativa la pendiente de -Hanseniaspora y +Saccharomyces (p < 0.001, Adj.R2 = 0.602)
## Glycerol El modelo ni siquiera sale significativo (p = 0.230, Adj.R2 = 0.040)
## Ethyl.acetate El modelo ni siquiera sale significativo (p = 0.718, Adj.R2 = -0.045)
#
## SIENNA MODULE - Acetic_acid, Tartaric_acid, EEFA, Fusel.alcohols

ko.sgm.n.mod2_df <- ko_df.t[,ko_df$KEGG_ko %in% subset(ko.sgm.node_df, membership == 3)$names]
ko.sgm.n.mod2_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.n.mod2_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Genus", "Acetic_acid", "Tartaric_acid", "EEFA", "Fusel.alcohols")])




summary(lm(Fusel.alcohols ~ Module.Ab:Genus, ko.sgm.n.mod2_df))

ggplot(ko.sgm.n.mod2_df) +
  geom_point(aes(x = Module.Ab, y = Ethyl.acetate, color = Genus)) +
  geom_smooth(aes(x = Module.Ab, y = Ethyl.acetate, color = Genus), method = "lm", se = FALSE)


## Con el Acetic_acid, solo sale significativa la pendiente de Lachancea (p < 0.001, Adj.R2 = 0.318)
## Con el Tartaric_acid, solo sale significativa la pendiente de Lachancea (p = 0.026, Adj.R2 = 0.146)
## Para el EEFA nada (p = 0.198, Adj.R2 = 0.053)
## Fusel.alcohols, solo sale significativa la pendiente de Kluyveromyces (p = 0.063, Adj.R2 = 0.115)
#
## YELLOWGREEN MODULE - Succinic_acid, Lactic_acid, Fusel.alcohols, MCFA

ko.sgm.n.mod3_df <- ko_df.t[,ko_df$KEGG_ko %in% subset(ko.sgm.node_df, membership == 4)$names]
ko.sgm.n.mod3_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.n.mod3_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Genus", "Succinic_acid", "Lactic_acid", "Fusel.alcohols", "MCFA")])




summary(lm(MCFA ~ Module.Ab:Genus, ko.sgm.n.mod3_df))

ggplot(ko.sgm.n.mod3_df) +
  geom_point(aes(x = Module.Ab, y = Fusel.alcohols, color = Genus)) +
  geom_smooth(aes(x = Module.Ab, y = Fusel.alcohols, color = Genus), method = "lm", se = FALSE)


## Con el Succinic_acid, solo sale significativa la pendiente de Lachancea (p = 0.001, Adj.R2 = 0.384)
## Con el Lactic_acid, solo NO sale significativa la pendiente de Lachancea (p < 0.001, Adj.R2 = 0.382)
## Con el Fusel.alcohols, solo NO sale significativa la pendiente de Lachancea (p = 0.051, Adj.R2 = 0.125)
## Con el MCFA, solo sale significativa la pendiente de Kluyveromyces (p = 0.211, Adj.R2 = 0.049)
#
## REDBROWN (#a52a2a) MODULE - Sugars

ko.sgm.n.mod4_df <- ko_df.t[,ko_df$KEGG_ko %in% subset(ko.sgm.node_df, membership == 5)$names]
ko.sgm.n.mod4_df <- cbind.data.frame(Module.Ab = rowSums(ko.sgm.n.mod4_df),
                                   sample_df[,colnames(sample_df) %in%
                                               c("Genus", "Sugars")])




summary(lm(Sugars ~ Module.Ab:Genus, ko.sgm.n.mod4_df))

ggplot(ko.sgm.n.mod4_df) +
  geom_point(aes(x = Module.Ab, y = Acetic_acid, color = Genus)) +
  geom_smooth(aes(x = Module.Ab, y = Acetic_acid, color = Genus), method = "lm", se = FALSE)


## Con el Sugars, solo sale significativa la pendiente de +Hanseniaspora y -Saccharomyces (p < 0.001, Adj.R2 = 0.383)
#
#### SUGARS ####
sugars_df <- subset(ko.sgm.tot_df, Positive == 7 | Negative == 5)
sugars_df$Interaction <- ifelse(sugars_df$Positive == 7, "Positive", "Negative")
sugars_df$Interaction[is.na(sugars_df$Interaction)] <- "Negative"

sugars_df <- sugars_df[,c(1,2,5)]








