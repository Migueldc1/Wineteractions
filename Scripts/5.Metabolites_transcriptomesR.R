# Author: M, de Celis Rodriguez
# Date: 04/07/2023
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(cowplot)
library(goseq)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - Wineteractions/GitHub/Wineteractions/")

#
#### LOAD DATA ####

## ORTHOLOG TABLES
#Normalized count table
ko.n_df <- readRDS("Data/Meta-transcriptomics/ko.n_df.rds")

kegg_df <- read.table("Data/Meta-transcriptomics/Annotation/KEGG_names.txt", sep = "\t", header = TRUE, quote = "")

go_df <- read.table("Data/Meta-transcriptomics/Annotation/go_df.txt", sep = "\t", header = TRUE)


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
################################################################################ FIGURE 5 ####
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
ko.sgm.node_df$module <- ko.sgm.node_df$membership

ko.sgm.node_df$membership <- ifelse(is.na(ko.sgm.node_df$cols), ko.sgm.node_df$membership, 
                                    max(ko.sgm.node_df$membership) + ko.sgm.node_df$cols)

ko.sgm.node_df$Group <- factor(ko.sgm.node_df$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))
row.names(ko.sgm.node_df) <- ko.sgm.node_df$names

ko.sgm.node_df <- merge(ko.sgm.node_df, l_df, by = "names")
ko.sgm.node_df$Type <- ifelse(is.na(ko.sgm.node_df$Group), "Ortholog", "Metabolite")
ko.sgm.node_df <- ko.sgm.node_df[order(ko.sgm.node_df$membership, decreasing = FALSE),]

## EDGE TABLE
ko.sgm.edge_df <- get.data.frame(ko.sgm_net)
ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,c(1,6,7)], by.x = "from", by.y = "names")
ko.sgm.edge_df <- merge(ko.sgm.edge_df, ko.sgm.node_df[,c(1,3,6,7)], by.x = "to", by.y = "names")
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
  scale_fill_manual(values = c("#9932cc", "#00bfff", "#ff8247", "#9acd32", "#910a25", 
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

cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Sugars, method = "spearman")
cor.test(ko.sgm.mod1_df$Module.Ab, ko.sgm.mod1_df$Fusel.alcohol.acetates, method = "spearman")

#
## YELLOWGREEN MODULE - Tartaric_acid, Acetic_acid, EEFA

ko.sgm.mod2_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 4)$names]
ko.sgm.mod2_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod2_df), Module.Ab = rowSums(ko.sgm.mod2_df),
                                   sample_df[row.names(ko.sgm.mod2_df),
                                             colnames(sample_df) %in% c("Tartaric_acid", "Acetic_acid", "EEFA")])

cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Acetic_acid, method = "spearman")
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$Tartaric_acid, method = "spearman")
cor.test(ko.sgm.mod2_df$Module.Ab, ko.sgm.mod2_df$EEFA, method = "spearman")

#
## REDBROWN (#910a25) MODULE - Ethanol, Glycerol

ko.sgm.mod3_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 5)$names]
ko.sgm.mod3_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod3_df), Module.Ab = rowSums(ko.sgm.mod3_df),
                                   sample_df[row.names(ko.sgm.mod3_df),
                                             colnames(sample_df) %in% c("Ethanol", "Glycerol")])

cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Ethanol, method = "spearman")
cor.test(ko.sgm.mod3_df$Module.Ab, ko.sgm.mod3_df$Glycerol, method = "spearman")

#
## SIENNA MODULE - Lactic_acid, Succinic_acid, Fusel.alcohols

ko.sgm.mod4_df <- ko.t_df[,colnames(ko.t_df) %in% subset(ko.sgm.node_df, membership == 3)$names]
ko.sgm.mod4_df <- cbind.data.frame(Sample_ID = row.names(ko.sgm.mod4_df), Module.Ab = rowSums(ko.sgm.mod4_df),
                                   sample_df[row.names(ko.sgm.mod4_df),
                                              colnames(sample_df) %in% c("Lactic_acid", "Succinic_acid", "Fusel.alcohols")])

cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Lactic_acid, method = "spearman")
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Succinic_acid, method = "spearman")
cor.test(ko.sgm.mod4_df$Module.Ab, ko.sgm.mod4_df$Fusel.alcohols, method = "spearman")

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
  geom_smooth(method = "scam",
              aes(x = Module.Ab, y = value),
              formula = y ~ s(x, k = 0), 
              se = FALSE) +
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
#### EXPORT FIGURE 5 ####

gg.figure5 <- plot_grid(gg.bipnet, gg.mod_gen, rel_heights = c(1.35, 1), ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figure5

ggsave("Figures/Inkscape/Figure_5.pdf", gg.figure5, bg = "white", width = 12.6, height = 14)


#
################################################################################ SUPPLEMENTARY TABLE S4 ####
#### ORTHOLOG - METABOLITE TABLE ####

ko.sgm.info_df <- merge(ko.sgm.node_df[,1:2], kegg_df, by.x = "names", by.y = "KEGG_ko")

ko.sgm.info_df <- merge(ko.sgm.edge_df[,c(1,2)], ko.sgm.info_df, by.x = "Ortholog", by.y = "names")

ko.sgm.info_df$Module <- ifelse(ko.sgm.info_df$membership == 2, "Module 1", 
                                ifelse(ko.sgm.info_df$membership == 4, "Module 2",
                                       ifelse(ko.sgm.info_df$membership == 5, "Module 3",
                                              ifelse(ko.sgm.info_df$membership == 3, "Module 4", NA))))

ko.sgm.info_df <- ko.sgm.info_df[complete.cases(ko.sgm.info_df), -3]
ko.sgm.info_df <- ko.sgm.info_df[order(ko.sgm.info_df$Metabolite, ko.sgm.info_df$Module),]

#
#### EXPORT TABLE S4 ####

write.table(ko.sgm.info_df, "Figures/Table_S4.txt", sep = "\t", row.names = FALSE, dec = ",")

#
################################################################################ SUPPLEMENTARY FIGURE S7 ####
#### BIOLOGICAL ENRICHMENT - MODULE 1 ####

count.bias <- rowSums(ko.n_df)

## Module 1
DEG_mod1 <- ifelse(names(count.bias) %in% subset(ko.sgm.node_df, membership == 2)$names, 1, 0)
names(DEG_mod1) <- names(count.bias)

pwf_mod1 <- nullp(DEgenes = DEG_mod1, bias.data = count.bias[names(DEG_mod1)])
go_mod1 <- goseq(pwf_mod1, gene2cat = go_df, test.cats = c("GO:BP"))
go_mod1 <- subset(go_mod1, ontology == "BP")

enrichGO_mod1 <- go_mod1[p.adjust(go_mod1$over_represented_pvalue, method = "fdr") < 0.05,]
enrichGO_mod1$term <- factor(enrichGO_mod1$term, 
                             levels = unique(enrichGO_mod1$term[order(enrichGO_mod1$numDEInCat, decreasing = TRUE)]))

gg.go_mod1 <- ggplot(enrichGO_mod1) +
  geom_bar(aes(x = term, y = numDEInCat), stat = "identity", position = "dodge", fill = "deepskyblue") + 
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_mod1

#
#### BIOLOGICAL ENRICHMENT - MODULE 2 ####

## Module 2
DEG_mod2 <- ifelse(names(count.bias) %in% subset(ko.sgm.node_df, membership == 4)$names, 1, 0)
names(DEG_mod2) <- names(count.bias)

pwf_mod2 <- nullp(DEgenes = DEG_mod2, bias.data = count.bias[names(DEG_mod2)])
go_mod2 <- goseq(pwf_mod2, gene2cat = go_df, test.cats = c("GO:BP"))
go_mod2 <- subset(go_mod2, ontology == "BP" & numDEInCat > 10)

enrichGO_mod2 <- go_mod2[p.adjust(go_mod2$over_represented_pvalue, method = "fdr") < 0.05,]
enrichGO_mod2$term <- factor(enrichGO_mod2$term, 
                             levels = unique(enrichGO_mod2$term[order(enrichGO_mod2$numDEInCat, decreasing = TRUE)]))

gg.go_mod2 <- ggplot(enrichGO_mod2) +
  geom_bar(aes(x = term, y = numDEInCat), stat = "identity", position = "dodge", fill = "#9acd32") + 
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_mod2

#
#### BIOLOGICAL ENRICHMENT - MODULE 3 ####

## Module 3
DEG_mod3 <- ifelse(names(count.bias) %in% subset(ko.sgm.node_df, membership == 5)$names, 1, 0)
names(DEG_mod3) <- names(count.bias)

pwf_mod3 <- nullp(DEgenes = DEG_mod3, bias.data = count.bias[names(DEG_mod3)])
go_mod3 <- goseq(pwf_mod3, gene2cat = go_df, test.cats = c("GO:BP"))
go_mod3 <- subset(go_mod3, ontology == "BP" & numDEInCat)

enrichGO_mod3 <- go_mod3[p.adjust(go_mod3$over_represented_pvalue, method = "fdr") < 0.05,]
enrichGO_mod3$term <- factor(enrichGO_mod3$term, 
                             levels = unique(enrichGO_mod3$term[order(enrichGO_mod3$numDEInCat, decreasing = TRUE)]))

gg.go_mod3 <- ggplot(enrichGO_mod3) +
  geom_bar(aes(x = term, y = numDEInCat), stat = "identity", position = "dodge", fill = "#910a25") + 
  scale_y_continuous(breaks = seq(2,10,4)) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_mod3

#
#### BIOLOGICAL ENRICHMENT - MODULE 4 ####

## Module 4
DEG_mod4 <- ifelse(names(count.bias) %in% subset(ko.sgm.node_df, membership == 3)$names, 1, 0)
names(DEG_mod4) <- names(count.bias)

pwf_mod4 <- nullp(DEgenes = DEG_mod4, bias.data = count.bias[names(DEG_mod4)])
go_mod4 <- goseq(pwf_mod4, gene2cat = go_df, test.cats = c("GO:BP"))
go_mod4 <- subset(go_mod4, ontology == "BP")

enrichGO_mod4 <- go_mod4[p.adjust(go_mod4$over_represented_pvalue, method = "fdr") < 0.05,]
enrichGO_mod4$term <- factor(enrichGO_mod4$term, 
                             levels = unique(enrichGO_mod4$term[order(enrichGO_mod4$numDEInCat, decreasing = TRUE)]))

gg.go_mod4 <- ggplot(enrichGO_mod4) +
  geom_bar(aes(x = term, y = numDEInCat), stat = "identity", position = "dodge", fill = "#ff8247") + 
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_mod4

#
#### EXPORT FIGURE S8 ####

gg.figureS8 <- plot_grid(plot_grid(gg.go_mod1, gg.go_mod3, ncol = 1, rel_heights = c(1,1), labels = c("A", "C"), label_size = 18), 
                         plot_grid(gg.go_mod2, gg.go_mod4, ncol = 1, rel_heights = c(0.5,1.5), labels = c("B", "D"), label_size = 18), ncol = 2)
gg.figureS8

ggsave("Figures/Figure_S8.png", gg.figureS8, bg = "white", width = 18, height = 13)


#