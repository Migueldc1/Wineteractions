# Author: M, de Celis Rodriguez
# Date: 18/10/2023
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
sample_sc <- readRDS("Data/Meta-transcriptomics/sample_sgm.rds")
row.names(sample_sc) <- sample_sc$Sample_ID


## DESeq2 OBJECTS
dds_sc <- readRDS("Data/Meta-transcriptomics/dds_sc.rds")
sample_sc <- sample_sc[colnames(dds_sc),]


## DEO TABLES
ress_sc.18C <- read.table("Data/Meta-transcriptomics/ress_sc.18C.txt", sep = "\t", header = TRUE)
ress_sc.NH4 <- read.table("Data/Meta-transcriptomics/ress_sc.NH4.txt", sep = "\t", header = TRUE)
ress_sc.SO2 <- read.table("Data/Meta-transcriptomics/ress_sc.SO2.txt", sep = "\t", header = TRUE)


## BIOLOGICAL ENRICHMENT
enrichGO_sc.18C <- read.table("Data/Meta-transcriptomics/enrichGO_sc.18C.txt", header = TRUE, sep = "\t")
enrichGO_sc.NH4 <- read.table("Data/Meta-transcriptomics/enrichGO_sc.NH4.txt", header = TRUE, sep = "\t")


## COLORS
col_cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")

#
################################################################################ FIGURE 4 ####
#### PERMANOVAs ####

vst_sc <- vst(dds_sc, blind = TRUE)
norm_smg.dis <- as.matrix(dists <- dist(t(assay(vst_sc))))

adonis2(norm_smg.dis ~ Condition, sample_sc[row.names(norm_smg.dis),])

#
#### DIFFERENTIAL EXPRESSION - VENN DIAGRAM ####

venn.df_sc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                     list(ress_sc.18C[,c(1,10)], ress_sc.NH4[,c(1,10)], ress_sc.SO2[,c(1,10)]))
colnames(venn.df_sc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sc[is.na(venn.df_sc)] <- 0

venn.plot_sc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 0)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,2] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,3] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,4] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,4] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,3] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,2] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 3)))

venn.plot_sc <- cbind.data.frame(venn.plot_sc, 
                                 x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                 y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col_cond[-1]) +
  labs(fill = NULL)

gg.venn_sc

#
#### DIFFERENTIAL EXPRESSION - HISTOGRAM ####

hist_sc <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                 cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                 cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_sc, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col_cond[-1]) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 14, color  = "black"),
        axis.title.x = element_text(size = 16, color  = "black"),
        axis.title.y = element_text(size = 16, color  = "black"),
        legend.text = element_text(size = 15, color  = "black"),
        legend.title = element_text(size = 16, color  = "black"),
        axis.text.x = element_text(size = 14, color  = "black")) 

gg.hist_sc

#
#### BIOLOGICAL ENRICHMENT vs CONTROL ####

go_sc <- rbind.data.frame(cbind(enrichGO_sc.18C, Comparison = "18C"),
                          cbind(enrichGO_sc.NH4, Comparison = "NH4"))

go_sc$Comparison <- factor(go_sc$Comparison, levels = c("18C", "NH4"))
go_sc$comp.cat <- ifelse(duplicated(go_sc$term) | duplicated(go_sc$term, fromLast = TRUE), "Both", go_sc$Comparison)
go_sc$comp.cat <- factor(go_sc$comp.cat, levels = c("1", "Both", "2"))

go_sc$term <- factor(go_sc$term, 
                     levels = unique(go_sc$term[order(go_sc$comp.cat, go_sc$numDEInCat, decreasing = TRUE)]))

gg.go_sc <- ggplot(go_sc) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = col_cond[2:3]) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_sc

#
#### EXPORT FIGURE 4 ####

gg.figure4 <- plot_grid(plot_grid(gg.venn_sc, gg.hist_sc, labels = c("A", "B"), label_size = 18),
                        gg.go_sc, rel_heights = c(1, 1.3), ncol = 1, labels = c("", "C"), label_size = 18)
gg.figure4

ggsave("Figures/Figure_4.png", gg.figure4, bg = "white", width = 9, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S7 ####
#### METABOLITE PRODUCTION - CONDITION ####

metabolite_sc <- melt(sample_sc[,c(4,6:23)])

metabolite_sc$variable <- factor(metabolite_sc$variable,
                                 levels = c("Glucose", "Fructose", "Ethanol", "Glycerol",  "pH", 
                                            "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                            "Citric_acid", "Succinic_acid", "EEFA", "SCFA", "MCFA", "Fusel.alcohols", 
                                            "Fusel.alcohol.acetates", "Ethyl.acetate"))

metabolite_sc <- metabolite_sc[complete.cases(metabolite_sc),]

gg.meta_sc <- ggplot(metabolite_sc) +
  geom_point(aes(x = Condition, y = value, color = Condition), size = 3, show.legend = FALSE) +
  scale_color_manual(values = col_cond) +
  facet_wrap(~ variable, ncol = 4, scales = "free_y") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 16, color = "black"))


## Set Strip background colors
gg.meta_sc <- ggplot_gtable(ggplot_build(gg.meta_sc))

strip_both <- which(grepl("strip-", gg.meta_sc$layout$name))
fills <- alpha(c("#82107a", "#82107a", "#82107a", "#82107a", "#30b3bf", "#30b3bf", "#82107a", "#82107a", "#30b3bf", 
                 "#30b3bf", "#30b3bf", "#30b3bf", "#81c236", "#81c236", "#81c236", "#bbcc06", "#bbcc06", "#30b3bf"), 
               alpha = 0.75)

for (i in 1:length(strip_both)) {
  
  gg.meta_sc$grobs[[strip_both[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fills[i]
  
}

grid::grid.draw(gg.meta_sc)

#
#### EXPORT FIGURE S7 ####

gg.figureS7 <- gg.meta_sc

ggsave("Figures/Figure_S7.png", gg.figureS7, bg = "white", width = 12.6, height = 6)

#