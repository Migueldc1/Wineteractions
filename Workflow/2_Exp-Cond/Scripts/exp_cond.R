# Author: M, de Celis Rodriguez
# Date: 17/07/2022
# Project: Wineteractions

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/2_Exp-Cond")

#
#### LIBRARIES ####
library(vegan)
library(ggplot2)
library(reshape2)
library(phyloseq)

load("exp_cond.RData")

#
#### FUNCTIONS ####
set_unid <- function(tax_df) {
  
  for (i in 1:ncol(tax_df)) {
    if (i != 7) {
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  paste("Unidentified", tax_df[,i-1], sep = " ")),
                           tax_df[,i])
    } else{
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  "sp."),
                           tax_df[,i])
      
    }
    
  }
  
  return(tax_df)
  
}

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("../0_ITS-GM/Inputs/sample_GM.txt", sep = "\t", header = TRUE)
sample_df <- subset(sample_df, Stage == "2_final")

sample_df <- sample_df[,-c(5,13,16:17,20:24)]
colnames(sample_df)[c(5,13:16)] <- c("Time.2", "Ammonia.2", "Sugars.2", "Malic_acid.2", "PAN.2" )

sample_df.f <- read.table("Inputs/sample_GM.end.txt", sep = "\t", header = TRUE)




row.names(sample_df) <- sample_df$Seq_ID




sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", 
                                                        "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha",
                                                        "Madrid", "Rioja"))






## Community data
asv_GM <- readRDS("../0_ITS-GM/Outputs/ASV_GM.rds")
asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

asv_GM.f <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM.f <- asv_GM.f[, colSums(asv_GM.f != 0) > 0]

tax_GM <- readRDS("../0_ITS-GM/Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM.f <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]

tax_GM.f.un <- set_unid(tax_GM.f)
tax_GM[is.na(tax_GM)] <- "Unidentified"
tax_GM.f[is.na(tax_GM.f)] <- "Unidentified"

#
#### TAXONOMIC EXPLORATION ####
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1:7)], by = "Seq_ID")

## Genus

asv.t_plot <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Sample_ID, asv.t_GM.p$Genus, asv.t_GM.p$Stage), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Stage", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus, asv.t_plot$Stage), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Stage", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other")]
orderG <- append(orderG, c("Other"))

asv.t_plot <- merge(asv.t_plot, sample_df[,2:7], by = c("Sample_ID", "Stage"))
asv.t_plot$Origin <- factor(asv.t_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

asv.t_plot$Sample_name <- paste(asv.t_plot$Farming, asv.t_plot$Condition, sep = "-")
asv.t_plot$Sample_name <- factor(asv.t_plot$Sample_name, levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                                                    "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))

ggplot(asv.t_plot, 
       aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ffff33", "#bc80bd", "#cab2d6", "#8da0cb",
                                               "#b15928", "#1b9e77", "#6a3d9a")) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        strip.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 1)) + facet_wrap(~Origin, nrow = 1)


#
#### SET DOMINANT sp. ####

dom_sp <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID), max)
colnames(dom_sp) <- c("Sample_ID", "value")
dom_sp <- merge(dom_sp, asv.t_plot[,c(1,3,4)], by = c("Sample_ID", "value"))

dom_sp$Dominant <- ifelse(dom_sp$Genus == "Saccharomyces" & dom_sp$value >= 0.5, "Saccharomyces", 
                          ifelse(dom_sp$Genus == "Lachancea" & dom_sp$value >= 0.5, "Lachancea", "Other"))

dom_sp$Dominant <- factor(dom_sp$Dominant, levels = c("Saccharomyces", "Lachancea", "Other"))

sample_df <- merge(sample_df, dom_sp[,c(1,4)], by = "Sample_ID")

#
#### FINAL WINES ####
cor.test(sample_df$Sugars_y15, sample_df$Glucose.Fructose)
cor.test(sample_df$pH_must, sample_df$pH)
cor.test(sample_df$Glycerol_y15, sample_df$Glycerol)
cor.test(sample_df$Malic_acid_y15, sample_df$Malic_acid)
cor.test(sample_df$Lactic_acid_y15, sample_df$Lactic_acid)
cor.test(sample_df$Acetic_acid_y15, sample_df$Acetic_acid)

#Teniendo el dato repetido en y15 e IR, me quedo con IR ¿?

wine_df <- sample_df[,c(1:2,16:20,3:15,24:25,30)]

## Ordination (PCA)

wine_pca <- prcomp(wine_df[,c(8:22)], scale = TRUE)

wine_pca.plot <- wine_pca$x
wine_pca.plot <- merge(wine_df, wine_pca.plot, by = "row.names")

ggplot(wine_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Condition, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  xlab(paste("PC1: ", round(((wine_pca$sdev)^2 / sum((wine_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((wine_pca$sdev)^2 / sum((wine_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 24, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"))

wine_pca

adonis2(as.matrix(dist(wine_df[,c(8:22)])) ~ Origin+Dominant, 
        wine_df[,c(3,4,5,7,23)], 
        permutations = 1000)


boxplot(wine_df$Ammonia_y15~wine_df$Dominant)

t.test(subset(wine_df, Dominant == "Saccharomyces" & Condition != "NH4")$PAN_y15,
       subset(wine_df, Dominant == "Lachancea" & Condition != "NH4")$PAN_y15)





# REMOVING NH4 SAMPLES

adonis2(as.matrix(dist(wine_df[wine_df$Condition != "NH4",c(8:22)])) ~ Farming, 
        wine_df[wine_df$Condition != "NH4", c(3,4,5,7,23)], 
        permutations = 1000)

#

## Metabolitos/Dominance

wine_df.p <- melt(sample_df[,c(1:2,16:20,3:15,24:25,30,21:22)])

ggplot(wine_df.p) + 
  geom_boxplot(aes(x = Dominant, y = value, color = Dominant), size = 1) +
  scale_color_manual(values = c("#1b9e77", "#8da0cb", "#6a3d9a")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x = element_text(size = 18, color = "black")) +
  facet_wrap(~variable, scales = "free_y", nrow = 3)

#
#### Sacch vs. Lachan t0 ####

sample_df.0 <- read.table("../0_ITS-GM/Inputs/sample_df.txt", sep = "\t", header = TRUE)
row.names(sample_df.0) <- sample_df.0$Seq_ID
sample_df.0$Condition <- factor(sample_df.0$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df.0$Origin <- factor(sample_df.0$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", 
                                                        "R3B", "R3C"))
sample_df.0$Region <- factor(sample_df.0$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha",
                                                        "Madrid", "Rioja"))

sample_df.0 <- subset(sample_df.0, Stage == "0_initial")
sample_df.0 <- merge(sample_df.0, dom_sp[,c(1,4)], by = "Sample_ID", all = TRUE)

sample_df.0 <- subset(sample_df.0, Dominant != "Other")
sample_df.0$Dominant <- factor(sample_df.0$Dominant, levels = c("Saccharomyces", "Lachancea"))

asv.t_GM.0 <- asv.t_GM[ ,colnames(asv.t_GM) %in% sample_df.0$Seq_ID]

## Differential abundance analysis

phy_GM.0 <- phyloseq(otu_table(asv.t_GM.0, taxa_are_rows = TRUE), tax_table(tax_GM), 
                     sample_data(sample_df.0))

phy_GM.0 <- subset_samples(phy_GM.0, !is.na(Dominant))
phy_GM.0 <- subset_samples(phy_GM.0, Dominant != "Other")

phy_GM.0.genus <- tax_glom(phy_GM.0, taxrank = "Genus")


asv_GM.0_genus <- phy_GM.0.genus@otu_table@.Data
asv_GM.0_genus <- melt(asv_GM.0_genus)
colnames(asv_GM.0_genus) <- c("Id", "Seq_ID", "value")
asv_GM.0_genus <- merge(asv_GM.0_genus, tax_GM[,c(6,8)], by = "Id")
asv_GM.0_genus <- merge(asv_GM.0_genus, sample_df.0[,c(2,26)], by = "Seq_ID")

t.test(subset(asv_GM.0_genus, Dominant == "Saccharomyces" & Genus == "Lachancea")$value,
       subset(asv_GM.0_genus, Dominant == "Lachancea" & Genus == "Lachancea")$value)


## Dominant ~ Global
library(DESeq2)

dds_GM_man.tot <- phyloseq_to_deseq2(phy_GM.0.genus, ~ Dominant)

keep <- rowSums(counts(dds_GM_man.tot)) >= 100
dds_GM_man.tot <- dds_GM_man.tot[keep,]

dds_GM_man.tot <- DESeq(dds_GM_man.tot, test = "Wald", fitType = "parametric")

res_GM_man.tot <- results(dds_GM_man.tot, cooksCutoff = FALSE)

sigtab_GM_man.tot <- as.data.frame(res_GM_man.tot[which(res_GM_man.tot$padj < 0.05), ])
sigtab_GM_man.tot <- merge(sigtab_GM_man.tot, tax_GM.un[,6], by = "row.names")
colnames(sigtab_GM_man.tot)[c(1,8)] <- c("Id", "Genus")

x <- tapply(sigtab_GM_man.tot$log2FoldChange, sigtab_GM_man.tot$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_GM_man.tot$Genus = factor(as.character(sigtab_GM_man.tot$Genus), levels=names(x))

ggplot(sigtab_GM_man.tot, aes(x = Genus, y = log2FoldChange)) + 
  geom_point(size = 4) + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"))


summary(res_GM_man.tot)

## Pheatmap

matrix <- t(as.matrix(data.frame(otu_table(phy_GM.0.genus))))
rownames(matrix) <- as.character(tax_table(phy_GM.0.genus)[, "Genus"])
matrix <- matrix[row.names(matrix) %in% sigtab_GM_man.tot$Genus,]

metadata_sub <- data.frame(sample_data(phy_GM.0.genus))

# Define the annotation color for columns and rows
annotation_col <- data.frame(Farming = as.factor(metadata_sub$Farming), 
                             `Dominant` = as.factor(metadata_sub$Dominant), 
                             check.names = FALSE)

rownames(annotation_col) <- rownames(metadata_sub)

# ann_color should be named vectors
ann_colors <- list(Farming = c(`CONV` = "red", `ECO` = "blue"),
                   `Dominant` = c(Saccharomyces = "#1b9e77", Lachancea = "#8da0cb"))


matrix <- matrix[,c(row.names(sample_df.0[sample_df.0$Dominant == "Saccharomyces",]),
                    row.names(sample_df.0[sample_df.0$Dominant == "Lachancea",]))]

pheatmap::pheatmap(matrix, scale = "row", cluster_cols = FALSE,
                         annotation_col = annotation_col, 
                         annotation_colors = ann_colors)
#
#### EFFECT OF EXPERIMENTAL CONDITION ON ANYTHING ####

#### SAVE OUTPUTs ####
save.image("exp_cond.RData")



#
