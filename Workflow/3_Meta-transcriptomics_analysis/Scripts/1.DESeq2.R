# Author: M, de Celis Rodriguez
# Date: 22/12/2022
# Project: Wineteractions - Metatranscriptomic RNAseq Analysis

library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggforce)
library(clusterProfiler)
library(reshape2)
library(goseq)

rm(list = ls())

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions")

# Load Environment
#load("Workflow/3_RNAseq-meta/Outputs/RData/RNAseq-meta.RData")
#save.image("Workflow/3_RNAseq-meta/Outputs/RData/RNAseq-meta.RData")

#
#### CUSTOM FUNCTION ####

## Load Count data + Emapper Annotations to get KEGG Ortholog counts
read.emapper <- function(path.count, path.orth) {
  
  count.list <- NULL
  
  for (smpl in list.files(path.count)) {
    sample <- gsub(".txt", "", smpl)
    
    count_df <- read.table(paste(path.count, smpl, sep = "/"), header = TRUE)[,c(1,7)]
    colnames(count_df) <- c("seed_ortholog", sample)
    
    orth_file <- list.files(paste(path.orth, sample, sep = "/"), full.names = TRUE)
    orth_df <- read.delim(orth_file, sep = "\t", skip = 4)
    orth_df <- unique(orth_df[1:(nrow(orth_df)-3), c(2,12)])
    
    count.deff_df <- merge(count_df, orth_df, by = "seed_ortholog")
    count.deff_df <- subset(count.deff_df, KEGG_ko != "-")
    
    count.deff_df <- aggregate(count.deff_df[,2], list(count.deff_df$KEGG_ko), sum)
    colnames(count.deff_df) <- c("KEGG_ko", sample)
    
    count.list[[sample]] <- count.deff_df
    
  }
  
  ko_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), count.list)
  ko_df[is.na(ko_df)] <- 0
  
  return(ko_df)
  
}

## Load Gene Ontology Data
read.GOdata <- function(file.ext, path.annot, cols) {
  
  kegg.go_df <- NULL
  
  for (smpl in list.files(path.annot)) {
    
    #Open Annotation files and select KO and GO columns
    annot_df <- read.delim(list.files(paste(path.annot, smpl, sep = "/"), full.names = TRUE), sep = "\t", skip = 4)[,cols]
    annot_df <- unique(annot_df[annot_df$KEGG_ko != "" & annot_df$KEGG_ko != "-", ])
    
    #Bind all the annotation files in a single table, removing duplicated info
    kegg.go_df <- rbind.data.frame(kegg.go_df, annot_df)
    kegg.go_df <- unique(kegg.go_df)
    
  }
  
  #Format the table into a useful form
  go_df <- NULL
  for (n in 1:nrow(annot_df)) {
    
    go_df <- rbind.data.frame(go_df, cbind(gene_id = annot_df$KEGG_ko[n],
                                           category = unlist(strsplit(as.character(annot_df$GOs[n]),",", fixed = TRUE))))
    
    go_df <- unique(go_df)
    
  }
  
  return(go_df)
  
}

#
#### LOAD DATA ####

ko_df <- read.emapper(path.count = "Data/RNAseq/Meta-transcriptomics/Count/eggnog", 
                      path.orth = "Data/RNAseq/Meta-transcriptomics/Annotation/eggnog")

ko_df.f <- ko_df[!grepl(",", ko_df$KEGG_ko),]

## ANNOTATION
kegg_df <- read.table("Data/RNAseq/Meta-transcriptomics/Annotation/KEGG_names.txt", 
                      sep = "\t", header = TRUE, quote = "")

go_df <- read.GOdata(".emapper.annotations", "Data/RNAseq/Meta-transcriptomics/Annotation/eggnog", c("KEGG_ko", "GOs"))

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Genus <- factor(sample_df$Genus, levels = c("Citeromyces", "Hanseniaspora", "Kluyveromyces", "Lachancea", 
                                                      "Saccharomyces", "Other"))

## COLORS
col.cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col.gen <- c("#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77", "gray70")

#
#### DIFFERENTIAL ANALYSIS - GLOBAL ####
dds_tot <- DESeqDataSetFromMatrix(countData = ko_df.f, 
                                  colData = sample_df, 
                                  design = ~ Condition + Genus + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_tot)) >= 10
dds_tot <- dds_tot[keep,]

dds_tot <- DESeq(dds_tot, fitType = "local")
resultsNames(dds_tot)

## PCA
vst_tot <- vst(dds_tot, blind = TRUE)

pcaData_tot <- plotPCA(vst_tot, intgroup = "Condition", returnData = TRUE)
percentVar_tot <- round(100 * attr(pcaData_tot, "percentVar"), 2)
pcaData_tot <- merge(pcaData_tot, sample_df[,-4], by.x = "name", by.y = "Sample_ID")

gg.pca_cond <- ggplot(pcaData_tot, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_tot[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_tot[2],"% variance")) + 
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")

# Samples are separated by dominant species
gg.pca_gen <- ggplot(pcaData_tot, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_tot[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_tot[2],"% variance")) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values = col.gen) +
  labs(color = "Condition")

gg.pca_gen

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/total_PCA.png", 
       plot_grid(gg.pca_cond, gg.pca_gen, align = "v", labels = c("a", "b"), label_size = 20), 
       width = 15.2, height = 5.4, dpi = 300, bg = "white")

#
################################################### ACABAR "DIFFERENTIAL ANALYSIS - GLOBAL" AQUÍ ####
## DIFFERENTIAL EXPRESSION
# Saccharomyces v Lachancea
res_lt.sc <- results(dds_tot, contrast = c("Genus", "Lachancea", "Saccharomyces"), alpha = 0.05)
summary(res_lt.sc)

ress_lt.sc <- as.data.frame(res_lt.sc)
ress_lt.sc$KEGG_ko <- row.names(ress_lt.sc)

ress_lt.sc <- merge(ress_lt.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.sc <- ress_lt.sc[order(ress_lt.sc$padj),]
ress_lt.sc$DEO <- ifelse(abs(ress_lt.sc$log2FoldChange) > 1 & ress_lt.sc$padj <= 0.05, 1, 0)

# Saccharomyces v Hanseniaspora
res_hs.sc <- results(dds_tot, contrast = c("Genus", "Hanseniaspora", "Saccharomyces"), alpha = 0.05)
summary(res_hs.sc)

ress_hs.sc <- as.data.frame(res_hs.sc)
ress_hs.sc$KEGG_ko <- row.names(ress_hs.sc)

ress_hs.sc <- merge(ress_hs.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.sc <- ress_hs.sc[order(ress_hs.sc$padj),]
ress_hs.sc$DEO <- ifelse(abs(ress_hs.sc$log2FoldChange) > 1 & ress_hs.sc$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_gen <- merge(ress_lt.sc[,c(1, 10)], ress_hs.sc[,c(1, 10)], by = "KEGG_ko")
colnames(venn.df_gen) <- c("KEGG_ko", "Lt.Sc", "Hs.Sc")
venn.df_gen[is.na(venn.df_gen)] <- 0

venn.plot_gen <- rbind.data.frame(cbind(Lt.Sc = 0, Hs.Sc = 0, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 0)),
                                  cbind(Lt.Sc = 0, Hs.Sc = 1, 
                                        Counts = sum(venn.df_gen[,2] == 0 & venn.df_gen[,3] == 1)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 0, 
                                        Counts = sum(venn.df_gen[,2] == 1 & venn.df_gen[,3] == 0)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 1, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 2)))

venn.plot_gen <- cbind.data.frame(venn.plot_gen, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out2 <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("Lachancea", "Hanseniaspora"))
venn.out2$labels <- factor(venn.out2$labels, levels = c("Lachancea", "Hanseniaspora"))

gg.venn_gen <- ggplot() +
  geom_circle(data = venn.out2, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_gen, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = NULL)

gg.venn_gen

# Accumulated LFC - Histogram
hist_gen <- rbind(cbind(subset(ress_lt.sc, DEO == 1), Genus = "Lachancea"),
                  cbind(subset(ress_hs.sc, DEO == 1), Genus = "Hanseniaspora"))

hist_gen$Genus <- factor(hist_gen$Genus, levels = c("Lachancea", "Hanseniaspora"))

gg.hist_gen <- ggplot(hist_gen, aes(x = abs(log2FoldChange), fill = Genus)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

gg.hist_gen

#
ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/total_summary.png", 
       plot_grid(gg.venn_gen, gg.hist_gen, nrow = 1, rel_widths = c(1, 2)), 
       width = 12, height = 8, dpi = 300, bg = "white")

##### PRUEBAS

plotCounts(dds_tot, "ko:K00474", intgroup = "Genus")

#
## BIOLOGICAL ENRICHMENT - KEGG
DEO_lt.sc <- gsub("ko:", "", subset(ress_lt.sc, DEO == 1)$KEGG_ko)
enrich_lt.sc <- enrichKEGG(DEO_lt.sc, organism = "ko", keyType = "kegg")
kk_lt.sc <- data.frame(enrich_lt.sc, Genus = "Lachancea")


DEO_hs.sc <- gsub("ko:", "", subset(ress_hs.sc, DEO == 1)$KEGG_ko)
enrich_hs.sc <- enrichKEGG(DEO_hs.sc, organism = "ko", keyType = "kegg")
kk_hs.sc <- data.frame(enrich_hs.sc, Genus = "Hanseniaspora")

kk_gen <- merge(kk_lt.sc[,c(2,6)], kk_hs.sc[,c(2,6)], by = "Description", all = TRUE)
kk_gen[is.na(kk_gen)] <- 1
colnames(kk_gen) <- c("Description", "Lachancea", "Hanseniaspora")
kk_gen <- melt(kk_gen)

kk_gen$variable <- factor(kk_gen$variable, levels = c("Lachancea", "Hanseniaspora"))
kk_gen <- kk_gen[order(kk_gen$value),]

to.rm <- c("Huntington disease", "Oocyte meiosis", "Amyotrophic lateral sclerosis",
           "Prion disease", "Progesterone-mediated oocyte maturation", "Fanconi anemia pathway",
           "Carbon fixation in photosynthetic organisms", "Coronavirus disease - COVID-19",
           "Longevity regulating pathway - multiple species", 
           "Human T-cell leukemia virus 1 infection", 
           "Tropane, piperidine and pyridine alkaloid biosynthesis", "Terpenoid backbone biosynthesis",
           "Pathogenic Escherichia coli infection", "Synaptic vesicle cycle", "Parkinson disease")

kk_gen <- kk_gen[!kk_gen$Description %in% to.rm,]
kk_gen 

ggplot(kk_gen, aes(x = reorder(Description, log(value)), y = -log(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = "Genus") + xlab("KEGG pathways") + ylab("-log adjusted p value") +
  coord_flip()

## BIOLOGICAL ENRICHMENT - Gene Ontology BP
count.bias <- rowSums(ko_df.f[,-1])
names(count.bias) <- ko_df.f[,1]

DEG_lt.sc <- as.integer(ress_lt.sc$padj < 0.05 & abs(ress_lt.sc$log2FoldChange) >= 1)
DEG_lt.sc[is.na(DEG_lt.sc)] <- 0
names(DEG_lt.sc) <- ress_lt.sc$KEGG_ko

pwf_lt.sc <- nullp(DEgenes = DEG_lt.sc, bias.data = count.bias[names(DEG_lt.sc)], )
go_lt.sc <- goseq(pwf_lt.sc, gene2cat = go_df)
go_lt.sc <- subset(go_lt.sc, ontology == "BP" & numDEInCat > 0)

enrichedGO_lt.sc <- go_lt.sc[p.adjust(go_lt.sc$over_represented_pvalue, method = "fdr") < 0.05]



DEG_hs.sc <- as.integer(ress_hs.sc$padj < 0.05 & abs(ress_hs.sc$log2FoldChange) >= 1)
names(DEG_hs.sc) <- ress_hs.sc$KEGG_ko

pwf_hs.sc <- nullp(DEgenes = DEG_hs.sc, bias.data = count.bias[names(DEG_hs.sc)], )
go_hs.sc <- goseq(pwf_hs.sc, gene2cat = go_df)
go_hs.sc <- subset(go_hs.sc, ontology == "BP" & numDEInCat > 0)


enrichedGO_lt.sc <- go_lt.sc[p.adjust(go_lt.sc$over_represented_pvalue, method = "fdr") < 0.05]


go_lt.sc$over_represented_pvalue[1:10]
p.adjust(go_lt.sc$over_represented_pvalue[1:10], method = "fdr")



go_lt.sc <- subset(go_lt.sc, ontology == "BP" & over_represented_pvalue <= 0.05)
enrichedGO_lt.sc <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                              + method="BH") < 0.05]


## PLOT
GO_plot.BP <- rbind.data.frame(cbind(GO_Lt.BP, Comparison = "Lachancea"),
                               cbind(GO_Hs.BP, Comparison = "Hanseniaspora"))

GO_plot.BP$Comparison <- factor(GO_plot.BP$Comparison, levels = c("Lachancea", "Hanseniaspora"))

ggplot(GO_plot.BP) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  coord_flip() +
  theme_bw()+
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("Gene count")  + xlab("")




#
#
#### DIFFERENTIAL ANALYSIS - GLOBAL.Genus ####
# KEEP SAMPLES DOMINATED BY ONE OF SAID SPECIES
sample_gen <- sample_df[sample_df$Genus %in% c("Saccharomyces", "Lachancea", "Hanseniaspora"),]
ko_df.gen <- ko_df.f[,c("KEGG_ko", sample_gen$Sample_ID)]

dds_gen <- DESeqDataSetFromMatrix(countData = ko_df.gen, 
                                  colData = sample_gen, 
                                  design = ~ Condition + Origin + Genus, tidy = TRUE)

dds_gen$Genus <- relevel(dds_gen$Genus, "Saccharomyces")

keep <- rowSums(counts(dds_gen) >= 10) >= 5
dds_gen <- dds_gen[keep,]

dds_gen <- DESeq(dds_gen, fitType = "local", betaPrior = FALSE)
resultsNames(dds_gen)

## PCA
vst_gen <- vst(dds_gen, blind = TRUE)

pcaData_gen <- plotPCA(vst_gen, intgroup = "Condition", returnData = TRUE)
percentVar_gen <- round(100 * attr(pcaData_gen, "percentVar"), 2)
pcaData_gen <- merge(pcaData_gen, sample_df[,-4], by.x = "name", by.y = "Sample_ID")

# Samples are separated by dominant species
gg.pca_gen <- ggplot(pcaData_gen, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_gen[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_gen[2],"% variance")) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.gen[c(2,4,5)]) +
  labs(color = "Condition")

gg.pca_gen

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/genus_PCA.png", gg.pca_gen, 
       width = 7.6, height = 5.4, dpi = 300, bg = "white")

#
## DIFFERENTIAL EXPRESSION
# Saccharomyces v Lachancea
res_lt.sc <- results(dds_gen, contrast = c("Genus", "Lachancea", "Saccharomyces"), alpha = 0.05)
res_lt.sc <- lfcShrink(dds_gen, res = res_lt.sc, type = "ashr")
summary(res_lt.sc)

ress_lt.sc <- as.data.frame(res_lt.sc)
ress_lt.sc$KEGG_ko <- row.names(ress_lt.sc)

ress_lt.sc <- merge(ress_lt.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.sc <- ress_lt.sc[order(ress_lt.sc$padj),]
ress_lt.sc$DEO <- ifelse(abs(ress_lt.sc$log2FoldChange) > 1 & ress_lt.sc$padj <= 0.05, 1, 0)

# Saccharomyces v Hanseniaspora
res_hs.sc <- results(dds_gen, contrast = c("Genus", "Hanseniaspora", "Saccharomyces"), alpha = 0.05)
res_hs.sc <- lfcShrink(dds_gen, res = res_hs.sc, type = "ashr")
summary(res_hs.sc)

ress_hs.sc <- as.data.frame(res_hs.sc)
ress_hs.sc$KEGG_ko <- row.names(ress_hs.sc)

ress_hs.sc <- merge(ress_hs.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.sc <- ress_hs.sc[order(ress_hs.sc$padj),]
ress_hs.sc$DEO <- ifelse(abs(ress_hs.sc$log2FoldChange) > 1 & ress_hs.sc$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_gen <- merge(ress_lt.sc[,c(1, 9)], ress_hs.sc[,c(1, 9)], by = "KEGG_ko")
colnames(venn.df_gen) <- c("KEGG_ko", "Lt.Sc", "Hs.Sc")

venn.plot_gen <- rbind.data.frame(cbind(Lt.Sc = 0, Hs.Sc = 0, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 0)),
                                  cbind(Lt.Sc = 0, Hs.Sc = 1, 
                                        Counts = sum(venn.df_gen[,2] == 0 & venn.df_gen[,3] == 1)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 0, 
                                        Counts = sum(venn.df_gen[,2] == 1 & venn.df_gen[,3] == 0)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 1, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 2)))

venn.plot_gen <- cbind.data.frame(venn.plot_gen, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out2 <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("Lachancea", "Hanseniaspora"))
venn.out2$labels <- factor(venn.out2$labels, levels = c("Lachancea", "Hanseniaspora"))

gg.venn_gen <- ggplot() +
  geom_circle(data = venn.out2, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_gen, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 15, color  = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = NULL)

gg.venn_gen

# Accumulated LFC - Histogram
hist_gen <- rbind(cbind(subset(ress_lt.sc, DEO == 1), Genus = "Lachancea"),
                  cbind(subset(ress_hs.sc, DEO == 1), Genus = "Hanseniaspora"))

hist_gen$Genus <- factor(hist_gen$Genus, levels = c("Lachancea", "Hanseniaspora"))

gg.hist_gen <- ggplot(hist_gen, aes(x = abs(log2FoldChange), fill = Genus)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        axis.text.y = element_text(size = 15, color  = "black"),
        axis.title.x = element_text(size = 17, color  = "black"),
        axis.title.y = element_text(size = 17, color  = "black"),
        legend.text = element_text(size = 15, color  = "black"),
        legend.title = element_text(size = 17, color  = "black"),
        axis.text.x = element_text(size = 15, color  = "black"),
        plot.margin = unit(c(0,0,0,0), "cm"))

gg.hist_gen

gg.summary_gen <- plot_grid(gg.venn_gen, gg.hist_gen, nrow = 1, rel_widths = c(2, 3), labels = c("a", "b"), label_size = 20)

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/genus_summary.png", gg.summary_gen, 
       width = 12, height = 7.2, dpi = 300, bg = "white")

#
##### PRUEBAS

plotCounts(dds_gen, "ko:K13076", intgroup = "Genus")

#
## BIOLOGICAL ENRICHMENT - Gene Ontology BP
count.bias <- rowSums(ko_df.gen[,-1])
names(count.bias) <- ko_df.gen[,1]

DEG_lt.sc <- as.integer(ress_lt.sc$padj < 0.05 & abs(ress_lt.sc$log2FoldChange) >= 1)
names(DEG_lt.sc) <- ress_lt.sc$KEGG_ko

pwf_lt.sc <- nullp(DEgenes = DEG_lt.sc, bias.data = count.bias[names(DEG_lt.sc)])
go_lt.sc <- goseq(pwf_lt.sc, gene2cat = go_df, test.cats = c("GO:BP"))
go_lt.sc <- subset(go_lt.sc, ontology == "BP" & numDEInCat > 5)

enrichGO_lt.sc <- go_lt.sc[p.adjust(go_lt.sc$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_hs.sc <- as.integer(ress_hs.sc$padj < 0.05 & abs(ress_hs.sc$log2FoldChange) >= 1)
names(DEG_hs.sc) <- ress_hs.sc$KEGG_ko

pwf_hs.sc <- nullp(DEgenes = DEG_hs.sc, bias.data = jitter(rep(1000, length(DEG_hs.sc))))
go_hs.sc <- goseq(pwf_hs.sc, gene2cat = go_df)
go_hs.sc <- subset(go_hs.sc, ontology == "BP" & numDEInCat > 5)

enrichGO_hs.sc <- go_hs.sc[p.adjust(go_hs.sc$over_represented_pvalue, method = "fdr") < 0.05,]

# Plot
go_gen <- rbind.data.frame(cbind(enrichGO_lt.sc, Comparison = "Lachancea"),
                           cbind(enrichGO_hs.sc, Comparison = "Hanseniaspora"))

go_gen$Comparison <- factor(go_gen$Comparison, levels = c("Lachancea", "Hanseniaspora"))
go_gen$comp.cat <- ifelse(duplicated(go_gen$term) | duplicated(go_gen$term, fromLast = TRUE), "Both", go_gen$Comparison)

go_gen$term <- factor(go_gen$term, 
                      levels = unique(go_gen$term[order(go_gen$comp.cat, go_gen$numDEInCat, decreasing = TRUE)]))

gg.go_gen <- ggplot(go_gen) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "top",
        aspect.ratio = 1.2,
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_gen

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/genus_enrichGO.png", gg.go_gen, 
       width = 10, height = 9, dpi = 300, bg = "white")

#
#### YO BORRARÍA LO DE ABAJO ####
## BIOLOGICAL ENRICHMENT - KEGG
DEO_lt.sc <- gsub("ko:", "", subset(ress_lt.sc, DEO == 1)$KEGG_ko)
enrichKEGG_lt.sc <- enrichKEGG(DEO_lt.sc, organism = "ko", keyType = "kegg")
kk_lt.sc <- data.frame(enrichKEGG_lt.sc, Genus = "Lachancea")


DEO_hs.sc <- gsub("ko:", "", subset(ress_hs.sc, DEO == 1)$KEGG_ko)
enrichKEGG_hs.sc <- enrichKEGG(DEO_hs.sc, organism = "ko", keyType = "kegg", )
kk_hs.sc <- data.frame(enrichKEGG_hs.sc, Genus = "Hanseniaspora")

kk_gen <- merge(kk_lt.sc[,c(2,6)], kk_hs.sc[,c(2,6)], by = "Description", all = TRUE)
kk_gen[is.na(kk_gen)] <- 1
colnames(kk_gen) <- c("Description", "Lachancea", "Hanseniaspora")
kk_gen <- melt(kk_gen)

kk_gen$variable <- factor(kk_gen$variable, levels = c("Lachancea", "Hanseniaspora"))
kk_gen <- kk_gen[order(kk_gen$value),]

to.rm <- c("Huntington disease", "Oocyte meiosis", "Amyotrophic lateral sclerosis",
           "Prion disease", "Progesterone-mediated oocyte maturation", "Fanconi anemia pathway",
           "Carbon fixation in photosynthetic organisms", "Coronavirus disease - COVID-19",
           "Longevity regulating pathway - multiple species", 
           "Human T-cell leukemia virus 1 infection", 
           "Tropane, piperidine and pyridine alkaloid biosynthesis", "Terpenoid backbone biosynthesis",
           "Pathogenic Escherichia coli infection", "Synaptic vesicle cycle", "Parkinson disease",
           "Alzheimer disease", "Autophagy - other", "Cardiac muscle contraction", "Cell cycle", "Diabetic cardiomyopathy",
           "Fc gamma R-mediated phagocytosis", "Insulin signaling pathway", "Non-alcoholic fatty liver disease",
           "One carbon pool by folate", "Pathways of neurodegeneration - multiple diseases", "Platinum drug resistance",
           "Spinocerebellar ataxia", "Thermogenesis", "Viral life cycle - HIV-1", "Choline metabolism in cancer")

kk_gen <- kk_gen[!kk_gen$Description %in% to.rm,]
kk_gen 

ggplot(kk_gen, aes(x = reorder(Description, log(value)), y = -log(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = "Genus") + xlab("KEGG pathways") + ylab("-log adjusted p value") +
  coord_flip()

logFC <- ress_hs.sc$log2FoldChange
names(logFC) <- gsub("ko:", "", ress_hs.sc$KEGG_ko)
logFC <- logFC[names(logFC) %in% DEO_lt.sc]

library(pathview)
pathview(gene.data = logFC, 
         pathway.id = "ko04113", 
         species = "ko",  
         gene.idtype = "kegg")

#
#### DIFFERENTIAL ANALYSIS - SACCHAROMYCES ####
sample_sc <- sample_df[sample_df$Genus == "Saccharomyces",]
sample_sc <- sample_sc[-c(5,6),]

ko_sc <- ko_df.f[,c("KEGG_ko", sample_sc$Sample_ID)]

dds_sc <- DESeqDataSetFromMatrix(countData = ko_sc, 
                                 colData = sample_sc, 
                                 design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc)) >= 10
dds_sc <- dds_sc[keep,]

dds_sc <- DESeq(dds_sc, fitType = "local")
resultsNames(dds_sc)

## PCA
vst_sc <- vst(dds_sc, blind = TRUE)

pcaData_sc <- plotPCA(vst_sc, intgroup = "Condition", returnData = TRUE)
percentVar_sc <- round(100 * attr(pcaData_sc, "percentVar"), 2)
pcaData_sc <- merge(pcaData_sc, sample_sc[,-4], by.x = "name", by.y = "Sample_ID")

gg.pca_sacc <- ggplot(pcaData_sc, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_sc[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sc[2],"% variance")) + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")

gg.pca_sacc

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/sacc_PCA.png", 
       gg.pca_sacc, width = 7.5, height = 6, dpi = 300)

#
## DIFFERENTIAL EXPRESSION

# Low Temperature
res_sc.18C <- results(dds_sc, contrast = c("Condition", "18C", "Control"), alpha = 0.05)
res_sc.18C <- lfcShrink(dds_sc, res = res_sc.18C, contrast = c("Condition", "18C", "Control"), type = "normal")
summary(res_sc.18C)

ress_sc.18C <- as.data.frame(res_sc.18C)
ress_sc.18C$KEGG_ko <- row.names(ress_sc.18C)

ress_sc.18C <- merge(ress_sc.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.18C <- ress_sc.18C[order(ress_sc.18C$padj),]
ress_sc.18C$DEO <- ifelse(abs(ress_sc.18C$log2FoldChange) > 1 & ress_sc.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_sc.NH4 <- results(dds_sc, contrast = c("Condition", "NH4", "Control"), alpha = 0.05)
res_sc.NH4 <- lfcShrink(dds_sc, res = res_sc.NH4, contrast = c("Condition", "NH4", "Control"), type = "normal")
summary(res_sc.NH4)

ress_sc.NH4 <- as.data.frame(res_sc.NH4)
ress_sc.NH4$KEGG_ko <- row.names(ress_sc.NH4)

ress_sc.NH4 <- merge(ress_sc.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.NH4 <- ress_sc.NH4[order(ress_sc.NH4$padj),]
ress_sc.NH4$DEO <- ifelse(abs(ress_sc.NH4$log2FoldChange) > 1 & ress_sc.NH4$padj <= 0.05, 1, 0)

# High Sulfites
res_sc.SO2 <- results(dds_sc, contrast = c("Condition", "SO2", "Control"), alpha = 0.05)
res_sc.SO2 <- lfcShrink(dds_sc, res = res_sc.SO2, contrast = c("Condition", "SO2", "Control"), type = "normal")
summary(res_sc.SO2)

ress_sc.SO2 <- as.data.frame(res_sc.SO2)
ress_sc.SO2$KEGG_ko <- row.names(ress_sc.SO2)

ress_sc.SO2 <- merge(ress_sc.SO2, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.SO2 <- ress_sc.SO2[order(ress_sc.SO2$padj),]
ress_sc.SO2$DEO <- ifelse(abs(ress_sc.SO2$log2FoldChange) > 1 & ress_sc.SO2$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_sacc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                       list(ress_sc.18C[,c(1,10)], ress_sc.NH4[,c(1,10)], ress_sc.SO2[,c(1,10)]))
colnames(venn.df_sacc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sacc[is.na(venn.df_sacc)] <- 0

venn.plot_sacc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 0)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,2] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,3] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,4] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,4] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,3] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,2] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 3)))

venn.plot_sacc <- cbind.data.frame(venn.plot_sacc, 
                                   x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                   y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sacc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = NULL)

gg.venn_sc

# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                   cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col.cond[-1]) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

gg.hist_sc

gg.summary_sc <- plot_grid(gg.venn_sc, gg.hist_sc, nrow = 1, rel_widths = c(2, 3), labels = c("a", "b"), label_size = 20)
gg.summary_sc

#
##### PRUEBAS

plotCounts(dds_sc, "ko:K06641", intgroup = "Condition")

## BIOLOGICAL ENRICHMENT - Gene Ontology BP
count.bias <- rowSums(ko_sc[,-1])
names(count.bias) <- ko_sc[,1]

DEG_sc.18C <- as.integer(ress_sc.18C$padj < 0.05 & abs(ress_sc.18C$log2FoldChange) >= 1)
names(DEG_sc.18C) <- ress_sc.18C$KEGG_ko

pwf_sc.18C <- nullp(DEgenes = DEG_sc.18C, bias.data = count.bias[names(DEG_sc.18C)])
go_sc.18C <- goseq(pwf_sc.18C, gene2cat = go_df, test.cats = c("GO:BP"))
go_sc.18C <- subset(go_sc.18C, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.18C <- go_sc.18C[p.adjust(go_sc.18C$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_sc.NH4 <- as.integer(ress_sc.NH4$padj < 0.05 & abs(ress_sc.NH4$log2FoldChange) >= 1)
names(DEG_sc.NH4) <- ress_sc.NH4$KEGG_ko

pwf_sc.NH4 <- nullp(DEgenes = DEG_sc.NH4, bias.data = count.bias[names(DEG_sc.NH4)])
go_sc.NH4 <- goseq(pwf_sc.NH4, gene2cat = go_df)
go_sc.NH4 <- subset(go_sc.NH4, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.NH4 <- go_sc.NH4[p.adjust(go_sc.NH4$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_sc.SO2 <- as.integer(ress_sc.SO2$padj < 0.05 & abs(ress_sc.SO2$log2FoldChange) >= 1)
names(DEG_sc.SO2) <- ress_sc.SO2$KEGG_ko
DEG_sc.SO2[is.na(DEG_sc.SO2)] <- 0

pwf_sc.SO2 <- nullp(DEgenes = DEG_sc.SO2, bias.data = jitter(rep(1000, length(DEG_sc.SO2))))
go_sc.SO2 <- goseq(pwf_sc.SO2, gene2cat = go_df)
go_sc.SO2 <- subset(go_sc.SO2, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.SO2 <- go_sc.SO2[p.adjust(go_sc.SO2$over_represented_pvalue, method = "fdr") < 0.05,]

# Plot
go_sc <- rbind.data.frame(cbind(enrichGO_sc.18C, Comparison = "18C"),
                          cbind(enrichGO_sc.NH4, Comparison = "NH4"))

go_sc$comp.cat <- ifelse(duplicated(go_sc$term) | duplicated(go_sc$term, fromLast = TRUE), "Both", go_sc$Comparison)

go_sc$term <- factor(go_sc$term, levels = unique(go_sc$term[order(go_sc$comp.cat, go_sc$numDEInCat, decreasing = TRUE)]))

gg.go_sc <- ggplot(go_sc) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#1e74eb", "#ebb249")) +
  coord_flip() +
  theme_bw()+
  theme(aspect.ratio = 0.62,
        legend.position = "bottom",
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_sc

ggsave("Workflow/3_Meta-transcriptomics_analysis/Outputs/Figures/sc_summary.png", 
       plot_grid(gg.summary_sc, gg.go_sc, nrow = 2, labels = c("", "c"), label_size = 20, 
                 rel_widths = c(1,2)), 
       width = 12, height = 12, dpi = 300, bg = "white")


#### LO QUE ESTABA ####
## DIFFERENTIAL EXPRESSION
# Low Temperature

ko_sc.18C <- ko_sc[,c("KEGG_ko", sample_sc[sample_sc$Condition %in% c("Control", "18C"), 1])]
sample_sc.18C <- sample_sc[sample_sc$Condition %in% c("Control", "18C"), ]

dds_sc.18C <- DESeqDataSetFromMatrix(countData = ko_sc.18C, 
                                 colData = sample_sc.18C, 
                                 design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc.18C) >= 10) >= 3
dds_sc.18C <- dds_sc.18C[keep,]

dds_sc.18C <- DESeq(dds_sc.18C, fitType = "local")
res_sc.18C <- lfcShrink(dds_sc, coef = 2, type = "ashr")
summary(res_sc.18C)

ress_sc.18C <- as.data.frame(res_sc.18C)
ress_sc.18C$KEGG_ko <- row.names(ress_sc.18C)

ress_sc.18C <- merge(ress_sc.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.18C <- ress_sc.18C[order(ress_sc.18C$padj),]
ress_sc.18C$DEO <- ifelse(abs(ress_sc.18C$log2FoldChange) > 1 & ress_sc.18C$padj <= 0.05, 1, 0)

# High Ammonia
ko_sc.NH4 <- ko_sc[,c("KEGG_ko", sample_sc[sample_sc$Condition %in% c("Control", "NH4"), 1])]
sample_sc.NH4 <- sample_sc[sample_sc$Condition %in% c("Control", "NH4"), ]

dds_sc.NH4 <- DESeqDataSetFromMatrix(countData = ko_sc.NH4, 
                                     colData = sample_sc.NH4, 
                                     design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc.NH4) >= 10) >= 3
dds_sc.NH4 <- dds_sc.NH4[keep,]

dds_sc.NH4 <- DESeq(dds_sc.NH4, fitType = "local")
res_sc.NH4 <- lfcShrink(dds_sc, coef = 3, type = "ashr")
summary(res_sc.NH4)

ress_sc.NH4 <- as.data.frame(res_sc.NH4)
ress_sc.NH4$KEGG_ko <- row.names(ress_sc.NH4)

ress_sc.NH4 <- merge(ress_sc.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.NH4 <- ress_sc.NH4[order(ress_sc.NH4$padj),]
ress_sc.NH4$DEO <- ifelse(abs(ress_sc.NH4$log2FoldChange) > 1 & ress_sc.NH4$padj <= 0.05, 1, 0)

# High Sulfites
ko_sc.SO2 <- ko_sc[,c("KEGG_ko", sample_sc[sample_sc$Condition %in% c("Control", "SO2"), 1])]
sample_sc.SO2 <- sample_sc[sample_sc$Condition %in% c("Control", "SO2"), ]

dds_sc.SO2 <- DESeqDataSetFromMatrix(countData = ko_sc.SO2, 
                                     colData = sample_sc.SO2, 
                                     design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc.SO2) >= 10) >= 3
dds_sc.SO2 <- dds_sc.SO2[keep,]

dds_sc.SO2 <- DESeq(dds_sc.SO2, fitType = "local", minReplicatesForReplace = Inf)
res_sc.SO2 <- lfcShrink(dds_sc, coef = 4, type = "ashr")
summary(res_sc.SO2)

ress_sc.SO2 <- as.data.frame(res_sc.SO2)
ress_sc.SO2$KEGG_ko <- row.names(ress_sc.SO2)

ress_sc.SO2 <- merge(ress_sc.SO2, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.SO2 <- ress_sc.SO2[order(ress_sc.SO2$padj),]
ress_sc.SO2$DEO <- ifelse(abs(ress_sc.SO2$log2FoldChange) > 1 & ress_sc.SO2$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_sacc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                       list(ress_sc.18C[,c(1,9)], ress_sc.NH4[,c(1,9)], ress_sc.SO2[,c(1,9)]))
colnames(venn.df_sacc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sacc[is.na(venn.df_sacc)] <- 0

venn.plot_sacc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 0)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,2] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,3] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,4] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,4] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,3] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,2] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 3)))

venn.plot_sacc <- cbind.data.frame(venn.plot_sacc, 
                                   x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                   y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sacc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = NULL)

gg.venn_sc
#
# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                   cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col.cond[-1]) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

gg.hist_sc
##### PRUEBAS

plotCounts(dds_sc, "ko:K01549", intgroup = "Condition")

## BIOLOGICAL ENRICHMENT
DEO_lt.sc <- gsub("ko:", "", subset(ress_lt.sc, DEO == 1)$KEGG_ko)
enrich_lt.sc <- enrichKEGG(DEO_lt.sc, organism = "ko", keyType = "kegg")
kk_lt.sc <- data.frame(enrich_lt.sc, Genus = "Lachancea")


DEO_hs.sc <- gsub("ko:", "", subset(ress_hs.sc, DEO == 1)$KEGG_ko)
enrich_hs.sc <- enrichKEGG(DEO_hs.sc, organism = "ko", keyType = "kegg")
kk_hs.sc <- data.frame(enrich_hs.sc, Genus = "Hanseniaspora")

kk_gen <- merge(kk_lt.sc[,c(2,6)], kk_hs.sc[,c(2,6)], by = "Description", all = TRUE)
kk_gen[is.na(kk_gen)] <- 1
colnames(kk_gen) <- c("Description", "Lachancea", "Hanseniaspora")
kk_gen <- melt(kk_gen)

kk_gen$variable <- factor(kk_gen$variable, levels = c("Lachancea", "Hanseniaspora"))
kk_gen <- kk_gen[order(kk_gen$value),]

to.rm <- c("Huntington disease", "Oocyte meiosis", "Amyotrophic lateral sclerosis",
           "Prion disease", "Progesterone-mediated oocyte maturation", "Fanconi anemia pathway",
           "Carbon fixation in photosynthetic organisms", "Coronavirus disease - COVID-19",
           "Longevity regulating pathway - multiple species", 
           "Human T-cell leukemia virus 1 infection", 
           "Tropane, piperidine and pyridine alkaloid biosynthesis", "Terpenoid backbone biosynthesis",
           "Pathogenic Escherichia coli infection", "Synaptic vesicle cycle", "Parkinson disease")

kk_gen <- kk_gen[!kk_gen$Description %in% to.rm,]
kk_gen 

ggplot(kk_gen, aes(x = reorder(Description, log(value)), y = -log(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = "Genus") + xlab("KEGG pathways") + ylab("-log adjusted p value") +
  coord_flip()


#


#
#### DIFFERENTIAL ANALYSIS - LACHANCEA ####
ko_lt <- ko_df[,c("KEGG_ko", sample_df[sample_df$Genus == "Lachancea", 1])]
sample_lt <- sample_df[sample_df$Genus == "Lachancea",]

dds_lt <- DESeqDataSetFromMatrix(countData = ko_lt, 
                                 colData = sample_lt, 
                                 design = ~ Condition + Origin + Condition:Origin, tidy = TRUE)

keep <- rowSums(counts(dds_lt) >= 10) >= 9
dds_lt <- dds_lt[keep,]

dds_lt <- DESeq(dds_lt, fitType = "local")
resultsNames(dds_lt)

## PCA
vst_lt <- vst(dds_lt, blind = TRUE)

pcaData_lt <- plotPCA(vst_lt, intgroup = "Condition", returnData = TRUE)
percentVar_lt <- round(100 * attr(pcaData_lt, "percentVar"), 2)
pcaData_lt <- merge(pcaData_lt, sample_lt[,-4], by.x = "name", by.y = "Sample_ID")

ggplot(pcaData_lt, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_lt[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_lt[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")


## DIFFERENTIAL EXPRESSION
# Low Temperature
res_lt.18C <- lfcShrink(dds_lt, coef = 2, type = "ashr")
summary(res_lt.18C)

ress_lt.18C <- as.data.frame(res_lt.18C)
ress_lt.18C$KEGG_ko <- row.names(ress_lt.18C)

ress_lt.18C <- merge(ress_lt.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.18C <- ress_lt.18C[order(ress_lt.18C$padj),]
ress_lt.18C$DEO <- ifelse(abs(ress_lt.18C$log2FoldChange) > 1 & ress_lt.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_lt.NH4 <- lfcShrink(dds_lt, coef = 3, type = "ashr")
summary(res_lt.NH4)

ress_lt.NH4 <- as.data.frame(res_lt.NH4)
ress_lt.NH4$KEGG_ko <- row.names(ress_lt.NH4)

ress_lt.NH4 <- merge(ress_lt.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.NH4 <- ress_lt.NH4[order(ress_lt.NH4$padj),]
ress_lt.NH4$DEO <- ifelse(abs(ress_lt.NH4$log2FoldChange) > 1 & ress_lt.NH4$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_lt <- merge(ress_lt.18C[,c(1, 9)], ress_lt.NH4[,c(1, 9)], by = "KEGG_ko")
colnames(venn.df_lt) <- c("KEGG_ko", "18C", "NH4")

venn.plot_lt <- rbind.data.frame(cbind(`18C` = 0, NH4 = 0, 
                                       Counts = sum(rowSums(venn.df_lt[,-1]) == 0)),
                                 cbind(`18C` = 0, NH4 = 1, 
                                       Counts = sum(venn.df_lt[,2] == 0 & venn.df_lt[,3] == 1)),
                                 cbind(`18C` = 1, NH4 = 0, 
                                       Counts = sum(venn.df_lt[,2] == 1 & venn.df_lt[,3] == 0)),
                                 cbind(`18C` = 1, NH4 = 1, 
                                       Counts = sum(rowSums(venn.df_lt[,-1]) == 2)))

venn.plot_lt <- cbind.data.frame(venn.plot_lt, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("18C", "NH4"))
venn.out$labels <- factor(venn.out$labels, levels = c("18C", "NH4"))

ggplot() +
  geom_circle(data = venn.out, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_lt, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = NULL)

# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_lt.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_lt.NH4, DEO == 1), Condition = "NH4"))

ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col.cond[-1]) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

ress_lt.NH4
d <- plotCounts(dds_lt, gene = "ko:K17987,ko:K21997", intgroup = "Condition", 
                returnData = TRUE)

ggplot(d, aes(x = Condition, y = count)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0))




#
#### DIFFERENTIAL ANALYSIS - HANSENIASPORA ####
ko_hs <- ko_df[,c("KEGG_ko", sample_df[sample_df$Genus == "Hanseniaspora", 1])]
sample_hs <- sample_df[sample_df$Genus == "Hanseniaspora",]

dds_hs <- DESeqDataSetFromMatrix(countData = ko_hs, 
                                 colData = sample_hs, 
                                 design = ~ Condition + Origin + Condition:Origin, tidy = TRUE)

dds_hs <- DESeq(dds_hs)
resultsNames(dds_hs)

## PCA
vst_hs <- vst(dds_hs, blind = TRUE)

pcaData_hs <- plotPCA(vst_hs, intgroup = "Condition", returnData = TRUE)
percentVar_hs <- round(100 * attr(pcaData_hs, "percentVar"), 2)
pcaData_hs <- merge(pcaData_hs, sample_hs[,-4], by.x = "name", by.y = "Sample_ID")

ggplot(pcaData_hs, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_hs[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_hs[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")


## DIFFERENTIAL EXPRESSION
# Low Temperature
res_hs.18C <- lfcShrink(dds_hs, coef = 2, type = "ashr")
summary(res_hs.18C)

ress_hs.18C <- as.data.frame(res_hs.18C)
ress_hs.18C$KEGG_ko <- row.names(ress_hs.18C)

ress_hs.18C <- merge(ress_hs.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.18C <- ress_hs.18C[order(ress_hs.18C$padj),]
ress_hs.18C$DEO <- ifelse(abs(ress_hs.18C$log2FoldChange) > 1 & ress_hs.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_hs.NH4 <- lfcShrink(dds_hs, coef = 3, type = "ashr")
summary(res_hs.NH4)

ress_hs.NH4 <- as.data.frame(res_hs.NH4)
ress_hs.NH4$KEGG_ko <- row.names(ress_hs.NH4)

ress_hs.NH4 <- merge(ress_hs.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.NH4 <- ress_hs.NH4[order(ress_hs.NH4$padj),]
ress_hs.NH4$DEO <- ifelse(abs(ress_hs.NH4$log2FoldChange) > 1 & ress_hs.NH4$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_hs <- merge(ress_hs.18C[,c(1, 9)], ress_hs.NH4[,c(1, 9)], by = "KEGG_ko")
colnames(venn.df_hs) <- c("KEGG_ko", "18C", "NH4")

venn.plot_hs <- rbind.data.frame(cbind(`18C` = 0, NH4 = 0, 
                                       Counts = sum(rowSums(venn.df_hs[,-1]) == 0)),
                                 cbind(`18C` = 0, NH4 = 1, 
                                       Counts = sum(venn.df_hs[,2] == 0 & venn.df_hs[,3] == 1)),
                                 cbind(`18C` = 1, NH4 = 0, 
                                       Counts = sum(venn.df_hs[,2] == 1 & venn.df_hs[,3] == 0)),
                                 cbind(`18C` = 1, NH4 = 1, 
                                       Counts = sum(rowSums(venn.df_hs[,-1]) == 2)))

venn.plot_hs <- cbind.data.frame(venn.plot_hs, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("18C", "NH4"))
venn.out$labels <- factor(venn.out$labels, levels = c("18C", "NH4"))

ggplot() +
  geom_circle(data = venn.out, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_hs, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = NULL)

# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_hs.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_hs.NH4, DEO == 1), Condition = "NH4"))

ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col.cond[-1]) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

ress_hs.NH4
d <- plotCounts(dds_hs, gene = "ko:K17987,ko:K21997", intgroup = "Condition", 
                returnData = TRUE)

ggplot(d, aes(x = Condition, y = count)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0))




#






################################

#### DIFFERENTIAL ANALYSIS - SACCHAROMYCES (ONLY) ####
ko_sc <- ko_df[,c("KEGG_ko", sample_df[sample_df$Genus == "Saccharomyces", 1])]
ko_sc <- ko_sc[,-c(6:8)]
sample_sc <- sample_df[sample_df$Sample_ID %in% colnames(ko_sc),]

dds_sc <- DESeqDataSetFromMatrix(countData = ko_sc, 
                                 colData = sample_sc, 
                                 design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc) >= 10) >= 5
dds_sc <- dds_sc[keep,]

dds_sc <- DESeq(dds_sc, fitType = "local")
resultsNames(dds_sc)

## PCA
vst_sc <- vst(dds_sc, blind = TRUE)

pcaData_sc <- plotPCA(vst_sc, intgroup = "Condition", returnData = TRUE, ntop = 500)
percentVar_sc <- round(100 * attr(pcaData_sc, "percentVar"), 2)
pcaData_sc <- merge(pcaData_sc, sample_sc[,-4], by.x = "name", by.y = "Sample_ID")

gg.pca_sacc <- ggplot(pcaData_sc, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_sc[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sc[2],"% variance")) + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")

gg.pca_sacc
ggsave("Workflow/3_RNAseq-meta/Outputs/Figures/sacc_PCA.png", 
       gg.pca_sacc, width = 7.5, height = 6, dpi = 300)

## DIFFERENTIAL EXPRESSION
# Low Temperature
res_sc.18C <- lfcShrink(dds_sc, coef = 2, type = "ashr")
summary(res_sc.18C)

ress_sc.18C <- as.data.frame(res_sc.18C)
ress_sc.18C$KEGG_ko <- row.names(ress_sc.18C)

ress_sc.18C <- merge(ress_sc.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.18C <- ress_sc.18C[order(ress_sc.18C$padj),]
ress_sc.18C$DEO <- ifelse(abs(ress_sc.18C$log2FoldChange) > 1 & ress_sc.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_sc.NH4 <- lfcShrink(dds_sc, coef = 3, type = "ashr")
summary(res_sc.NH4)

ress_sc.NH4 <- as.data.frame(res_sc.NH4)
ress_sc.NH4$KEGG_ko <- row.names(ress_sc.NH4)

ress_sc.NH4 <- merge(ress_sc.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.NH4 <- ress_sc.NH4[order(ress_sc.NH4$padj),]
ress_sc.NH4$DEO <- ifelse(abs(ress_sc.NH4$log2FoldChange) > 1 & ress_sc.NH4$padj <= 0.05, 1, 0)

# High Sulfites
res_sc.SO2 <- lfcShrink(dds_sc, coef = 4, type = "ashr")
summary(res_sc.SO2)

ress_sc.SO2 <- as.data.frame(res_sc.SO2)
ress_sc.SO2$KEGG_ko <- row.names(ress_sc.SO2)

ress_sc.SO2 <- merge(ress_sc.SO2, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.SO2 <- ress_sc.SO2[order(ress_sc.SO2$padj),]
ress_sc.SO2$DEO <- ifelse(abs(ress_sc.SO2$log2FoldChange) > 1 & ress_sc.SO2$padj <= 0.05, 1, 0)

## SUMMARY
# DE Orthologs - Venn diagram
venn.df_sacc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                       list(ress_sc.18C[,c(1,9)], ress_sc.NH4[,c(1,9)], ress_sc.SO2[,c(1,9)]))
colnames(venn.df_sacc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sacc[is.na(venn.df_sacc)] <- 0

venn.plot_sacc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 0)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,2] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,3] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,4] == 1 & rowSums(venn.df_sacc[,-1]) == 1)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                         Counts = sum(venn.df_sacc[,4] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,3] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(venn.df_sacc[,2] == 0 & rowSums(venn.df_sacc[,-1]) == 2)),
                                   cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                         Counts = sum(rowSums(venn.df_sacc[,-1]) == 3)))

venn.plot_sacc <- cbind.data.frame(venn.plot_sacc, 
                                   x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                   y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sacc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = NULL)

gg.venn_sc
#
# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                   cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = col.cond[-1]) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

gg.hist_sc
##### PRUEBAS

plotCounts(dds_sc, "ko:K01487", intgroup = "Condition")

## BIOLOGICAL ENRICHMENT
DEO_lt.sc <- gsub("ko:", "", subset(ress_lt.sc, DEO == 1)$KEGG_ko)
enrich_lt.sc <- enrichKEGG(DEO_lt.sc, organism = "ko", keyType = "kegg")
kk_lt.sc <- data.frame(enrich_lt.sc, Genus = "Lachancea")


DEO_hs.sc <- gsub("ko:", "", subset(ress_hs.sc, DEO == 1)$KEGG_ko)
enrich_hs.sc <- enrichKEGG(DEO_hs.sc, organism = "ko", keyType = "kegg")
kk_hs.sc <- data.frame(enrich_hs.sc, Genus = "Hanseniaspora")

kk_gen <- merge(kk_lt.sc[,c(2,6)], kk_hs.sc[,c(2,6)], by = "Description", all = TRUE)
kk_gen[is.na(kk_gen)] <- 1
colnames(kk_gen) <- c("Description", "Lachancea", "Hanseniaspora")
kk_gen <- melt(kk_gen)

kk_gen$variable <- factor(kk_gen$variable, levels = c("Lachancea", "Hanseniaspora"))
kk_gen <- kk_gen[order(kk_gen$value),]

to.rm <- c("Huntington disease", "Oocyte meiosis", "Amyotrophic lateral sclerosis",
           "Prion disease", "Progesterone-mediated oocyte maturation", "Fanconi anemia pathway",
           "Carbon fixation in photosynthetic organisms", "Coronavirus disease - COVID-19",
           "Longevity regulating pathway - multiple species", 
           "Human T-cell leukemia virus 1 infection", 
           "Tropane, piperidine and pyridine alkaloid biosynthesis", "Terpenoid backbone biosynthesis",
           "Pathogenic Escherichia coli infection", "Synaptic vesicle cycle", "Parkinson disease")

kk_gen <- kk_gen[!kk_gen$Description %in% to.rm,]
kk_gen 

ggplot(kk_gen, aes(x = reorder(Description, log(value)), y = -log(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = "Genus") + xlab("KEGG pathways") + ylab("-log adjusted p value") +
  coord_flip()


#


#
#### TRIES ####


plot(res_sc.18C$baseMean+1, -log10(res_sc.18C$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))


use <- res_sc.18C$baseMean > metadata(res_sc.18C)$filterThreshold
h1 <- hist(res_sc.18C$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_sc.18C$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


resNorm <- lfcShrink(dds_sc, name = "Condition_18C_vs_Control", coef = 2, type = "normal")
resAsh <- lfcShrink(dds_sc, contrast = c("Condition", "18C", "Control"), coef = 2, type = "ashr")


DESeq2::plotMA(res_sc.18C)
DESeq2::plotMA(resIHW)
DESeq2::plotMA(resNorm)
DESeq2::plotMA(resAsh)


ntd <- normTransform(dds_sc)
library("vsn")
meanSdPlot(assay(ntd))

vst_sc <- vst(dds_sc, blind = TRUE)
meanSdPlot(assay(vst_sc))

rld_sc <- rlog(dds_sc, blind = TRUE)
meanSdPlot(assay(rld_sc))

plotPCA(vst_sc, intgroup = "Condition")
plotPCA(rld_sc, intgroup = "Condition")


plotDispEsts(dds_sc)


plot(metadata(res_sc.18C)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_sc.18C)$lo.fit, col="red")
abline(v=metadata(res_sc.18C)$filterTheta)
