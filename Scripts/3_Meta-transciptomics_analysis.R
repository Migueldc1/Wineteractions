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

## RNA RELATED DATA

#KO Table
ko_df <- read.emapper(path.count = "Data/Meta-transcriptomics/Count/eggnog", 
                      path.orth = "Data/Meta-transcriptomics/Annotation/eggnog")
ko_df.f <- ko_df[!grepl(",", ko_df$KEGG_ko),]

#Annotation Tables
kegg_df <- read.table("Data/Meta-transcriptomics/Annotation/KEGG_names.txt", sep = "\t", header = TRUE, quote = "")
go_df <- read.GOdata(".emapper.annotations", "Data/Meta-transcriptomics/Annotation/eggnog", c("KEGG_ko", "GOs"))


## TAXONOMY DATA
tax_rna <- read.table("Data/Meta-transcriptomics/bracken.report.txt", check.names = FALSE)

tax_KO <- read.table("Data/Meta-transcriptomics/contig_tax.txt", check.names = FALSE)

tax_ITS <- as.data.frame(readRDS("Data/Sequencing/Outputs/tax_SGM.rds"))
tax_ITS[is.na(tax_ITS)] <- "Unidentified"
asv_ITS <- readRDS("Data/Sequencing/Outputs/ASV_SGM.rds")
asv_ITS <- apply(asv_ITS, 1, function(x) x/sum(x))


## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))


## COLORS
col_cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col_genus <- c("#caf55d", "#5df5cc", "#cc3939", "#cab2d6", "#8da0cb", "#d5eb26", "#1b9e77",
                        "#f78e4d", "#e6d8bd", "#6a3d9a", "#c05b17")

#
#### TAXONOMIC COMPOSITION ####

tax_SGM <- aggregate(asv_ITS, list(tax_ITS$Genus), sum)
colnames(tax_SGM)[1] <- "Genus"
tax_SGM$Genus <- gsub("g__", "", tax_SGM$Genus)

tax_rna <- cbind.data.frame(Genus = colsplit(row.names(tax_rna), pattern = " ", names = c("Genus", "to.rm"))[1], tax_rna)
tax_rna <- aggregate(. ~ Genus, tax_rna, sum)

tax_KO <- cbind.data.frame(Genus = row.names(tax_KO), tax_KO)


tax_plot <- rbind(cbind.data.frame(melt(tax_SGM), Assay = "ITS"),
                  cbind.data.frame(melt(tax_rna), Assay = "RNA reads"),
                  cbind.data.frame(melt(tax_KO), Assay = "Contigs"))

tax_plot[tax_plot$value < 0.025, "Genus"] <- "Other"
tax_plot <- aggregate(value ~ Genus + variable + Assay, tax_plot, sum)

orderG <- levels(factor(tax_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

tax_plot <- cbind.data.frame(tax_plot, colsplit(tax_plot$variable, "-", c("Origin", "Farming", "Condition")))
tax_plot$Sample_ID <- paste(tax_plot$Farming, tax_plot$Condition, sep = " ")

tax_plot$Origin <- factor(tax_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
tax_plot$Assay <- factor(tax_plot$Assay, levels = c("ITS", "RNA reads", "Contigs"))
tax_plot$Sample_ID <- factor(tax_plot$Sample_ID, levels = c("CONV Control", "CONV 18C", "CONV NH4", "CONV SO2", 
                                                            "ECO Control", "ECO 18C", "ECO NH4", "ECO SO2"))

gg.tax <- ggplot(tax_plot, 
       aes(x = Sample_ID, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = col_genus) +
  facet_grid(Assay ~ Origin) +
  guides(fill=guide_legend(nrow = 2)) +
  xlab("") + ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.tax

#
#### PCA - GLOBAL ####

tax_max <- cbind.data.frame(Sample_ID = colnames(tax_rna[,-1]), Genus = tax_rna$Genus[apply(tax_rna[,-1], 2, which.max)],
                            value = apply(tax_rna[,-1], 2, max))
tax_max$Genus <- ifelse(tax_max$value >= 0.9, tax_max$Genus, "Other")
tax_max$Genus <- ifelse(tax_max$Genus %in% c("Saccharomyces", "Hanseniaspora", "Lachancea"), tax_max$Genus, "Other")

sample_df <- merge(sample_df, tax_max[,-3], by = "Sample_ID")

dds_tot <- DESeqDataSetFromMatrix(countData = ko_df.f, 
                                  colData = sample_df, 
                                  design = ~ Condition + Genus + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_tot)) >= 10
dds_tot <- dds_tot[keep,]

dds_tot <- DESeq(dds_tot, fitType = "local")
resultsNames(dds_tot)

## PCA
vst_tot <- vst(dds_tot, blind = TRUE)

pcaData_tot <- plotPCA(vst_tot, intgroup = "Genus", returnData = TRUE)
pcaData_tot$Genus <- factor(pcaData_tot$Genus, levels = c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"))
percentVar_tot <- round(100 * attr(pcaData_tot, "percentVar"), 2)

# Samples are separated by dominant species
gg.pca_gen <- ggplot(pcaData_tot, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_tot[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_tot[2],"% variance")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#cc3939", "#8da0cb", "#1b9e77", "gray70"))

gg.pca_gen

#
#### DIFFERENTIAL ANALYSIS - GLOBAL ####
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
  theme(legend.position = "none", 
        legend.text = element_text(size = 13, color  = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = NULL)

gg.venn_gen

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
  theme(legend.position = "none",
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_gen

#
#### EXPORT FIGURE 3 ####

gg.figure3 <- plot_grid(gg.tax, plot_grid(gg.pca_gen, plot_grid(gg.venn_gen, gg.go_gen, ncol = 1, rel_heights = c(0.4, 1),
                                                                labels = c("C", "D"), label_size = 18), labels = "B", label_size = 18), 
                        labels = "A", label_size = 18, ncol = 1)
gg.figure3

ggsave("Figures/Figure_3.png", gg.figure3, bg = "white", width = 16, height = 16)

#