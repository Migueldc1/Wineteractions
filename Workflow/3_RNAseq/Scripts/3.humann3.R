# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - RNAseq analysis

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(reshape2)
library(Maaslin2)

#load("huamnn.Rdata")

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df <- cbind.data.frame(Sample_ID = sample_df[,1],
                              colsplit(sample_df[,1], "-",
                                       names = c("Origin", "Farming", "Condition")),
                              sample_df[,-1])
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

## Bracken Taxonomy
asv_RNA <- readRDS("Outputs/ASV_t1-RNA.rds")
asv_RNA <- apply(asv_RNA, 1, function(x) x/sum(x))

tax_RNA <- readRDS("Outputs/tax_t1-RNA.rds")
tax_RNA <- cbind.data.frame(tax_RNA, Id = row.names(tax_RNA))

asv_RNA <- melt(asv_RNA)
colnames(asv_RNA) <- c("Id", "Sample_ID", "value")
asv_RNA <- merge(asv_RNA, tax_RNA[,c(1,3)])

## Dominant Genus
dom_sp <- aggregate(asv_RNA$value, list(asv_RNA$Sample_ID), max)
colnames(dom_sp) <- c("Sample_ID", "value")
dom_sp <- merge(dom_sp, asv_RNA[,2:4], by = c("Sample_ID", "value"))

sample_df <- merge(sample_df, dom_sp[,c(1,3)], by = "Sample_ID")
sample_df$Genus <- relevel(as.factor(sample_df$Genus), ref = "Saccharomyces")

## HUMANN3 PAHTABUNDANCES
path_hu <- read.delim("Inputs/RNA_data/wineteractions_pathabundance_relab.tsv", header = TRUE, 
                      sep = "\t", row.names = 1)
colnames(path_hu) <- gsub("\\.", "-", colnames(path_hu))
colnames(path_hu) <- gsub("_Abundance", "", colnames(path_hu))

path_hu.tot <- path_hu[!grepl("\\|", row.names(path_hu)),]
path_hu.tot <- apply(path_hu.tot, 2, function(x) x/sum(x))

#
#### ANOVAS ####


datos <- t(as.matrix(path_hu.tot))
datos <- cbind.data.frame(Sample_ID = row.names(datos), datos)
datos <- merge(sample_df[,c(1,17)], datos, by = "Sample_ID")


result <- data.frame()

for (n.path in 3:ncol(datos)) {
  
  path <- colnames(datos)[n.path]
  dato_fresco <- cbind.data.frame(datos[,1:2],
                                  value = datos[,n.path])
  
  fitmrna <- lm(value ~ Genus, data = dato_fresco)
  stat <- anova(fitmrna)
  
  result <- rbind.data.frame(result,
                             cbind(Path = path,
                                   Saccharomyces_AvgExp = mean(datos[datos$Genus == "Saccharomyces", path]),
                                   Lachancea_AvgExp = mean(datos[datos$Genus == "Lachancea", path]),
                                   Hanseniaspora_AvgExp = mean(datos[datos$Genus == "Hanseniaspora", path]),
                                   p.value = stat$`Pr(>F)`[1]))
  
  
}







path
fitmrna <- lm("1CMET2-PWY: folate transformations III (E. coli)" ~ Genus, data = datos)
stat <- anova(fitmrna)

result <- cbind.data.frame(Path = "VALSYN-PWY: L-valine biosynthesis",
                           Saccharomyces_AvgExp = mean(datos[datos$Genus == "Saccharomyces", "VALSYN-PWY: L-valine biosynthesis"]),
                           Lachancea_AvgExp = mean(datos[datos$Genus == "Lachancea", "VALSYN-PWY: L-valine biosynthesis"]),
                           Hanseniaspora_AvgExp = mean(datos[datos$Genus == "Hanseniaspora", "VALSYN-PWY: L-valine biosynthesis"]),
                           p.value = stat$`Pr(>F)`[1])






















result[linea, "Condition_1_AvgExp"] <- mean(as.numeric(data[linea, Condition_1]))
result[linea, "Condition_2_AvgExp"] <- mean(as.numeric(data[linea, Condition_2]))
result[linea, "p.value"] <- stat$`Pr(>F)`[3]

rownames(result)[linea] <- rownames(data)[linea]
boxplot(as.numeric(data[linea,Condition_1]), as.numeric(data[linea,Condition_2]),
        col = c("salmon", "lightblue"), main = paste(strtrim(rownames(data)[linea], 25),"...", sep= ""),
        ylab = "Relative abundance", ylim = c(-0.05, 0.1))

legend("topright", paste("Pvalue: ", round(stat$`Pr(>F)`[3], digits = 3)), bty = "n")















length(unique(datos$path))







for(linea in 1:nrow(path_hu.tot)){
  datos <- NA
  datos <- data.frame(exp = as.numeric(data[linea,]),
                      Farming = NA,
                      Condition = NA)
  
  datos[Condition_1, 2] <- "Condition_1"
  datos[Condition_2, 2] <- "Condition_2"
  
  datos[Condition_A, 3] <- "Condition_A"
  datos[Condition_B, 3] <- "Condition_B"
  datos[Condition_C, 3] <- "Condition_C"
  datos[Condition_D, 3] <- "Condition_D"
  
  datos <- na.omit(datos)
  fitmrna <- lm(exp~Farming*Condition, data = datos)
  stat <- anova(fitmrna)
  
  result[linea, "Condition_1_AvgExp"] <- mean(as.numeric(data[linea, Condition_1]))
  result[linea, "Condition_2_AvgExp"] <- mean(as.numeric(data[linea, Condition_2]))
  result[linea, "p.value"] <- stat$`Pr(>F)`[3]
  
  rownames(result)[linea] <- rownames(data)[linea]
  boxplot(as.numeric(data[linea,Condition_1]), as.numeric(data[linea,Condition_2]),
          col = c("salmon", "lightblue"), main = paste(strtrim(rownames(data)[linea], 25),"...", sep= ""),
          ylab = "Relative abundance", ylim = c(-0.05, 0.1))
  
  legend("topright", paste("Pvalue: ", round(stat$`Pr(>F)`[3], digits = 3)), bty = "n")
  
}


#### MaAsLin2 ####

row.names(sample_df) <- sample_df$Sample_ID
path_hu.tot <- t(path_hu.tot)

sample_df.m2 <- subset(sample_df, Genus == "Saccharomyces" | Genus == "Lachancea" | 
                         Genus == "Hanseniaspora")

sample_df.m2$Genus <- relevel(as.factor(as.character(sample_df.m2$Genus)), ref = "Saccharomyces")

path_hu.m2 <- path_hu.tot[row.names(path_hu.tot) %in% row.names(sample_df.m2),]
path_hu.m2 <- path_hu.m2[,colMeans(path_hu.m2) >= 0.0001]

fit_m2 <- Maaslin2(input_data = path_hu.m2,
                   input_metadata = sample_df.m2,
                   output = "MaAsLin2",
                   min_abundance = 0.0001,
                   min_prevalence = 0.1,
                   normalization = "NONE",
                   transform = "LOG",
                   analysis_method = "LM",
                   max_significance = 0.05,
                   fixed_effects = c("Genus"),
                   reference = c("Genus,Saccharomyces"))

#
### Saccharomyces vs. Lachancea
sample_df.sl <- subset(sample_df.m2, Genus == "Saccharomyces" | Genus == "Lachancea")

path_hu.sl <- path_hu.tot[row.names(path_hu.tot) %in% row.names(sample_df.sl),]
path_hu.sl <- path_hu.sl[,colMeans(path_hu.sl) >= 0.0001]

fit_sl <- Maaslin2(input_data = path_hu.sl,
                   input_metadata = sample_df.sl,
                   output = "Sacc_vs_Lach",
                   min_abundance = 0.0001,
                   min_prevalence = 0.1,
                   normalization = "NONE",
                   transform = "NONE",
                   analysis_method = "LM",
                   max_significance = 0.05,
                   fixed_effects = c("Genus"),
                   reference = c("Genus,Saccharomyces"))









