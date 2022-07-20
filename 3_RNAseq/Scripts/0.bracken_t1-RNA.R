# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T1-GM samples

library(reshape2)

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### MAKE ASV AND TAX TABLES ####
brack_df <- t(read.table("Inputs/bracken.report.txt"))
row.names(brack_df) <- gsub("\\.", "-", row.names(brack_df))

tax_df <- cbind.data.frame(Genus = colsplit(colnames(brack_df), " ", c("Genus", "species"))[,1],
                           Species = colnames(brack_df))

row.names(tax_df) <- paste("Sp", 1:97, sep = "-")
colnames(brack_df) <- paste("Sp", 1:97, sep = "-")

#
#### SAVE DATA ####
save.image("dada2_t1-RNA.RData")

saveRDS(tax_df, "Outputs/tax_t1-RNA.rds")
saveRDS(brack_df, "Outputs/ASV_t1-RNA.rds")

#