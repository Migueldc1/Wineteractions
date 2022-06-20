# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T1-SGM samples

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#### LIBRARIES ####
library(dada2)
library(ShortRead)
library(Biostrings)

#
#### CUSTOM FUNCTIONS ####
getN <- function(x) sum(getUniques(x))

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#
#### GETTING READY ####
path <- "Inputs/ITS_reads/t1-SGM/"

fnFs <- sort(list.files(file.path(path, "R1"), pattern = ".fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(file.path(path, "R2"), pattern = ".fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- gsub("NGS025-21-RUN-2-", "", sample.names)

#Identify Primers
FWD <- "TCCTCCGCTTATTGATATGC"
REV <- "GTGARTCATCGAATCTTTG"

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#
#### REMOVE Ns ####
fnFs.filtN <- file.path(path, "Remove.primers/R1/filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "Remove.primers/R2/filtN", basename(fnRs))

filtN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, verbose = TRUE)

#
#### REMOVE PRIMERS ####
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

#Set cutadapt path
cutadapt <- "C:/Users/Laboratorio14/AppData/Local/Programs/Python/Python310/Scripts/cutadapt.exe"
system2(cutadapt, args = "--version")

fnFs.cut <- file.path(path, "Remove.primers/R1/cutadapt", basename(fnFs))
fnRs.cut <- file.path(path, "Remove.primers/R2/cutadapt", basename(fnRs))

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-g", REV, "-G", FWD,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i]))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]]))

#
#### INSPECT QUALITY PROFILES ####
plotQualityProfile(fnFs.cut[sample(1:length(sample.names), 3)])
plotQualityProfile(fnRs.cut[sample(1:length(sample.names), 3)])


#
#### FILTER AND TRIM ####
filtFs <- file.path(path, "Remove.primers/R1/filtered", basename(fnFs))
filtRs <- file.path(path, "Remove.primers/R2/filtered", basename(fnRs))

out.cut <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, truncLen = c(240,200),
                         maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                         compress = TRUE, verbose = TRUE) 
head(out.cut)

#
#### LEARN ERROR RATES ####
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ = TRUE)

#
#### SAMPLE INFERENCE ####
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

#
#### MERGE PAIRED READS ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

#
#### CONSTRUCT SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
hist(nchar(getSequences(seqtab)))

row.names(seqtab) <- sapply(strsplit(row.names(seqtab), "_"), `[`, 1)
row.names(seqtab) <- gsub("NGS025-21-RUN-2-", "", row.names(seqtab))

#
#### REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#
#### TRACK READS ####
track.cut <- cbind.data.frame(filtN[,1], out.cut[,2], sapply(dadaFs, getN), 
                              sapply(dadaRs, getN), sapply(mergers, getN), 
                              rowSums(seqtab.nochim))

colnames(track.cut) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.cut) <- sample.names

track.cut <- cbind.data.frame(track.cut, perc = track.cut[,6]*100/track.cut[,1])

#
#### ASSIGN TAXONOMY ####
taxa.cut <- assignTaxonomy(seqtab.nochim, 
                           "Inputs/Databases/sh_general_release_dynamic_10.05.2021.fasta")

#
#### SAVE DATA ####
save.image("dada2_t1-SGM.RData")

saveRDS(taxa.cut, "Outputs/tax_t1-SGM.rds")
saveRDS(seqtab.nochim, "Outputs/ASV_t1-SGM.rds")
write.table(track.cut, "Outputs/track_t1-SGM.txt", sep = "\t", dec = ",")

