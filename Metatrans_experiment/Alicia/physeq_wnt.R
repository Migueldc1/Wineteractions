#### LIBRARIES ####
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(vegan)
library(hillR)
library(igraph)
#library(dada2) #Get nucleotidic sequences
#library(phangorn)
#library(msa) #Perfom multiple alignment

setwd("C:/Users/Laboratorio14/Desktop/Miguel/Wineteractions/Fungi")
load("physeq_fun.RData")

#
#### FUNCTIONS ####
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#
#### DATA PREPROCESSING ####
asv_fun <- readRDS("ASV_cut.rds")
asv_fun <- asv_fun[row.names(asv_fun) %in% meta.ns$Sample_ID, ]
tasv_fun <- asv_fun/rowSums(asv_fun)

tax_fun <- readRDS("tax_cut.rds")
tax_fun <- cbind.data.frame(tax_fun, Id = row.names(tax_fun))
tax_fun$ASV <- paste("ASV", 1:nrow(tax_fun), sep = "-")

meta <- read.table("meta.txt", sep = "\t", header = TRUE)
row.names(meta) <- meta$Sample_ID
meta$name <- paste(meta$Origin, meta$Sample_name)

meta$Time <- factor(meta$Time, levels = c("Initial", "Middle", "Final"))
meta$Origin <- factor(meta$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B",
                                              "R3C"))
meta$Region <- factor(meta$Region, levels = c("RiberaGuadiana", "Madrid", "LaMancha", 
                                              "Valdepe?as", "LaRioja"))

meta.ns <- subset(meta, Condition != "sintetico")

phy_fun <- phyloseq(otu_table(asv_fun, taxa_are_rows = FALSE), tax_table(tax_fun), 
                    sample_data(meta))

phy_fun@tax_table[is.na(phy_fun@tax_table)] <- "Unidentified"

phy_fun <- subset_samples(phy_fun, order > 54)
phyt_fun <- transform_sample_counts(phy_fun, function(x) x / sum(x))

#
#### Phylogenetic tree ####
seqs_fun <- getSequences(colnames(asv_fun))
names(seqs_fun) <- seqs_fun

mult_fun <- msa(seqs_fun, method = "ClustalW", type = "dna", order = "input")

phang.align_fun <- as.phyDat(mult_fun, type = "DNA", names = seqs_fun)
dm_fun <- dist.ml(phang.align_fun)
treeNJ_fun <- NJ(dm_fun) # Note, tip order != sequence order
fit_fun <- pml(treeNJ_fun, data = phang.align_fun)

fitGTR_fun <- update(fit_fun, k = 4, inv = 0.2)
fitGTR_fun <- optim.pml(fitGTR_fun, model = "GTR", optInv = TRUE, optGamma = TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

#
#### Get most abundant genus at final stage ####
max <- apply(tasv_fun, 1, function(x) which.max(x))
max <- data.frame(max)
max$Sample_ID <- row.names(max)
max$Id <- row.names(xx)[max$max]

max <- merge(max, tax_fun[,c(6,8)])
max <- merge(max, meta, by = "Sample_ID")
max$name <- paste(max$Origin, max$Sample_name)

max <- subset(max, Time == "Final")[,c(14, 4)]

max$most.ab <- ifelse(max$Genus == "g__Lachancea" | max$Genus == "g__Saccharomyces", 
                     as.character(max$Genus), 
                     "Other")

meta <- merge(meta, max[,c(1,3)], by = "name", all.x = TRUE)

#
#### TAXONOMIC ASSIGNATION ####
##Phylum
glomP_fun <- tax_glom(phyt_fun, taxrank = 'Phylum', NArm = FALSE)
dataP_fun <- psmelt(glomP_fun)
dataP_fun$Phylum <- as.character(dataP_fun$Phylum)
dataP_fun$Phylum[dataP_fun$Abundance < 0.01] <- "Other"

orderP <- levels(factor(dataP_fun$Phylum))
orderP <- orderP[! orderP %in% c("Unidentified")]
orderP <- append(orderP, c("Unidentified"))

set.seed(123)
ggplot(dataP_fun, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = orderP))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = sample(col_vector, 16)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust=1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Sample") + ylab("Abundance")

##Genus
glomG_fun <- tax_glom(phyt_fun, taxrank = 'Genus', NArm = FALSE)
dataG_fun <- psmelt(glomG_fun)
dataG_fun$Genus <- as.character(dataG_fun$Genus)
dataG_fun$Genus[dataG_fun$Abundance < 0.05] <- "Other"
dataG_fun$Genus[is.na(dataG_fun$Genus)] <- "Unidentified"

orderG <- levels(factor(dataG_fun$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

dataG_fun.plot <- aggregate(dataG_fun$Abundance, list(dataG_fun$Sample, dataG_fun$Genus), sum)
colnames(dataG_fun.plot) <- c("Sample_ID", "Genus", "Abundance")

dataG_fun.plot <- merge(dataG_fun.plot, meta, by = "Sample_ID")

set.seed(123)
ggplot(dataG_fun.plot, 
       aes(x = Sample_name, y = Abundance, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = sample(col_vector, 18)) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("Sample") + ylab("Abundance") + facet_wrap(~Time + Origin, nrow = 3) +
  guides(fill = guide_legend(nrow = 3))

#
#### ALPHA DIVERISTY: Dynamics ####
alpha_fun <- cbind.data.frame(# Hill based taxonomic alpha diversity
                              t.q0 = hill_taxa(tasv_fun, q = 0),
                              t.q1 = hill_taxa(tasv_fun, q = 1),
                              t.q2 = hill_taxa(tasv_fun, q = 2),
                              
                              # Hill based phylogenetic alpha diversity
                              p.q0 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 0),
                              p.q1 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 1),
                              p.q2 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 2))

alpha_fun$Sample_ID <- row.names(alpha_fun)

alpha_fun <- merge(alpha_fun, subset(meta, order > 54)[,-c(7,8)], by = "Sample_ID")

alpha_fun.plot <- melt(alpha_fun)
alpha_fun.plot$diversity <- ifelse(startsWith(as.character(alpha_fun.plot$variable), "t."),
                                   "Taxonomic", "Phylogenetic")
alpha_fun.plot$diversity <- factor(alpha_fun.plot$diversity, 
                                   levels = c("Taxonomic", "Phylogenetic"))

alpha_fun.plot$q <- as.factor(substr(alpha_fun.plot$variable, 
                                     nchar(as.character(alpha_fun.plot$variable)),
                                     nchar(as.character(alpha_fun.plot$variable))))

alpha_fun.plot$name <- paste(alpha_fun.plot$Origin, alpha_fun.plot$Sample_name, 
                             alpha_fun.plot$q)

xx <- subset(alpha_fun.plot, diversity == "Taxonomic")
xx <- cbind.data.frame(xx, 
                       value2 = rbind.data.frame(cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value)))
xxx <- cbind.data.frame(xx[,-14], 
                        value2 = rbind.data.frame(cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value)))


ggplot(subset(alpha_fun.plot, diversity == "Taxonomic"), 
       aes(x = Time, y = value, shape = Management, color = Condition, group = name)) + 
  geom_point(position = position_dodge(width = 0.3), size = 2) + 
  geom_line(position = position_dodge(width = 0.3), size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~q + Origin, nrow = 3)

#
#### ALPHA DIVERISTY: Normalized Plots ####
alpha_fun.norm <- alpha_fun[,-c(5,6,7)]

alpha_fun.norm$t.q0 <- ave(alpha_fun.norm$t.q0, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)
alpha_fun.norm$t.q1 <- ave(alpha_fun.norm$t.q1, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)
alpha_fun.norm$t.q2 <- ave(alpha_fun.norm$t.q2, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)

alpha_fun.norm_plot <- melt(alpha_fun.norm)

ggplot(alpha_fun.norm_plot, 
       aes(x = Time, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~variable + Condition, nrow = 3)

#
#### ALPHA DIVERSITY: Normalized dynamics plot ####
ggplot(alpha_fun.norm_plot, 
       aes(x = Time, y = value, shape = Management, color = Condition, group = Sample_name)) + 
  geom_point(position = position_dodge(width = 0.3), size = 2) + 
  geom_line(position = position_dodge(width = 0.3), size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~variable + Origin, nrow = 3)



#
#### ALPHA DIVERISTY: Diversity loss: Raw / Raw.prop ####
alpha_fun.plot <- subset(alpha_fun.plot, diversity == "Taxonomic")

alpha_fun.plot$name2 <- paste(alpha_fun.plot$Origin, alpha_fun.plot$Sample_name,
                              alpha_fun.plot$variable)

alpha_fun.dyn <- subset(alpha_fun.plot, Time == "Initial")
colnames(alpha_fun.dyn)[10] <- "Initial"

alpha_fun.dyn <- merge(alpha_fun.dyn, subset(alpha_fun.plot, Time == "Middle")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_fun.dyn)[15] <- "Middle"

alpha_fun.dyn <- merge(alpha_fun.dyn, subset(alpha_fun.plot, Time == "Final")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_fun.dyn)[16] <- "Final"

alpha_fun.dyn$`Initial-Final` <- (alpha_fun.dyn$Initial-alpha_fun.dyn$Final)/alpha_fun.dyn$Initial
alpha_fun.dyn$`Initial-Middle` <- (alpha_fun.dyn$Initial-alpha_fun.dyn$Middle)/alpha_fun.dyn$Initial

alpha_fun.dyn <- alpha_fun.dyn[,c(3,4,6,8,10,12,13,17,18)]
alpha_fun.dyn.plot <- melt(alpha_fun.dyn)
colnames(alpha_fun.dyn.plot)[5] <- "div.index"

alpha_fun.dyn.plot$variable <- factor(alpha_fun.dyn.plot$variable, 
                                      levels = c("Initial-Middle", "Initial-Final"))

ggplot(subset(xx, diversity == "Taxonomic"), 
       aes(x = Condition, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Diversity Change") + facet_wrap(~div.index + variable, nrow = 3)


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q0"),
                       cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q1"),
                       cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### ALPHA DIVERISTY: Diversity loss: Normalized ####
alpha_fun.norm_plot$name <- paste(alpha_fun.norm_plot$Origin, alpha_fun.norm_plot$Sample_name,
                                  alpha_fun.norm_plot$variable)

alpha_fun.dyn_norm <- subset(alpha_fun.norm_plot, Time == "Initial")
colnames(alpha_fun.dyn_norm)[10] <- "Initial"

alpha_fun.dyn_norm <- merge(alpha_fun.dyn_norm, 
                            subset(alpha_fun.norm_plot, Time == "Middle")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_fun.dyn_norm)[12] <- "Middle"

alpha_fun.dyn_norm <- merge(alpha_fun.dyn_norm, 
                            subset(alpha_fun.norm_plot, Time == "Final")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_fun.dyn_norm)[13] <- "Final"

alpha_fun.dyn_norm$`Initial-Final` <- (alpha_fun.dyn_norm$Initial-alpha_fun.dyn_norm$Final)
alpha_fun.dyn_norm$`Initial-Middle` <- (alpha_fun.dyn_norm$Initial-alpha_fun.dyn_norm$Middle)

alpha_fun.dyn_norm <- alpha_fun.dyn_norm[,c(1:10,14,15)]
alpha_fun.dyn_norm.plot <- melt(alpha_fun.dyn_norm)
colnames(alpha_fun.dyn_norm.plot)[10] <- "div.index"

alpha_fun.dyn_norm.plot$variable <- factor(alpha_fun.dyn_norm.plot$variable, 
                                      levels = c("Initial-Middle", "Initial-Final"))

ggplot(alpha_fun.dyn_norm.plot, 
       aes(x = Condition, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Diversity Change") + facet_wrap(~div.index + variable, nrow = 3)


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q0"),
                        cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q1"),
                        cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### BETA DIVERSITY ####
phyt_bray <- phyloseq::distance(phyt_fun, method = "bray")
nMDS_fun <- metaMDS(phyt_bray)
nMDS_fun$stress

nMDS_fun.plot <- as.data.frame(scores(nMDS_fun))
nMDS_fun.plot$Sample_ID <- rownames(nMDS_fun.plot)
nMDS_fun.plot <- merge(nMDS_fun.plot, meta, by = "Sample_ID")
nMDS_fun.plot$group <- paste(nMDS_fun.plot$Origin, nMDS_fun.plot$Sample_name)
nMDS_fun.plot$order <- as.numeric(nMDS_fun.plot$Time)
nMDS_fun.plot <- nMDS_fun.plot[order(nMDS_fun.plot$order),]

ggplot(nMDS_fun.plot) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Time, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


pcoa_fun <- wcmdscale(phyt_bray)

pcoa_fun.plot <- as.data.frame(scores(pcoa_fun))
pcoa_fun.plot$Sample_ID <- rownames(pcoa_fun.plot)
pcoa_fun.plot <- merge(pcoa_fun.plot, meta, by = "Sample_ID")
pcoa_fun.plot$group <- paste(pcoa_fun.plot$Origin, pcoa_fun.plot$Sample_name)
pcoa_fun.plot$order <- as.numeric(pcoa_fun.plot$Time)
pcoa_fun.plot <- pcoa_fun.plot[order(pcoa_fun.plot$order),]

ggplot(pcoa_fun.plot) + 
  geom_point(aes(x = Dim1, y = Dim2, color = Time, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


set.seed(1)
adonis(phyt_bray ~ Time + Region + Management:Region, data = meta.ns)

meta.ns <- meta.ns[match(sample_names(phyt_fun), meta.ns$Sample_ID),]

meta.phy <- data.frame(sample_data(phy_fun))


beta <- betadisper(phyt_bray, meta.ns$Time, type = "centroid")
permutest(beta)
plot(beta)

meta.ns$Sample_ID

#
#### BETA DIVERSITY: CONSTRAINED ####
adonis2(as.matrix(phyt_bray) ~ Region + Management + Time + Condition,
        meta.phy, permutations = 1000)

fun_cap <- capscale(as.matrix(phyt_bray) ~ ., meta.phy[,c(9,2,3)], na.action = na.omit)

anova(fun_cap)
RsquareAdj(fun_cap)

arrowmat <- rbind.data.frame(scores(fun_cap, display = "bp")[-c(6:7),],
                             scores(fun_cap, display = "cn"))
arrowmat$Sample_ID <- row.names(arrowmat)

sitemat <- as.data.frame(scores(fun_cap, display = "sites"))
sitemat$Sample_ID <- row.names(sitemat)
sitemat <- merge(sitemat, meta.phy, by = "Sample_ID")

summary(fun_cap)$cont

ggplot() + 
  geom_point(data = sitemat, aes(x = CAP1, y = CAP2, color = Time), size = 3) +
  geom_segment(data = arrowmat[-c(6:8),], aes(x = 0, y = 0, xend = CAP1*2, yend = CAP2*2), 
               arrow = arrow(length = unit(0.2,"cm")), size = 1, alpha = 0.8) + 
  geom_text(data = arrowmat, aes(x = CAP1*2, y = CAP2*2, label = Sample_ID), size = 5) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black")) +
  xlab("CAP1 (38.43%)") + ylab("CAP2 (13.06%)")

#
#### BETA DIVERSITY: Initial ####
phyt_fun.filt <- subset_samples(phyt_fun, Time == "Initial")
meta.filt <- data.frame(sample_data(phyt_fun.filt))

phyt_bray.filt <- phyloseq::distance(phyt_fun.filt, method = "bray")

nMDS_fun <- metaMDS(phyt_bray.filt)
nMDS_fun$stress

nMDS_fun.plot <- as.data.frame(scores(nMDS_fun))
nMDS_fun.plot$Sample_ID <- rownames(nMDS_fun.plot)
nMDS_fun.plot <- merge(nMDS_fun.plot, meta, by = "Sample_ID")
nMDS_fun.plot$group <- paste(nMDS_fun.plot$Origin, nMDS_fun.plot$Sample_name)
nMDS_fun.plot$order <- as.numeric(nMDS_fun.plot$Time)
nMDS_fun.plot <- nMDS_fun.plot[order(nMDS_fun.plot$order),]

ggplot(nMDS_fun.plot) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Management, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

#
## Beta-dispersion boxplots
bd_fun <- betadisper(phyt_bray.filt, meta.filt$Management, type = "centroid")

bd_fun <- data.frame(Distance_to_centroid = bd_fun$distances)
bd_fun$Sample_ID <- row.names(bd_fun)

bd_fun <- merge(bd_fun, meta, by = "Sample_ID")

#BETADISPER BOXPLOT
ggplot(bd_fun, aes(x = Management, y = Distance_to_centroid, color = Management)) + 
  geom_boxplot(size = 1)

summary(aov(Distance_to_centroid ~ Management, data = bd_fun))

set.seed(1)
adonis(phyt_bray ~ Region*Management, data = meta.ns)



#
##PLOTS CON LAS LINEAS
ggplot(nMDS_fun.plot[order(nMDS_fun.plot$order),]) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Condition, shape = Management, size = Time)) +
  geom_line(aes(x = NMDS1, y = NMDS2, group = group, color = Condition)) +
  scale_size_manual(values = c(5, 3, 2)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

ggplot(subset(nMDS_fun.plot, Time == "Initial")) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Condition, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

#
#### SAVE DATA ####
save.image("physeq_fun.RData")

#