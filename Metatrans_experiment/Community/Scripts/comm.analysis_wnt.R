#### LIBRARIES ####
library(RColorBrewer)
library(dada2)
library(msa)
library(phangorn)
library(reshape2)
library(ggplot2)

setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/Metatrans_experiment/Community/")
#load("physeq_fun.RData")

#
#### FUNCTIONS ####
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

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
#### DATA PREPROCESSING ####
asv_wnt <- readRDS("Inputs/ASV_wnt.rds")
asv_wnt <- asv_wnt[-1,]
asv.t_wnt <- apply(asv_wnt, 1, function(x) x/sum(x))

tax_wnt <- readRDS("Inputs/tax_wnt.rds")
tax_wnt <- cbind.data.frame(tax_wnt, Id = row.names(tax_wnt))
tax_wnt <- gsub("^[a-z]__", "", as.matrix(tax_wnt))

tax_wnt <- set_unid(tax_wnt)

sample_df <- read.table("Inputs/sample_df.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Sample_ID
sample_df <- sample_df[sample_df$Sample_ID %in% row.names(asv_wnt),]

#
#### Phylogenetic tree ####
seqs_wnt <- dada2::getSequences(colnames(asv_wnt))
names(seqs_wnt) <- seqs_wnt

mult_wnt <- msa(seqs_wnt, method = "ClustalW", type = "dna", order = "input")

phang.align_wnt <- as.phyDat(mult_wnt, type = "DNA", names = seqs_wnt)
dm_wnt <- dist.ml(phang.align_wnt)
treeNJ_wnt <- NJ(dm_wnt) # Note, tip order != sequence order
fit_wnt <- pml(treeNJ_wnt, data = phang.align_wnt)

fitGTR_wnt <- update(fit_wnt, k = 4, inv = 0.2)
fitGTR_wnt <- optim.pml(fitGTR_wnt, model = "GTR", optInv = TRUE, optGamma = TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

#
#### TAXONOMIC EXPLORATION ####
##Genus
dataG_wnt <- merge(asv.t_wnt, tax_wnt, by = "row.names")[,c(2:60,66)]
dataG_wnt <- aggregate(dataG_wnt[,1:59], list(dataG_wnt$Genus), sum)

dataG_wnt <- melt(dataG_wnt)
colnames(dataG_wnt) <- c("Genus", "Sample_ID", "value")
dataG_wnt$Genus[dataG_wnt$value < 0.05] <- "Other"

glomG_wnt <- tax_glom(phyt_wnt, taxrank = 'Genus', NArm = FALSE)
dataG_wnt <- psmelt(glomG_wnt)
dataG_wnt$Genus <- as.character(dataG_wnt$Genus)
dataG_wnt$Genus[dataG_wnt$Abundance < 0.05] <- "Other"
dataG_wnt$Genus[is.na(dataG_wnt$Genus)] <- "Unidentified"

orderG <- levels(factor(dataG_wnt$Genus))
orderG <- orderG[! orderG %in% "Other"]
orderG <- append(orderG, "Other")

dataG_wnt <- aggregate(dataG_wnt$value, list(dataG_wnt$Sample_ID, dataG_wnt$Genus), sum)
colnames(dataG_wnt) <- c("Sample_ID", "Genus", "value")

dataG_wnt <- merge(dataG_wnt, sample_df, by = "Sample_ID")

set.seed(123)
ggplot(dataG_wnt, 
       aes(x = Sample_ID, y = value, fill = factor(Genus, levels = orderG))) + 
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
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3))


#
#### ALPHA DIVERISTY: Dynamics ####
alpha_wnt <- cbind.data.frame(# Hill based taxonomic alpha diversity
                              t.q0 = hill_taxa(tasv_wnt, q = 0),
                              t.q1 = hill_taxa(tasv_wnt, q = 1),
                              t.q2 = hill_taxa(tasv_wnt, q = 2),
                              
                              # Hill based phylogenetic alpha diversity
                              p.q0 = hill_phylo(tasv_wnt, fitGTR_wnt$tree, q = 0),
                              p.q1 = hill_phylo(tasv_wnt, fitGTR_wnt$tree, q = 1),
                              p.q2 = hill_phylo(tasv_wnt, fitGTR_wnt$tree, q = 2))

alpha_wnt$Sample_ID <- row.names(alpha_wnt)

alpha_wnt <- merge(alpha_wnt, subset(meta, order > 54)[,-c(7,8)], by = "Sample_ID")

alpha_wnt.plot <- melt(alpha_wnt)
alpha_wnt.plot$diversity <- ifelse(startsWith(as.character(alpha_wnt.plot$variable), "t."),
                                   "Taxonomic", "Phylogenetic")
alpha_wnt.plot$diversity <- factor(alpha_wnt.plot$diversity, 
                                   levels = c("Taxonomic", "Phylogenetic"))

alpha_wnt.plot$q <- as.factor(substr(alpha_wnt.plot$variable, 
                                     nchar(as.character(alpha_wnt.plot$variable)),
                                     nchar(as.character(alpha_wnt.plot$variable))))

alpha_wnt.plot$name <- paste(alpha_wnt.plot$Origin, alpha_wnt.plot$Sample_name, 
                             alpha_wnt.plot$q)

xx <- subset(alpha_wnt.plot, diversity == "Taxonomic")
xx <- cbind.data.frame(xx, 
                       value2 = rbind.data.frame(cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value)))
xxx <- cbind.data.frame(xx[,-14], 
                        value2 = rbind.data.frame(cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value)))


ggplot(subset(alpha_wnt.plot, diversity == "Taxonomic"), 
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
alpha_wnt.norm <- alpha_wnt[,-c(5,6,7)]

alpha_wnt.norm$t.q0 <- ave(alpha_wnt.norm$t.q0, as.factor(alpha_wnt.norm$Origin), 
                           FUN = scale)
alpha_wnt.norm$t.q1 <- ave(alpha_wnt.norm$t.q1, as.factor(alpha_wnt.norm$Origin), 
                           FUN = scale)
alpha_wnt.norm$t.q2 <- ave(alpha_wnt.norm$t.q2, as.factor(alpha_wnt.norm$Origin), 
                           FUN = scale)

alpha_wnt.norm_plot <- melt(alpha_wnt.norm)

ggplot(alpha_wnt.norm_plot, 
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
ggplot(alpha_wnt.norm_plot, 
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
alpha_wnt.plot <- subset(alpha_wnt.plot, diversity == "Taxonomic")

alpha_wnt.plot$name2 <- paste(alpha_wnt.plot$Origin, alpha_wnt.plot$Sample_name,
                              alpha_wnt.plot$variable)

alpha_wnt.dyn <- subset(alpha_wnt.plot, Time == "Initial")
colnames(alpha_wnt.dyn)[10] <- "Initial"

alpha_wnt.dyn <- merge(alpha_wnt.dyn, subset(alpha_wnt.plot, Time == "Middle")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_wnt.dyn)[15] <- "Middle"

alpha_wnt.dyn <- merge(alpha_wnt.dyn, subset(alpha_wnt.plot, Time == "Final")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_wnt.dyn)[16] <- "Final"

alpha_wnt.dyn$`Initial-Final` <- (alpha_wnt.dyn$Initial-alpha_wnt.dyn$Final)/alpha_wnt.dyn$Initial
alpha_wnt.dyn$`Initial-Middle` <- (alpha_wnt.dyn$Initial-alpha_wnt.dyn$Middle)/alpha_wnt.dyn$Initial

alpha_wnt.dyn <- alpha_wnt.dyn[,c(3,4,6,8,10,12,13,17,18)]
alpha_wnt.dyn.plot <- melt(alpha_wnt.dyn)
colnames(alpha_wnt.dyn.plot)[5] <- "div.index"

alpha_wnt.dyn.plot$variable <- factor(alpha_wnt.dyn.plot$variable, 
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


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q0"),
                       cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q1"),
                       cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### ALPHA DIVERISTY: Diversity loss: Normalized ####
alpha_wnt.norm_plot$name <- paste(alpha_wnt.norm_plot$Origin, alpha_wnt.norm_plot$Sample_name,
                                  alpha_wnt.norm_plot$variable)

alpha_wnt.dyn_norm <- subset(alpha_wnt.norm_plot, Time == "Initial")
colnames(alpha_wnt.dyn_norm)[10] <- "Initial"

alpha_wnt.dyn_norm <- merge(alpha_wnt.dyn_norm, 
                            subset(alpha_wnt.norm_plot, Time == "Middle")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_wnt.dyn_norm)[12] <- "Middle"

alpha_wnt.dyn_norm <- merge(alpha_wnt.dyn_norm, 
                            subset(alpha_wnt.norm_plot, Time == "Final")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_wnt.dyn_norm)[13] <- "Final"

alpha_wnt.dyn_norm$`Initial-Final` <- (alpha_wnt.dyn_norm$Initial-alpha_wnt.dyn_norm$Final)
alpha_wnt.dyn_norm$`Initial-Middle` <- (alpha_wnt.dyn_norm$Initial-alpha_wnt.dyn_norm$Middle)

alpha_wnt.dyn_norm <- alpha_wnt.dyn_norm[,c(1:10,14,15)]
alpha_wnt.dyn_norm.plot <- melt(alpha_wnt.dyn_norm)
colnames(alpha_wnt.dyn_norm.plot)[10] <- "div.index"

alpha_wnt.dyn_norm.plot$variable <- factor(alpha_wnt.dyn_norm.plot$variable, 
                                      levels = c("Initial-Middle", "Initial-Final"))

ggplot(alpha_wnt.dyn_norm.plot, 
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


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q0"),
                        cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q1"),
                        cbind.data.frame(subset(alpha_wnt.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### BETA DIVERSITY ####
phyt_bray <- phyloseq::distance(phyt_wnt, method = "bray")
nMDS_wnt <- metaMDS(phyt_bray)
nMDS_wnt$stress

nMDS_wnt.plot <- as.data.frame(scores(nMDS_wnt))
nMDS_wnt.plot$Sample_ID <- rownames(nMDS_wnt.plot)
nMDS_wnt.plot <- merge(nMDS_wnt.plot, meta, by = "Sample_ID")
nMDS_wnt.plot$group <- paste(nMDS_wnt.plot$Origin, nMDS_wnt.plot$Sample_name)
nMDS_wnt.plot$order <- as.numeric(nMDS_wnt.plot$Time)
nMDS_wnt.plot <- nMDS_wnt.plot[order(nMDS_wnt.plot$order),]

ggplot(nMDS_wnt.plot) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Time, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


pcoa_wnt <- wcmdscale(phyt_bray)

pcoa_wnt.plot <- as.data.frame(scores(pcoa_wnt))
pcoa_wnt.plot$Sample_ID <- rownames(pcoa_wnt.plot)
pcoa_wnt.plot <- merge(pcoa_wnt.plot, meta, by = "Sample_ID")
pcoa_wnt.plot$group <- paste(pcoa_wnt.plot$Origin, pcoa_wnt.plot$Sample_name)
pcoa_wnt.plot$order <- as.numeric(pcoa_wnt.plot$Time)
pcoa_wnt.plot <- pcoa_wnt.plot[order(pcoa_wnt.plot$order),]

ggplot(pcoa_wnt.plot) + 
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

meta.ns <- meta.ns[match(sample_names(phyt_wnt), meta.ns$Sample_ID),]

meta.phy <- data.frame(sample_data(phy_wnt))


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
phyt_wnt.filt <- subset_samples(phyt_wnt, Time == "Initial")
meta.filt <- data.frame(sample_data(phyt_wnt.filt))

phyt_bray.filt <- phyloseq::distance(phyt_wnt.filt, method = "bray")

nMDS_wnt <- metaMDS(phyt_bray.filt)
nMDS_wnt$stress

nMDS_wnt.plot <- as.data.frame(scores(nMDS_wnt))
nMDS_wnt.plot$Sample_ID <- rownames(nMDS_wnt.plot)
nMDS_wnt.plot <- merge(nMDS_wnt.plot, meta, by = "Sample_ID")
nMDS_wnt.plot$group <- paste(nMDS_wnt.plot$Origin, nMDS_wnt.plot$Sample_name)
nMDS_wnt.plot$order <- as.numeric(nMDS_wnt.plot$Time)
nMDS_wnt.plot <- nMDS_wnt.plot[order(nMDS_wnt.plot$order),]

ggplot(nMDS_wnt.plot) + 
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
bd_wnt <- betadisper(phyt_bray.filt, meta.filt$Management, type = "centroid")

bd_wnt <- data.frame(Distance_to_centroid = bd_wnt$distances)
bd_wnt$Sample_ID <- row.names(bd_wnt)

bd_wnt <- merge(bd_wnt, meta, by = "Sample_ID")

#BETADISPER BOXPLOT
ggplot(bd_wnt, aes(x = Management, y = Distance_to_centroid, color = Management)) + 
  geom_boxplot(size = 1)

summary(aov(Distance_to_centroid ~ Management, data = bd_wnt))

set.seed(1)
adonis(phyt_bray ~ Region*Management, data = meta.ns)



#
##PLOTS CON LAS LINEAS
ggplot(nMDS_wnt.plot[order(nMDS_wnt.plot$order),]) + 
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

ggplot(subset(nMDS_wnt.plot, Time == "Initial")) + 
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
save.image("comm.analysis_wnt.Rdata")

#