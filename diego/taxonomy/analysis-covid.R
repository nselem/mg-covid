library("phyloseq")
library("ggplot2")
library("edgeR")
library("DESeq2")
library("pheatmap")
library("readr")
library("tidyr")
library("purrr")
library("ape")
library("plyr")
library("dplyr")
library(RColorBrewer)
library(ggpubr)

#setwd("C:/Users/cairo/Documents/sideprojects/covid/")
setwd("D:/documentos/sideprojects/covid/from-reads")

### ---------- Lets begin with the assembly of the phyloseq object ---------- ###

#Lets import the biom object where the data is stored
covid <- import_biom("covid-taxon.biom")
covid@tax_table@.Data <- substring(covid@tax_table@.Data, 4)
# Defining the tax table where the OTU names information is hoarded
covid.tax <- covid@tax_table@.Data
rownames(covid.tax) <- rownames(covid@otu_table@.Data)
colnames(covid.tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
covid.tax <- tax_table(covid.tax)
# Load the csv file where the metadata is located
meta <- read.csv("metadata-covid3.csv", row.names = 1)
## Get rid of the library that is underrepresented
meta <- meta[-13,]
meta <- sample_data(meta)
# Load the out table where the OTU abundance in each sample is stored
covid.otu <- covid@otu_table@.Data
##colnames(covid.otu) <-  as.character(map(strsplit(colnames(covid.otu@.Data), "_"),1))
covid.otu <- as.data.frame(covid.otu)
## Get rid of the library that is underrepresented
covid.otu <- select(covid.otu, -13)
colnames(covid.otu) <- row.names(meta)
covid.otu <- otu_table(covid.otu, taxa_are_rows = TRUE)
# Now, lets do the phylogenetic tree of the samples
#-#covid.tree <- rtree(ntaxa(covid.tax), rooted = TRUE, tip.label=taxa_names(covid.tax))
# Next, lets merge the objects for the data treatment steps
#-#covid <- merge_phyloseq(covid.otu, covid.tax, covid.tree, meta)
covid <- merge_phyloseq(covid.otu, covid.tax,meta)

### -------- Removing non bacterial OTUs --------- ###

bacteria <- subset_taxa(covid, Kingdom == "Bacteria")
covid <- prune_taxa(taxa_names(bacteria), covid)

no_cont <- subset_taxa(covid, Family != "mitocondria" & Class != "Chloroplast")
covid <- prune_taxa(taxa_names(no_cont), covid)

### -------- Analyzing the deepness of the sequenciation --------- ####

deepn <- data.frame(
  samples=as.character(map(strsplit(colnames(covid@otu_table@.Data), "_"),1)),
  reads= sample_sums(covid))
ggplot(data = deepn, aes(y= reads, x = samples))+
  geom_bar(stat="identity", fill="violetred4") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept = 4242198, col= "cyan3", size = 1.5, alpha = 0.5) +
  xlab("Samples") + ylab("Number of reads") +
  ggtitle("Deepness of covid libraries") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))

### --------- Before any normalization, the alpha diversity indexes -------- ###

# In order to generate the plots, we need to generate an object "p" from where the data will be retrieved
p <- plot_richness(covid, color = "Pacient",
                   measures =c( "Observed","Chao1","Shannon","Simpson")) 

simp <- p$data$variable == "Simpson"
simp <- p$data[simp,]

shan <- p$data$variable == "Shannon"
shan <- p$data[shan,]

chao <- p$data$variable == "Chao1"
chao <- p$data[chao, ]

simp$Pacient <- as.factor(simp$Pacient)
shan$Pacient <- as.factor(shan$Pacient)
chao$Pacient <- as.factor(chao$Pacient)

simp$Saliva <- as.factor(simp$Saliva) 
shan$Saliva <- as.factor(shan$Saliva)
chao$Saliva <- as.factor(chao$Saliva)

simp$Symptoms <- as.factor(simp$Symptoms)
shan$Symptoms <- as.factor(shan$Symptoms)
chao$Symptoms <- as.factor(chao$Symptoms)

#By means of the next code we can explore some of the variables in the data 
ggplot(data= simp, aes(x= Symptoms, y= value)) +
  #geom_violin(trim = FALSE) + scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2, aes(color = Symptoms))
#geom_point(size = 2, aes(color = Pacient))

## Lets create the mixed plot ##
g1 <-ggplot(data= chao, aes(x= Pacient, y= value)) +
  #geom_violin(trim = FALSE) + scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 3.5, aes(color = Pacient)) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  labs(title = "Chao1 index", y ="Chao1 diversity index"  )+
  theme_bw()

g2 <-ggplot(data= simp, aes(x= Pacient, y= value)) +
  #geom_violin(trim = FALSE) + scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 3.5, aes(color = Pacient)) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  labs(title = "Simpson index", y ="Simpson diversity index"  )+
  theme_bw()

g3 <-ggplot(data= shan, aes(x= Pacient, y= value)) +
  #geom_violin(trim = FALSE) + scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 3.5, aes(color = Pacient)) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  labs(title = "Shannon index", y ="Shannon diversity index" )+
  theme_bw()

ggarrange(g1 + rremove("legend"), g2 + rremove("legend"), g3 + rremove("legend"),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1) 

### ---------- Plotting the beta-diversity indexes -------- ###

covid.ord <- ordinate(covid, method = "PCoA", distance = "bray")#, weighted= TRUE)
gc1 <- plot_ordination(covid, covid.ord, color = "Pacient") +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "PCoA - Bray-Curtis")

covid.ord <- ordinate(covid, method = "PCoA", distance = "jaccard")#, weighted= TRUE)
gc2 <- plot_ordination(covid, covid.ord, color = "Pacient") +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "PCoA - Jaccard")

covid.ord <- ordinate(covid, method = "PCoA", distance = "euclidean")
gc3 <- plot_ordination(covid, covid.ord, color = "Pacient") +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "PCoA - Euclidean")

#covid.ord <- ordinate(covid, method = "NMDS", distance = "bray", weighted= FALSE)
#gc4 <- plot_ordination(covid, covid.ord, color = "Pacient") +
#  geom_point(size=3, alpha=0.5) +
#  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
#  theme_bw()+
#  labs(title = "PCoA - Unifrac")

gc4 <- barplot(covid.ord$values$Relative_eig)

ggarrange(gc1 + rremove("legend"), gc2 + rremove("legend"),
          gc3 ,
          ncol = 2, nrow = 2) 


### ---------- Defining the normalization method ---------- ###

proportion = function(physeq) {
  # Normalize total sequences represented
  normf = function(x, tot = max(sample_sums(physeq))) {
    tot * x/sum(x)
  }
  physeq = transform_sample_counts(physeq, normf)
  # Scale by dividing each variable by its standard deviation. physeq =
  # transform_sample_counts(physeq, function(x) x/sd(x)) Center by subtracting
  # the median physeq = transform_sample_counts(physeq, function(x)
  # (x-median(x)))
  return(physeq)
}
covid <- proportion(covid)

## !!!!! The normalization method was by variance stabilization with edgeR !!!!! ##

edgeRnorm = function(physeq, ...) {
  require("edgeR")
  require("phyloseq")
  # physeq = simlist[['55000_1e-04']] z0 = simlisttmm[['55000_1e-04']] physeq
  # = simlist[['1000_0.2']] z0 = simlisttmm[['1000_0.2']] Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # See if adding a single observation, 1, everywhere (so not zeros) prevents
  # errors without needing to borrow and modify calcNormFactors (and its
  # dependent functions) It did. This fixed all problems.  Can the 1 be
  # reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  # Don't need the following additional steps, which are also built-in to some
  # of the downstream distance methods. z1 = estimateCommonDisp(z) z2 =
  # estimateTagwiseDisp(z1)
  return(z)
}
    z<- edgeRnorm(covid, method = "TMM")
    z1<-edgeRnorm(covid, method = "RLE")
    z2<-edgeRnorm(covid, method = "upperquartile")

    covid.otu2 <- otu_table(z@.Data[[1]], taxa_are_rows = TRUE)
    covid2 <- merge_phyloseq(covid.otu2, covid.tax,meta)

### ---------- Exploring the data with a heatmap ---------- ###
    
#In order to see if the bacterial OTUs identified in other research projects
covid.top <- filter_taxa(covid2, function(x) mean(x) > 50000, TRUE)

rownames(covid.top@otu_table@.Data) <- paste(#rownames(covid.top@otu_table@.Data),
  covid.top@tax_table@.Data[,"Genus"],
  covid.top@tax_table@.Data[,"Species"],
  sep = '-')
rownames(covid.top@tax_table@.Data) <- paste(#rownames(covid.top@tax_table@.Data),
  covid.top@tax_table@.Data[,"Genus"],
  covid.top@tax_table@.Data[,"Species"], 
  sep = '-')
plot_heatmap(covid2, method = "PCoA", distance = "bray" )

### --------- Creating specific databases for bacteria of interest

## Those that appear when bacterial complications ##
veillo <- subset_taxa(covid2, Genus == "Veillonella")
prevo <- subset_taxa(covid2, Genus == "Prevotella")
kleb <- subset_taxa(covid2, Genus == "Klebsiella")
rothia <- subset_taxa(covid2, Genus == "Rothia")
## Those whose abundance seems to increase
strep <- subset_taxa(covid2, Genus == "Streptococcus")
pseudo <- subset_taxa(covid2, Genus == "Pseudomonas")
bacil <- subset_taxa(covid2, Genus == "Bacillus")
## The combined object ##
g.taxa <- merge_phyloseq(veillo, prevo, kleb, rothia,
                         strep, pseudo, bacil)

rownames(g.taxa@otu_table@.Data) <- paste(#rownames(covid.top@otu_table@.Data),
  g.taxa@tax_table@.Data[,"Genus"],
  g.taxa@tax_table@.Data[,"Species"],
  sep = '-')
rownames(g.taxa@tax_table@.Data) <- paste(#rownames(covid.top@tax_table@.Data),
  g.taxa@tax_table@.Data[,"Genus"],
  g.taxa@tax_table@.Data[,"Species"], 
  sep = '-')

g.taxa.gen <- tax_glom(g.taxa ,taxrank=rank_names(g.taxa)[6])

plot_heatmap(bacil, method = "NMDS", distance = "bray" )


.### ---------- Assembling the data for the pheatmap ---------- ###

covid.trim <- transform_sample_counts(covid, function(x) x / sum(x))
## In order to see the data at phylum level
covid.phyla <- tax_glom(covid.trim,taxrank = rank_names(covid.trim)[2])
covid.phyla.top <- filter_taxa(covid.phyla, function(x) mean(x) > 0.005, TRUE)

covid.phyla <- tax_glom(covid,taxrank = rank_names(covid)[2])
covid.phyla.top <- filter_taxa(covid.phyla, function(x) mean(x) > 100000, TRUE)
## So as to see the data at genus level 
covid.gene <- tax_glom(covid.trim,taxrank = rank_names(covid.trim)[6])
covid.gene.top <- filter_taxa(covid.gene, function(x) mean(x) > 0.05, TRUE)

covid.gene <- tax_glom(covid2,taxrank = rank_names(covid)[6])
covid.gene.top <- filter_taxa(covid.gene, function(x) mean(x) > 50000, TRUE)

## So as to see the dominant OTUs

covid.frame <- as.data.frame(covid.gene.top@otu_table@.Data)
rownames(covid.frame) <- covid.gene.top@tax_table@.Data[,6]

top.log <- log10(covid.frame)
top.log[top.log=="Inf"] <- 0
top.log[top.log=="-Inf"] <- 0

breaksList = seq(3, 7, by = .4)

my.col <- data.frame(Pacient=meta$Pacient, row.names = rownames(meta))
my.col$Age <- meta$Pacient.Age

my.col <- data.frame(Saliva=meta$Saliva, row.names = rownames(meta))
my.col$Family <- meta$Family
my.col$Age <- meta$Pacient.Age

my_color <- list(Saliva = c(`Negative`= "#a6cee3", `Positive`= "#33a02c", 
                            `Negative_Control`= "#1f78b4"),
                 Family = c(`Fam1`="#1b9e77" , `Fam2`="#7570b3", 
                            `Fam3`="#e7298a", `Na`="#bababa"),
                 Age = c(`<18` = "#4d9221", `20-30`= "#a1d76a",
                         `30-40` = "#fde0ef", `40-55`= "#e9a3c9",
                         `>55` = "#c51b7d"))

pheatmap(top.log,
         color = (c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014",
                    "#cc4c02","#993404","#662506")), 
         breaks = breaksList , cluster_cols = TRUE, 
         cutree_cols = 4, cutree_rows = 5, border_color ="#000000",
         annotation_col = my.col, annotation_colors = my_color,)

### --------- In order to do the before/after figure --------- ###

bef.aft <- subset_samples(covid2, Same != "0")
bef.aft <- tax_glom(bef.aft,taxrank = rank_names(covid)[6])
bef.aft <- filter_taxa(bef.aft, function(x) mean(x) > 50000, TRUE)

ba.frame <- as.data.frame(bef.aft@otu_table@.Data)
rownames(ba.frame) <- bef.aft@tax_table@.Data[,6]

top.log <- log10(ba.frame)
top.log[top.log=="Inf"] <- 0
top.log[top.log=="-Inf"] <- 0

breaksList = seq(3, 7, by = .4)

my.col <- data.frame(Pacient=bef.aft@sam_data@.Data[[7]], 
                     row.names = bef.aft@sam_data@row.names)
my.col$Sample <- bef.aft@sam_data@.Data[[4]]
my.col$Time <- c("Before","After","Before","After","Before","After",
                 "Before","After","Before","After")

my_color <- list(Pacient = c(`Negative`= "#a6cee3", `Positive`= "#33a02c"),
                 Sample = c(`1`="#e41a1c" , `2`="#377eb8", 
                            `3`="#4daf4a", `4`="#984ea3", `5`="#ff7f00"),
                 Time = c(`Before` = "#018571", `After`= "#a6611a"))

pheatmap(top.log,
         color = (c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014",
                    "#cc4c02","#993404","#662506")), 
         breaks = breaksList , 
         cluster_cols= FALSE, border_color ="#000000",
         cutree_cols = 4, gaps_col = c(2,4,6,8), cutree_rows = 4,  
         annotation_col = my.col, annotation_colors = my_color,)

### --------- In order to do the species specific figure --------- ###

prevo1 <- subset_taxa(covid2, Genus == "Prevotella" & Species == "intermedia" )
prevo2 <- subset_taxa(covid2, Genus == "Prevotella" & Species == "melaninogenica" )
staph1 <- subset_taxa(covid2, Genus == "Staphylococcus" & Species == "aureus" )
strep1 <- subset_taxa(covid2, Genus == "Streptococcus" & Species == "pneumoniae" )
veillo1 <- subset_taxa(covid2, Genus == "Veillonella" & Species == "parvula" )
actino1 <- subset_taxa(covid2, Genus == "Acinetobacter" & Species == "baumannii")
kleb1 <- subset_taxa(covid2, Genus == "Klebsiella" & Species == "pneumoniae")
lepto1 <- subset_taxa(covid2, Genus == "Klebsiella" & Species == "pneumoniae")
strep1 <- subset_taxa(covid2, Genus == "Streptococcus" & Species == "oralis")
strep2 <- subset_taxa(covid2, Genus == "Streptococcus" & Species == "mitis")
roth1 <- subset_taxa(covid2, Genus == "Rothia" & Species == "mucilaginosa")
roth2 <- subset_taxa(covid2, Genus == "Rothia" & Species == "dentocariosa")
int.abun <- merge_phyloseq(prevo1, prevo2, staph1, strep1, veillo1, actino1,
                           kleb1, lepto1, strep1, strep2, roth1, roth2)

rownames(int.abun@otu_table@.Data) <- paste(int.abun@tax_table@.Data[,"Genus"],
  int.abun@tax_table@.Data[,"Species"],
  sep = '-')
rownames(int.abun@tax_table@.Data) <- paste(int.abun@tax_table@.Data[,"Genus"],
  int.abun@tax_table@.Data[,"Species"], 
  sep = '-')

sel.taxa <- as.data.frame(int.abun@otu_table@.Data)


top.log <- log10(sel.taxa)
top.log[top.log=="Inf"] <- 0
top.log[top.log=="-Inf"] <- 0

breaksList = seq(3, 7, by = .4)

my.col <- data.frame(Saliva=meta$Saliva, row.names = rownames(meta))
my.col$Family <- meta$Family
my.col$Age <- meta$Pacient.Age

my_color <- list(Saliva = c(`Negative`= "#a6cee3", `Positive`= "#33a02c", 
                            `Negative_Control`= "#1f78b4"),
                 Family = c(`Fam1`="#1b9e77" , `Fam2`="#7570b3", 
                            `Fam3`="#e7298a", `Na`="#bababa"),
                 Age = c(`<18` = "#4d9221", `20-30`= "#a1d76a",
                         `30-40` = "#fde0ef", `40-55`= "#e9a3c9",
                         `>55` = "#c51b7d"))

pheatmap(top.log,
         color = (c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014",
                    "#cc4c02","#993404","#662506")), 
         breaks = breaksList , cluster_cols = TRUE, 
         cutree_cols = 3, cutree_rows = 4, border_color ="#000000",
         annotation_col = my.col, annotation_colors = my_color,)
