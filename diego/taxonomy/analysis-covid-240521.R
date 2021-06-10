## Covid analysis of taxonomic assignation data:

#-- Let's load the packages needed for the analysis --#
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
library(vegan)
library("picante")

setwd("C:/Users/cairo/Documents/sideprojects/covid/from-reads")
#-- First, we need to import our data with phyloseq command --#                                                    --#

covid <- import_biom("non-human/covid.biom")
#-- Now, we will get trim our data for the analysis: --#
covid@tax_table@.Data <- substring(covid@tax_table@.Data, 4)
colnames(covid@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#-- Also, we need to load the metadata from our file --#
meta <- read.csv("metadata-covid3.csv", row.names = 1)
colnames(covid@otu_table@.Data) <- row.names(meta)
#-- Let's tell R that the `meta` object will be part of a phyloseq object, and merge them --#
meta <- sample_data(meta)
covid <- merge_phyloseq(covid,meta)

#-- Removing non bacterial OTUs --#

bacteria <- subset_taxa(covid, Kingdom == "Bacteria")
covid <- prune_taxa(taxa_names(bacteria), covid)

no_cont <- subset_taxa(covid, Family != "mitocondria" & Class != "Chloroplast")
covid <- prune_taxa(taxa_names(no_cont), covid)

#-- Analyzing the deepness of the sequenciation --#

deepn <- data.frame(
  samples=as.character(map(strsplit(colnames(covid@otu_table@.Data), "_"),1)),
  reads= sample_sums(covid))
ggplot(data = deepn, aes(y= reads, x = samples))+
  geom_bar(stat="identity", fill="violetred4") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept = 4242198, col= "cyan3", size = 1.5, alpha = 0.5) +
  xlab("Samples") + ylab("Number of reads") +
  ggtitle("Deepness of covid libraries") +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))
#¡# By this result, we see that the sample SS09 is underepresented so we will
#¡# leave it out of the analysis
#-- Redefining the objets to comply with this result: --#
meta <- meta[-13,]
covid.otu <- covid@otu_table@.Data
covid.otu <- as.data.frame(covid.otu)
covid.otu <- select(covid.otu, -13)
colnames(covid.otu) <- row.names(meta)
covid.otu <- otu_table(covid.otu, taxa_are_rows = TRUE)
covid.tax <- covid@tax_table@.Data
covid.tax <- tax_table(covid.tax)
covid <- merge_phyloseq(covid.otu, covid.tax,meta)

#-- Defining the normalization method --#
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

#-- Defining the new phyloseq object --#
covid.otu2 <- otu_table(z@.Data[[1]], taxa_are_rows = TRUE)
#-- Merging all the objects in the new normalized phyloseq object --#
covid2 <- merge_phyloseq(covid.otu2, covid.tax,meta)

#-- Let's transform the read counts in relative abundances --#
covid2 <- transform_sample_counts(physeq = covid2, function(x) x*100/sum(x))

#-- In order to use vegan for the multivariate analysis, we will extract our data from the phyloseq objects --#
d.covid <- t(covid2@otu_table@.Data)

#-- Obtainment of beta diversity with NMDS --#
covid.ord <- ordinate(physeq = covid2,method = "NMDS", distance = "bray")
#-- Let's plot the beta diversity analysis --#
plot_ordination(covid2, covid.ord, color = "Covid") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Covid),
               type = "norm", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  theme_bw()+
  geom_point(size=4, alpha=0.8) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis") 

#-- Getting the data needed for the multivariate analysis --#
meta.covid <- data.frame(Age = as.factor(covid2@sam_data@.Data[[5]]),
                         Sex = as.factor(covid2@sam_data@.Data[[2]]),
                         Pacient = as.factor(covid2@sam_data@.Data[[3]]),
                         Symt = as.factor(covid2@sam_data@.Data[[6]]),
                         Diag = as.factor(covid2@sam_data@.Data[[7]]),
                         Covid = as.factor(covid2@sam_data@.Data[[10]]))
#-- Multivariate analysis --#
adonis(d.covid ~ Diag * Symt, data = meta.covid, permutations = 999, strata = meta.covid$Covid)
adonis(d.covid ~ Diag , data = meta.covid, permutations = 999)
adonis(d.covid ~ Symt , data = meta.covid, permutations = 999)
adonis(d.covid ~ Age , data = meta.covid, permutations = 999)
adonis(d.covid ~ Covid , data = meta.covid, permutations = 999)
## A combined analysis of multivariance ##
adonis(d.covid ~ Age*Pacient, data = meta.covid, permutations = 999)


#N According to the MAGs annotation from Professor Nelly Selem, we have some taxa that made
#N the majority of the read counts, we will extract this Genera of microbes, and do the
#N same analysis as with the entire dataset

#-- Extracting the bacteria at Genus level that where identified at MAG level --#
genera <- c("Streptococcus","Actinomyces","Eubacterium",
            "Gemella","Rothia","Prevotella","Megasphaera")

t.covid <- tax_glom(physeq = covid2, taxrank = "Genus")
#-- Defining the `tax_table` for the phyloseq object --#
tax.covid <- as.data.frame(covid2@tax_table@.Data)
tax.covid <- tax.covid[(tax.covid$Genus %in% genera),]
row.names(tax.covid)
tax.covid <- as.matrix(tax.covid)
tax.covid <- tax_table(tax.covid)
#-- Defining the `otu_table` for the phyloseq object --#
otu.covid <- as.data.frame(covid2@otu_table@.Data)
otu.covid <- otu.covid[row.names(tax.covid),]
otu.covid <- otu_table(otu.covid, taxa_are_rows = TRUE)

#-- Merging all into a trimmed phyloseq object --#
trim.covid <- merge_phyloseq(otu.covid, tax.covid,meta)

#-- Obtainment of beta diversity with NMDS --#
covid.ord <- ordinate(physeq = trim.covid,method = "NMDS", distance = "bray")
#-- Let's plot the new beta diversity analysis --#
plot_ordination(trim.covid, covid.ord, color = "Covid") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Covid),
               type = "t", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  theme_bw()+
  geom_point(size=4, alpha=0.8) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis") 

#-- Extraction of the data from the phyloseq object for the multivariate analysis -##
t.covid <- t(trim.covid@otu_table@.Data)

#-- Multivariate analysis --#
adonis(t.covid ~ Diag * Symt, data = meta.covid, permutations = 999)
adonis(t.covid ~ Diag , data = meta.covid, permutations = 999)
adonis(t.covid ~ Symt , data = meta.covid, permutations = 999)
adonis(t.covid ~ Age , data = meta.covid, permutations = 999)
adonis(t.covid ~ Covid , data = meta.covid, permutations = 999)
## A combined analysis of multivariance ##
adonis(t.covid ~ Covid*Diag, data = meta.covid, permutations = 999)
