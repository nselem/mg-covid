#This script is useful to manipulate taxonomic assignment databases that have been created with kreken. 
#It is important to create a biom file first with all the *kraken.report files that are desirable to analyze,
# in this example 'kraken-biom' is used, but other malgorithms are available that generate files in json format:

# $ kraken-biom file1-kraken.report file2-kraken.report file3-kraken.report file4-kraken.report --fmt json -o samp.biom

#It is important to change the word "samp" to any name that the user want to gave to its database 



library("phyloseq")
library("ggplot2")
library("ape")
setwd(#dicertory where the files are)

samp <- import_biom("samp.biom")

##First, it is important to trim the tax table
samp@tax_table@.Data <- substring(samp@tax_table@.Data, 4)
##This is to replace any missing assingment a Unknown value
samp@tax_table@.Data[samp@tax_table@.Data==""]<- "Unknown"

##Next, we are going to create a table were the abundances are going to be stores and manipulated
samp.otu <- samp@otu_table@.Data
rownames(samp.otu) <- paste(#samp@tax_table@.Data[,"Rank1"],samp@tax_table@.Data[,"Rank2"],
  #samp@tax_table@.Data[,"Rank3"],
  #samp@tax_table@.Data[,"Rank4"],
  samp@tax_table@.Data[,"Rank5"],
  samp@tax_table@.Data[,"Rank6"],samp@tax_table@.Data[,"Rank7"], rownames(samp@otu_table@.Data),
  sep = '-')
colnames(samp.otu) <- c(#names of the samples that are represented as columns)
samp.otu <- otu_table(samp.otu, taxa_are_rows = TRUE)

##Now we need to read the taxonomy information from the biom file
samp.tax <- samp@tax_table@.Data
rownames(samp.tax) <- rownames(samp.otu)
colnames(samp.tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
samp.tax <- tax_table(samp.tax)

##The next step is to write the metadata of the samples. Here 3 categories are stablished, but this depends on the metadata the user is 
## managing
samp.meta <- sample_data(data.frame(
  Category1 <- c(#Category1 information for each sample),
  Category2 <- c(#Category2 information for each sample),
  Category3 <- c(#Category3 information for each sample),
  row.names=sample_names(samp.otu),
  stringsAsFactors=FALSE
))

##Proceeding to make a tree on the data with ape package
samp.tree <- rtree(ntaxa(samp.tax), rooted = TRUE, tip.label=taxa_names(samp.tax))

##The data needs to be merged:
samp.fin <- merge_phyloseq(samp.otu, samp.tax, samp.tree, samp.meta)
#----------------------------------------------------------------------------------------------------------------------------------------
##If the user needs to trim or manipulate the data in some way, it is possible by means of different commands that can be found
## in the phyloseq information page(https://joey711.github.io/phyloseq/). In this example, the first line is to just leave those OTUs
## that in average have an abundance greater than 50000; the second one is to transform the data into an relative abundance.
samp.trim <- filter_taxa(samp.fin, function(x) mean(x) > 50000, TRUE)
samp.trim <- transform_sample_counts(samp.fin, function(x) x / sum(x))
#----------------------------------------------------------------------------------------------------------------------------------------

##This next section is to make a PCoA analysis using the different distance methods that are enlisted in the variable 'samp_methods' 
## (jsd, manhattan, euclidean , bray, jaccard, kulczynski)
samp_methods <- c("jsd","manhattan","euclidean","bray","jaccard","kulczynski")
nlist <- vector("list", length(samp_methods))
names(nlist) = samp_methods
for( i in samp_methods ){
  ##Calculate distance matrix
  idist <- distance(samp.trim, method=i)
  ##Calculate ordisampion
  iPCoA  <- ordisampe(samp.trim, "PCoA", distance=i)
  ##Make plot
  ##Don't carry over previous plot (if error, p will be blank)
  graf <- NULL
  ##Create plot, store as temp variable, p
  graf <- plot_ordisampion(samp.trim, iPCoA, color="Population", shape="State")
  ##Add title to each plot
  graf <- graf + ggtitle(paste("PCoA using distance method ", i, sep=""))
  ##Save the graphic to file.
  nlist[[i]] = graf
}

df = ldply(nlist, function(x) x$data)
names(df)[1] <- "distance"
graf = ggplot(df, aes(Axis.1, Axis.2, color=Population, shape=State))+
  geom_point(size=4, alpha=0.8) + 
  facet_wrap(~distance, scales="free") + 
  ggtitle("PCoA on various distance metrics for sampural microbiome on Kraken db")+
  #The next parameters are esthetical only and are described in the ggplot package
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=15), legend.text=element_text(size=12), 
        legend.title=element_text(size=15),
        strip.text=element_text(size=15))
graf

##If the user wants to generate only one plot, the next command can be used:
samp.ord <- ordisampe(samp.trim, method ="PCoA" , distance = "euclidean")
plot_ordisampion(samp.fin, samp.ord) +
  geom_point(size=5, alpha= 0.9) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=15))+
  ggtitle("PCoA with euclidean distance metrics for sampural microbiome on Kraken db")

##To generate a heatmap on the samples use the next series of lines:
plot_heatmap(samp.fin, method = "PCoA", distance = "euclidean",  
             taxa.order = "Genus") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))

##Alpha diversity measures can gave useful information concerning the community diversity. Here, we plot the Chao1, Shannon and Simpson
## diversity measures, reasons for this can be found in bibliography.
plot_richness(samp.fin,
              measures = c("Chao1","Shannon","Simpson"))+
  geom_point(size=5, alpha=0.8) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=15), legend.text=element_text(size=12), 
        legend.title=element_text(size=15),
        strip.text=element_text(size=15))+
  xlab("Samples")
