library("phyloseq")
library("ggplot2")
library("ape")
setwd("/media/diego/IYAVVANA/maestria/genomas/Metagenomas/taxonomia")

nat <- import_biom("nat.biom")
##First, it is important to trim the tax table
nat@tax_table@.Data <- substring(nat@tax_table@.Data, 4)
nat@tax_table@.Data[nat@tax_table@.Data==""]<- "Unknown"


rownames(nat.otu) <- paste(#nat@tax_table@.Data[,"Rank1"],nat@tax_table@.Data[,"Rank2"],
  #nat@tax_table@.Data[,"Rank3"],
  #nat@tax_table@.Data[,"Rank4"],
  nat@tax_table@.Data[,"Rank5"],
  nat@tax_table@.Data[,"Rank6"],nat@tax_table@.Data[,"Rank7"], rownames(nat@otu_table@.Data),
  sep = '-')
colnames(nat.otu) <- c("qcuatro","qpocitos-1","qpocitos-2","qarlegundo",
                       "slimones", "scarrizal-1","scarrizal-2")
nat.otu <- otu_table(nat.otu, taxa_are_rows = TRUE)

##Now we need to read the taxonomy information from the biom file
nat.tax <- nat@tax_table@.Data
rownames(nat.tax) <- rownames(nat.otu)
colnames(nat.tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
nat.tax <- tax_table(nat.tax)

##We need to introduce the metadata from the samples
nat.meta <- sample_data(data.frame(
  State= c("Querétaro","Querétaro","Querétaro","Querétaro",
           "SanLuis","SanLuis","SanLuis"),
  Population = c("Cuatro","Pocitos","Pocitos","Arlegundo",
                "Limones","Carrizal","Carrizal"),
  Sample = c("qcuatro","qpocitos-1","qpocitos-2","qarlegundo",
             "slimones", "scarrizal-1","scarrizal-2"),
  row.names=sample_names(nat.otu),
  stringsAsFactors=FALSE
))

##Proceeding to do a tree
nat.tree <- rtree(ntaxa(nat.tax), rooted = TRUE, tip.label=taxa_names(nat.tax))

##Now, we need to merge the data
nat.fin <- merge_phyloseq(nat.otu, nat.tax, nat.tree, nat.meta)
nat.fin@otu_table@.Data = nat.fin@otu_table@.Data[!rownames(nat.fin@otu_table@.Data) %in% "Micrococcales-Microbacteriaceae-Clavibacter-michiganensis", ]
nat.fin@tax_table@.Data = nat.fin@tax_table@.Data[!rownames(nat.fin@tax_table@.Data) %in% "Micrococcales-Microbacteriaceae-Clavibacter-michiganensis", ]
nat.fin@otu_table@.Data = nat.fin@otu_table@.Data[!rownames(nat.fin@otu_table@.Data) %in% "Microbacteriaceae-Clavibacter-michiganensis-28447", ]
nat.fin@tax_table@.Data = nat.fin@tax_table@.Data[!rownames(nat.fin@tax_table@.Data) %in% "Microbacteriaceae-Clavibacter-michiganensis-28447", ]

nat.trim <- filter_taxa(nat.fin, function(x) mean(x) > 50000, TRUE)
nat.trim <- transform_sample_counts(nat.fin, function(x) x / sum(x))

##This next section is to make a PCoA analysis using the different distance methods that are enlisted in the variable 'samp_methods' 
## (jsd, manhattan, euclidean , bray, jaccard, kulczynski)
nat_methods <- c("jsd","manhattan","euclidean","bray","jaccard","kulczynski")
nlist <- vector("list", length(nat_methods))
names(nlist) = nat_methods
for( i in nat_methods ){
  # Calculate distance matrix
  idist <- distance(nat.trim, method=i)
  # Calculate ordination
  iPCoA  <- ordinate(nat.trim, "PCoA", distance=i)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  kr <- NULL
  # Create plot, store as temp variable, p
  kr <- plot_ordination(nat.trim, iPCoA, color="Population", shape="State")
  # Add title to each plot
  kr <- kr + ggtitle(paste("PCoA using distance method ", i, sep=""))
  # Save the graphic to file.
  nlist[[i]] = kr
}

df = ldply(nlist, function(x) x$data)
names(df)[1] <- "distance"
kr = ggplot(df, aes(Axis.1, Axis.2, color=Population, shape=State))+
  geom_point(size=4, alpha=0.8) + 
  facet_wrap(~distance, scales="free") + 
  ggtitle("PCoA on various distance metrics for Natural microbiome on Kraken db")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=15), legend.text=element_text(size=12), 
        legend.title=element_text(size=15),
        strip.text=element_text(size=15))
kr



nat.ord <- ordinate(nat.trim, method ="PCoA" , distance = "euclidean")
plot_ordination(nat.fin, nat.ord, color="Sample", shape = "State") +
  geom_point(size=5, alpha= 0.9) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=15))+
  ggtitle("PCoA with euclidean distance metrics for Natural microbiome on Kraken db")


plot_heatmap(nat.fin, method = "PCoA", distance = "euclidean",  
             taxa.order = "Genus", sample.order = "State") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))


plot_richness(nat.fin, color = "Population", shape="State",
              measures = c("Chao1","Shannon","Simpson"))+
  geom_point(size=5, alpha=0.8) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=15), legend.text=element_text(size=12), 
        legend.title=element_text(size=15),
        strip.text=element_text(size=15))+
  xlab("Samples")

        

