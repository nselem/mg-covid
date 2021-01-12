library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
library("plyr")

#------------------------------------------------------------------------
# Datos
samples<-c("AC1ME2SS01_kraken","AC1ME2SS02_kraken","AC1ME2SS03_kraken")
#AC1ME2SS29_kraken.ranked-wc  
#AC1ME2SS42_kraken.lineage_table-wc 
names(samples)<-c("M1","M2","M3")

#-------------------------------------------------------------------------
# Function read data7
read_data <- function(file) {
  ranked_file<-paste("table/",file,".ranked-wc", sep = "")
  lineaged_file<-paste("table/",file,".lineage_table-wc", sep = "")
  
  #print(ranked_file)
  #print(lineaged_file)
  
  OTUS<-read_delim(ranked_file,"\t", escape_double = FALSE, trim_ws = TRUE)
  TAXAS<-read_delim(lineaged_file, "\t", 
                    escape_double = FALSE,  col_types = cols(subspecies = col_character(), 
                                                             subspecies_2 = col_character()), trim_ws = TRUE)
  
  names1 = OTUS$OTU
  names2 = TAXAS$OTU
  
  OTUS$OTU = NULL
  TAXAS$OTU = NULL
  
  abundances = as.matrix(OTUS)
  lineages = as.matrix(TAXAS)
  
  row.names(abundances) = names1
  row.names(lineages) = names2
  
  OTU = otu_table(abundances, taxa_are_rows = TRUE)
  TAX = tax_table(lineages)
  
  metagenome = phyloseq(OTU, TAX)
  Bacteria <- subset_taxa(metagenome, superkingdom == "Bacteria")
  
  metagenome <- prune_taxa(taxa_sums(metagenome)>10,metagenome)
  
  #no_contam <- subset_taxa(metagenomes, family != "mitochondria" & class != "Chloroplast" & genus != "Escherichia" & genus != "Staphylococcus", genus != "Wolbachia") 
  no_nullo <- subset_taxa(metagenome, phylum != "NA") 
  metagenome <- prune_taxa(c(taxa_names(no_nullo)),metagenome)
  
  
  return(metagenome)    
  
}    

#---------------------------------------------------------------------------
samples<-c("AC1ME2SS01_kraken","AC1ME2SS02_kraken","AC1ME2SS03_kraken")
names(samples)<-c("M1","M2","M3")

#samples<-c("pKV1MC6SS01_S3_L002","pKV1MC6SS02_S4_L002","pKV1MC6SS03_S6_L003")
#names(samples)<-c("S3_L002","S4_L002","S6_L003")

#metagenomes<-apply(as.matrix(samples),1, read_data)
metagenomes<-lapply(samples, read_data)

m1<-metagenomes[[1]]
m2<-metagenomes[[2]]
m3<-metagenomes[[3]]

merged_metagenomes = lapply(metagenomes, merge_phyloseq)
merged_metagenomes1 = merge_phyloseq(m1,m2,m3)
#merge_phyloseq(vector)
p = plot_richness(merged_metagenomes1, measures = c("Observed", "Chao1", "Shannon")) 
p + geom_point(size=5, alpha=0.7)  
#__________________________________________________________________________________________________

percentages  = transform_sample_counts(merged_metagenomes1, function(x) x*100 / sum(x) )
absolute_count = plot_bar(merged_metagenomes1, fill="family")
absolute_count = absolute_count + geom_bar(aes(color=family, fill=family), stat="identity", position="stack") + ggtitle("Absolute abundance")

percentages = plot_bar(percentages, fill="family")
percentages = percentages + geom_bar(aes(color=family, fill=family), stat="identity", position="stack") + ggtitle("Relative abundance")

absolute_count | percentages  
