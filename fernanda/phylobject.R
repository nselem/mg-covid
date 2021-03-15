##Metagenomics R script
#by Nelly Selem, Diego Garfias and 
#The code/lesson was retrieved from 
#https://carpentries-incubator.github.io/metagenomics/09-diversity-analysis/index.html

#Libraries
library("ggplot2")
library("readr")
library("phyloseq")
library("patchwork")
library("vegan")


#Avoid replacing mannually
#The code wrote in the next paragraph only is useful in the case of how my data is named
samnum <- "32"
path <- paste(c("/home/mfcg/Descargas/covid/microbiome/AC1ME2SS", samnum, "_S", samnum, "/KRAKEN/results"), collapse = "")
rankedwc <- paste(c("AC1ME2SS", samnum, "_kraken.ranked-wc"), collapse = "")
ltablewc <- paste(c("AC1ME2SS", samnum, "_kraken.lineage_table-wc"), collapse = "")
psfile <- paste(c("AC1ME2SS", samnum, ".rds"), collapse = "")  

#Path
setwd(path)

#Reading data
OTUS <- read_delim(rankedwc,"\t", escape_double = FALSE, trim_ws = TRUE)
TAXAS <- read_delim(ltablewc, "\t", escape_double = FALSE,
                    col_types = cols(subspecies = col_character(),
                                     subspecies_2 = col_character()), trim_ws = TRUE)

abundance <- as.matrix(OTUS[ , -1]) #To avoid that the OTU column is taken as a sample the first column is omitted
lineages <- as.matrix(TAXAS)

row.names(abundance) <- OTUS$OTU #Their identity is not lost cause their names are saved as rownames before
row.names(lineages) <- TAXAS$OTU

OTU <- otu_table(abundance, taxa_are_rows = TRUE)
TAX <- tax_table(lineages)

psmetagenome <- phyloseq(OTU, TAX)
#psmetagenome <- prune_taxa(taxa_sums(psmetagenome)>10, psmetagenome) #filters OTUs with less than 10 reads, avoid if you want untrimmed data
write_rds(psmetagenome, file = psfile)

#Merging

setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq") #set path

#The following steps are to avoid automatic replacing of names
#This steps generate a vector with the names for the samples 1 to 9
acme0 <- rep("AC1ME2SS0",9) #adding 0 in the prefix of the name cause this is how it is saved
unums <- as.character(c(1:9)) #numbers as characters
punto <- ".rds" #the extension
unombres <- paste(acme0, unums, rep(punto, 9), sep = "") #vector saving the names for samples 1 to 9
initial_mgps <- read_rds(unombres[1]) #A initial phyloseq object to start merging
#This steps generate a vector with the names of the following samples 10 to 33, but are useful for bigger numbers, at least up to 99
acme <- rep("AC1ME2SS", 23) #prefix, can be replaced for any prefix, just repeat it up to the number of samples you have
nums <- as.character(c(10:32)) #the numbers as characters, if you have a prefix and it is followed by numeration, you can use this
nombres <- paste(acme, nums, rep(punto, 23), sep = "") #vector saving the file names, it pastes the aforementioned vectors in order, repeat the extension for as many samples you have
#These steps merge all the phyloseq objects with unit names
initial_mgps <- read_rds(unombres[1]) #A initial phyloseq object to start merging
merged_mg_units <- initial_mgps #we need a initial phyloseq object to start merging
for(i in 2:length(unombres)){ #here we start in number 2 cause we already have the phyloseq object of the first element in unombres
  pssolo <- read_rds(unombres[i])
  merged_mg_units <- merge_phyloseq(merged_mg_units, pssolo)
}
#These steps merge the phyloseq objects with unit names with the rest of the files
merged_all_mg <- merged_mg_units #we need a starting phyloseq object, to acelerate the process we will use the merged phyloseq object of the units
for(i in 1:length(nombres)){ #here we start in one cause none of the elements have their phyloseq number yet
  pssolo <- read_rds(nombres[i])
  merged_all_mg <- merge_phyloseq(merged_all_mg, pssolo)
}
merged_all_mg <- subset_taxa(merged_all_mg, superkingdom == "Bacteria") #Only Bacteria OTUs are left, this is done in the merged file
write_rds(merged_all_mg, file = "AC1MESS_untrim.rds")#saving the phyloseq object with all the samples as RDS

#Reading our merged file
merged_all_mg <- read_rds("AC1MESS_untrim.rds")
abs_count_plot <- plot_bar(merged_all_mg, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat = "identity", position = "stack") +
  ggtitle("Absolute abundance")
percentages <- transform_sample_counts(merged_all_mg, function(x) x*100 / sum (x))
perc_plot <- plot_bar(percentages, fill = "phylum") +
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack")

#Retrieving the tables from a phyloseq object and saving them as a csv file
OTUtable <- as.data.frame(otu_table(merged_all_mg))
TAXtable <- as.data.frame(tax_table(merged_all_mg))

#Reading the tables 
SAMtable <- read.csv("metadata-covid3.csv", header = TRUE, row.names = 1)
OTUtable <- read.csv("covid-otu-table.csv", header = TRUE, row.names = 1)
TAXtable <- read.csv("covid-taxa-table.csv", header = TRUE, row.names = 1)

#to save any changes use the following commands
#write.csv(OTUtable, "covid-otu-table.csv")
#write.csv(TAXtable, "covid-taxa-table.csv")
#write.csv(SAMtable, "metadata-covid3.csv")

SAMtable<- sample_data(SAMtable)
OTUtable<- otu_table(OTUtable, taxa_are_rows = TRUE)
TAXtable<- tax_table(as.matrix(TAXtable))
covid_ps <- phyloseq(SAMtable, OTUtable, TAXtable) #integrating all the tables
covid_ps <- read_rds("covid_ups.rds")
#write_rds(covid_ps, file = "covid_ups.rds") #saves phyloseq object

#Alpha diversity plots
plot_richness(covid_ps, x = "Saliva", measures = "Simpson") +
  geom_boxplot()
plot_richness(covid_ps, x = "Family", color = "Symptoms", measures = "Simpson")
plot_richness(covid_ps, x = "Symptoms", measures = "Simpson") +  
  geom_boxplot()

#Beta diversity
#This needs two things: 
#distances, those can vary due to the method used to calculate, the distances can be calculated directly in the ordination function, but they are required for the PERMANOVA
jc_dist <- phyloseq::distance(covid_ps, method = "jaccard") #Jaccard is a presence/absence based distance
bc_dist <- phyloseq::distance(covid_ps, method = "bray") #Bray-Curtis is an abundance based distance
#jc_coa_ord <- phyloseq::distance(covid_ps, method = "unifrac") #Unifrac unweighted is a phylogeny based distance, it requires a tree in the phyloseq object
#and an ordination method, the method here is PCoA, Principal Coordinates Analysis
jc_coa_ord <- phyloseq::ordinate(covid_ps, method = "PCoA", distance = "jaccard") 
bc_coa_ord <- phyloseq::ordinate(covid_ps, method = "PCoA", distance = "bray")
#uf_coa_ord <- phyloseq::ordinate(covid_ps, method = "DPCoA", distance = "unifrac")  #DPCoA can also be used but requires a tree in the phyloseq file
#PCoA plots grouped by saliva PCR results
jc_saliva <- plot_ordination(covid_ps, jc_coa_ord, color = "Saliva") + stat_ellipse()
bc_saliva <- plot_ordination(covid_ps, bc_coa_ord, color = "Saliva") + stat_ellipse()
#PCoA plots grouped by the symptoms presence absence
jc_symptoms <- plot_ordination(covid_ps, jc_coa_ord, color = "Symptoms") + stat_ellipse()
bc_symptoms <- plot_ordination(covid_ps, bc_coa_ord, color = "Symptoms") + stat_ellipse()
#PCoa plots grouped by family
jc_fam <- plot_ordination(covid_ps, jc_coa_ord, color = "Family") + stat_ellipse()
bc_fam <- plot_ordination(covid_ps, bc_coa_ord, color = "Family") + stat_ellipse()
#Combining info in PCoA
jc_fs <- plot_ordination(covid_ps, jc_coa_ord, color = "Family", shape = "Symptoms")
bc_fs <- plot_ordination(covid_ps, bc_coa_ord, color = "Family", shape = "Symptoms")


#PERMANOVA, Permutational multivariate analysis of variance 
adonis(jc_dist ~ sample_data(covid_ps)$Saliva) #Permanova using jaccard distances for Saliva results
adonis(bc_dist ~ sample_data(covid_ps)$Saliva) #Permanova using Bray-Curtis distances for Saliva results
adonis(jc_dist ~ sample_data(covid_ps)$Symptoms) #Permanova using jaccard distances for the symptoms presence
adonis(bc_dist ~ sample_data(covid_ps)$Symptoms) #Permanova using Bray-Curtis distances for the symptoms presence
adonis(jc_dist ~ sample_data(covid_ps)$Family) #Analysis on Family using Jaccard distances
adonis(bc_dist ~ sample_data(covid_ps)$Family) #Analysis on Family using Bray-Curtis distances

