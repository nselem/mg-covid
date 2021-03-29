###TRANSFORM TO PHYLOSEQ AND MERGE###

#The following script shows how to transform to phyloseq objects, 
#the tables produced by the abundance.sh script
#The abundace.sh script can be found here https://github.com/carpentries-incubator/metagenomics/blob/gh-pages/files/abundance.sh

#This script is a modified version from the one found here https://carpentries-incubator.github.io/metagenomics/08-automating_abundance/index.html
#if you wish to learn more, please visit the link above

#This script also shows how to merge the phyloseq objectw produced by different samples, into one.



#LIBRARIES
library("readr")
library("phyloseq")



#PATH AND PREFIX TO THE SAMPLE TABLES
sampath <- "PATH/TO/YOUR/FILES" #write the path to the directory that contains the files produced by the abundance.sh script
setwd(sampath) #set the path as the working directory
sampname <- "sample_kraken" #write the prefix of your files, or the name that is before every extension, it should look similar to the example

#FROM ABUNDANCE.SH TABLES TO PHYLOSEQ OBJECT

psmetagenome <- list() #make an empty list to save the phyloseq object that will be produced by the function ab_to_ps

#The following function automates the steps to transform from the tables done by abundace.sh to a phyloseq object
ab_to_ps <- function(sampname){ #this function uses the sampname object to work, sampname = sampname
  #File names
  rankedwc <- paste0(sampname, ".ranked-wc") #add the name of your ranked-wc file 
  ltablewc <- paste0(sampname, ".lineage_table-wc") #add the name of your lineage_table-wc file 
  
  #Read tables
  OTUS <- read_delim(rankedwc,"\t", escape_double = FALSE, trim_ws = TRUE)
  TAXAS <- read_delim(ltablewc, "\t", escape_double = FALSE,
                      col_types = cols(subspecies = col_character(),
                                       subspecies_2 = col_character()), trim_ws = TRUE)
  #Matrix format
  abundance <- as.matrix(OTUS[ , -1]) #To avoid that the OTU column is taken as a sample the first column is omitted
  lineages <- as.matrix(TAXAS[ , -1]) #To avoid that the OTU column is repeated the first column is omitted
  
  #Retrieving the OTUs' identity as rownames using the the original file as reference
  row.names(abundance) <- OTUS$OTU 
  row.names(lineages) <- TAXAS$OTU
  
  #Phyloseq format
  OTU <- otu_table(abundance, taxa_are_rows = TRUE)
  TAX <- tax_table(lineages)
  psmetagenome <<- phyloseq(OTU, TAX)
  
}

#Save phyloseq object to a rds file
write_rds(psmetagenome, file = "filename.rds") 

#Below are some modifications that can be done to your phyloseq object
#Although, these modifications were not done individually to each sample phyloseq object
#Neither at this point of the analysis

#Prune data
#Prune taxa with reads below a treshold 
#Not recommended before an alpha diversity test
#prunt_meta <- prune_taxa(taxa_sums(psmetagenome)>10, psmetagenome) #filters OTUs with less than 10 reads, avoid if you want untrimmed data

#Subset to certain taxa
#Subset data to a certain taxonomic group
#Not necessary this will be done in the merged file
#subs_meta <- subset_taxa(psmetagenome, superkingdom == "Bacteria")  #Only Bacteria OTUs are left



#PHYLOSEQ FILES FROM ALL THE SAMPLES
#To make the work easier move all the phyloseq objects that you will merge to the same directory
#The directory must be exclusive for the sample phyloseq objects
psfpath <- "/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq" #write the path to the directory that contain the phyloseq objects
setwd(psfpath)
all_samp_ps <- list.files(path = ".", pattern = ".rds") #make a character vector with the names of the phyloseq objects
merged_all_ps <- list()
for(i in 1:length(all_samp_ps)){
  nps <- read_rds(all_samp_ps[i])
  merged_all_ps <<- merge_phyloseq(merged_all_ps, nps)
}

#write_rds(merged_all_ps, file="all_samp_ps.rds") #The merged phyloseq object can be saved as it is

#Phyloseq object that only contain bacteria OTUs
bacteria_pss <- subset_taxa(merged_all_ps, superkingdom == "Bacteria")  #Only Bacteria OTUs are left in the merged file
write_rds(bacteria_pss, file = "metagenomes.rds") #Save the phyloseq object that contains the Bacteria OTUs from all the samples

#RETRIEVE THE TABLES

#From phyloseq to matrix
OTUtable <- as.data.frame(otu_table(bacteria_pss))
TAXtable <- as.data.frame(tax_table(bacteria_pss))

#Save as a csv file
write.csv(OTUtable, "otu-table.csv")
write.csv(TAXtable, "taxa-table.csv")



##ADD SAMPLE DATA

#Read a csv file containing a table with the sample data or metadata
SAMtable <- read.csv("metada.csv", header = TRUE, row.names = 1)

#Ensure that the sample names are the same in the sample data table as in the OTU table
#Change the name of the OTUtable columns, that correspond to the sample names, they are already sorted due to the list.files function
#Sort in the same manner, alphabetically, the SAMtable row or sample names, and use the to replace the names in the OTUtable
colnames(OTUtable) <- sort(rownames(SAMtable), decreasing = FALSE)

#Tables to phyloseq format
TAXtable <- tax_table(as.matrix(TAXtable))
OTUtable <- otu_table(as.matrix(OTUtable), taxa_are_rows = TRUE)
SAMtable <- sample_data(SAMtable)
bacteria_pss <- phyloseq(TAXtable, OTUtable, SAMtable) #Add the metadata from the SAMtable and update the names from the OTUtable in the new bacteria pss

#Save changes
write_rds(bacteria_pss, file = "metagenomes.rds") #This file is ready to be used in the analysis done in phylobject.R
