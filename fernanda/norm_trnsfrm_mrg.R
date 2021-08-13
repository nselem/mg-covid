###TRANSFORM TO PHYLOSEQ AND MERGE###

#The following script shows how to transform to phyloseq objects, 
#the tables produced by the abundance.sh script
#The abundace.sh script can be found here https://github.com/carpentries-incubator/metagenomics/blob/gh-pages/files/abundance.sh

#This script is a modified version from the one found here https://carpentries-incubator.github.io/metagenomics/08-automating_abundance/index.html
#if you wish to learn more, please visit the link above

#This script also shows how to merge the phyloseq object produced by different samples, into one.

#Finally it helps you to prepare your phyloseq object for the following analyses by leaving only bacterial reads, 
#samples with a similar read depth and normalized data

#LIBRARIES
library("readr")
library("phyloseq")
library("edgeR")


#PATH AND PREFIX TO THE SAMPLE TABLES
muestra <- "AC1ME2SS10"
sampath <- paste0("/home/mfcg/Descargas/covid/microbiome/taxonomia/", muestra, "/results") #write the path to the directory that contains the files produced by the abundance.sh script
setwd(sampath) #set the path as the working directory
sampname <- paste0(muestra, "_kraken") #write the prefix of your files, or the name that is before every extension, it should look similar to the example

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

ab_to_ps(sampname = sampname)

#Save phyloseq object to a rds file
filename <- paste0(muestra, ".rds")
write_rds(psmetagenome, file = filename) 

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

#PHYLOSEQ OBJECT WITH ALL THE SAMPLES
#To make the work easier move all the phyloseq objects that you will merge to the same directory
#The directory must be exclusive for the sample phyloseq objects
psfpath <- "/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq/files" #write the path to the directory that contain the phyloseq objects
setwd(psfpath)
all_samp_ps <- list.files(path = ".", pattern = ".rds") #make a character vector with the names of the phyloseq objects
merged_all_ps <- list()
for(i in 1:length(all_samp_ps)){
  nps <- read_rds(all_samp_ps[i])
  merged_all_ps <<- merge_phyloseq(merged_all_ps, nps)
}

#Save phhyloseq object with all the taxas present
#write_rds(merged_all_ps, file="allphyla.rds") #The merged phyloseq object can be saved as it is

##PHYLOSEQ OBJECT WITH ONLY BACTERIA TAXA

setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq/diversity") #Set a new path to work if necessary

#Some of these steps were taken from the R script of Diego Garfias which can be found here https://github.com/nselem/mg-covid/blob/master/diego/taxonomy/analysis-covid-240521.R 
bacteria_ps <- subset_taxa(merged_all_ps, superkingdom == "Bacteria")  #Only Bacteria OTUs are left in the merged file

leaveout <- subset_taxa(bacteria_ps, family != "mitocondria" & class != "Chloroplast") #Subset from the bacteria object everything that is not a mitochondrial or chloroplast sequence
bacteria_ps <- prune_taxa(taxa_names(leaveout), bacteria_ps) #Leave only the reads of bacteria without chloroplast and mitochonria ones

#write_rds(bacteria_ps, file = "bacteria_cov.rds") #Save it to the rds file
#bacteria_ps <- read_rds("bacteria_cov.rds") #Read the rds file with the phyloseq object

#Retrieve the tables from the phyloseq object that only contains bacterial OTUs, bacteria_ps
otustab <- otu_table(bacteria_ps) #Table with the OTU abundance matrix
taxatab <- tax_table(bacteria_ps) #Table that contain the taxonomic classification of the OTUs

##ADD SAMPLE DATA TO THE PHYLOSEQ OBJECT
#setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq/diversity") #Set a new path if needed

#Read a csv file containing a table with the sample data or metadata
samptab <- read.csv("covid_isam_table.csv", header = TRUE, row.names = 1) #Metadata with interpretation, this file can be found in the same github

#Ensure that the sample names are the same in the sample data table as in the OTU table
#Change the name of the OTUtable columns, that correspond to the sample names, they are already sorted due to the list.files function
#Sort in the same manner, alphabetically, the SAMtable row or sample names, and use the to replace the names in the OTUtable
colnames(otustab) <- sort(rownames(samptab), decreasing = FALSE)

#Format the table with the sample data to phyloseq format
samptab <- sample_data(samptab)

#Update the phyloseq object with the metadata
bacteria_ps <- phyloseq(taxatab, otustab, samptab) #Add the metadata from the SAMtable and update the names from the OTUtable in the new bacteria pss

#Save changes
write_rds(bacteria_ps, file = "bacteria_cov.rds") #This file is ready to be used in the analysis done in phylobject.R

##NORMALIZATION
#These steps are a summary of the R script of Diego Garfias which can be found here https://github.com/nselem/mg-covid/blob/master/diego/taxonomy/analysis-covid-240521.R 

#Overview of the read depth in each sample
sample_sums(bacteria_ps) 

#Leaving out the samples below a certain treshold
bacteria_ps <- prune_samples(sample_sums(bacteria_ps) > 10000, bacteria_ps) #Many samples start at this order of magnitude, 10000, so any sample with this ammount of reads is kept

#Retrieve the tables from the phyloseq object with the prunned samples
otustab <- otu_table(bacteria_ps) #Table with the OTU abundance matrix
taxatab <- tax_table(bacteria_ps) #Table that contain the taxonomic classification of the OTUs
samptab <- sample_data(bacteria_ps) #Table that contain the sample data or metadata

#Function with the normalization method 
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

#Use the function on the phyloseq object that only contains bacterial OTUs and comparable samples
normbacps <- edgeRnorm(bacteria_ps, method = "TMM") #DGElist with your normalized data

#Make a new phyloseq object with normalized data
nOTUtable <- otu_table(normbacps@.Data[[1]], taxa_are_rows = TRUE) #Making a new otu table from the normalized data
nbacteria_ps <- merge_phyloseq(nOTUtable, taxatab, samptab) #Update the phyloseq object with the new otu table
nbacteria_ps <- write_rds(nbacteria_ps, file = "normbac.rds") #Save the phyloseq object in a rds file

#Save the tables from the normalized phyloseq object
write.csv(nOTUtable, "norm-otab.csv")
write.csv(taxatab, "norm-ttab.csv")
write.csv(samptab, "norm-stab.csv")
