###METAGENOMIC ANALYSIS 

#Some parts of this analysis were modified from the code/lesson retrieved from https://carpentries-incubator.github.io/metagenomics/09-diversity-analysis/index.html

##LIBRARIES
library("ggplot2")
library("readr")
library("phyloseq")
library("patchwork")
library("vegan")
library("dplyr")
library("ANCOMBC")
library("tidyverse")
library("microbiome")

##INITIAL EXPLORATION OF THE MERGED FILE

#Read merged file
#The merged file can be created with the transform_merge.R script found here https://github.com/nselem/mg-covid/blob/master/fernanda/transform_merge.R
setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq") #Set path to the merged file
covid_ps <- read_rds("covid_ups.rds") #Read the merged rds file


#PHYLA ABUNDANCE PLOTS

#Plot the phyla present in your samples
#Absolute counts
abs_count_plot <- plot_bar(covid_ps, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat = "identity", position = "stack") +
  ggtitle("Absolute abundance")

#Percent
percentages <- transform_sample_counts(covid_ps, function(x) x*100 / sum (x))
perc_plot <- plot_bar(percentages, fill = "phylum") +
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack")



##ALPHA DIVERSITY PLOTS

#The following function plot a alpha diversity plot grouping by a variable, vgroup
#vgroup must be a character vector of length 1
#the indexes the function will use are chao1, shannon and simpson
alpha_plot <- function(vgroup){
  plot_richness(covid_ps, x = vgroup, 
                measures = c("Chao1", "Shannon", "Simpson"))
}

#PLOTS GROUP BY A VARIABLE

#Some variables are interpretations, for more info about them, check this script https://github.com/nselem/mg-covid/blob/master/fernanda/interpretation_sam_data.R
#and the data https://github.com/nselem/mg-covid/blob/master/fernanda/covid_isam_table.csv

PacAlpha <- alpha_plot("Pacient") #By the patient initial test to COVID
SymAlpha <- alpha_plot("Symptoms") #By the pressence of symptoms 
#PCRalpha <- alpha_plot("Saliva") #By the result of the PCR result to COVID at the moment of the microbiome. Replaced by "Stage"
StaAlpha <- alpha_plot("Stage") #By the patient stage of the COVID disease (not infected, infected, post infected). Replaced PCR
FamAlpha <- alpha_plot("Family") #By each family
ConAlpha <- alpha_plot("Condition") #By the condition of the COVID disease (control, symtomatyc, asymptomatyc)
SexAlpha <- alpha_plot("Sex") #By the patient sex

##STATISTICAL TESTS ON ALPHA INDEXES 

#Richness data
richness <- plot_richness(covid_ps, measures = c("Chao1", "Shannon", "Simpson"))
richness <- richness$data

#Subset richness data by alpha diversity index
H <- richness %>% filter(variable == "Shannon") #Just Shannon values
S <- richness %>% filter(variable == "Simpson") #Just Simpson values
C <- richness %>% filter(variable == "Chao1") #Just Shannon values

#KRUSKAL-WALLIS

#Kruskal Wallis test for each index subset by a chosen variable
#Create a vector that contain the variables to group the elements for the Kruskal-Wallis test
samtable <- sample_data(covid_ps)
samcols <- colnames(samtable) #This vector contains all the variables from the metadata, even the ones of interpretation

#A function that will display a list with the results of each variable to the Kruskal-Wallis test
#for a given set of index values, indata
#This will help explore if there is any significative difference in a variable
kt_indata <- function(indata) {
  klist <- vector("list", length(samcols)) #list to save the results of the for cycle
  names(klist) <- samcols #to preserve the identity of each variable tested in the list
  for(i in 1:length(samcols)){
    groupi <- indata[,i]
    kti <- kruskal.test(indata$value, groupi)
    klist[[i]] <- as.vector(kti)
  } 
  return(klist)
}

#Kruskal-Wallis results for each variable in the sample data of the phyloseq object, covid_ps
H_kw_res <- kt_indata(indata = H) #For the Shannon indexes
S_kw_res <- kt_indata(indata = S) #For the Simpson indexes
C_kw_res <- kt_indata(indata = C) #For the Chao1 indexes

#WILCOXON

#Wilcoxon test between some groups of interest
#Only "Saliva" variable, also called "Stage", was chosen. Being the only one with a difference <0.1 in the Kruskal test
Ssplit <- split(S, S$Saliva) #Split the groups by the PCR result 
NvC_simsal <- wilcox.test(Ssplit$Negative$value, Ssplit$Negative_Control$value) #Negative(or post-covid) vs control
PvC_simsal <- wilcox.test(Ssplit$Positive$value, Ssplit$Negative_Control$value) #Positive(or infected) vs control
PvN_simsal <- wilcox.test(Ssplit$Positive$value, Ssplit$Negative$value) #Positive(infected) vs negative(post-covid)



##BETA DIVERSITY

#To obtain the beta diversity two things are needed: distances and ordination method.

#DISTANCES 
#These can vary due to the method used to calculate them.
#Distances can be calculated directly in the ordination function, but they are required as individual object for the PERMANOVA
#Two methods were chosen for the distances:
jc_dist <- phyloseq::distance(covid_ps, method = "jaccard") #Jaccard is a presence/absence based distance
bc_dist <- phyloseq::distance(covid_ps, method = "bray") #Bray-Curtis is an abundance based distance
#Other useful distance method is Unifrac unweighted/weighted, a phylogeny based one, it requires a tree in the phyloseq object

#ORDINATION METHOD
#The chosen method for this work is Principal Coordinates Analysis or PCoA
jc_coa_ord <- phyloseq::ordinate(covid_ps, method = "PCoA", distance = "jaccard") #Jaccard distances
bc_coa_ord <- phyloseq::ordinate(covid_ps, method = "PCoA", distance = "bray") #Bray distances
#Double Principal Coordinate Analysis or DPCoA can also be used but requires a tree in the phyloseq file

#PCoA PLOTS BY 

#Saliva PCR results
jc_saliva <- plot_ordination(covid_ps, jc_coa_ord, color = "Saliva") + stat_ellipse() #Jaccard distances
bc_saliva <- plot_ordination(covid_ps, bc_coa_ord, color = "Saliva") + stat_ellipse() #Bray-Curtis distances

#Family
jc_fam <- plot_ordination(covid_ps, jc_coa_ord, color = "Family") + stat_ellipse() #Jaccard distances
bc_fam <- plot_ordination(covid_ps, bc_coa_ord, color = "Family") + stat_ellipse() #Bray-Curtis distances

#Condition
jc_con <- plot_ordination(covid_ps, jc_coa_ord, color = "Condition") + stat_ellipse() #Jaccard distances
bc_con <- plot_ordination(covid_ps, bc_coa_ord, color = "Condition") + stat_ellipse() #Bray-Curtis distances

#Pacient
jc_pac <- plot_ordination(covid_ps, jc_coa_ord, color = "Pacient") + stat_ellipse() #Jaccard distances
bc_pac <- plot_ordination(covid_ps, bc_coa_ord, color = "Pacient") + stat_ellipse() #Bray-Curtis distances

#Sex
jc_sex <- plot_ordination(covid_ps, jc_coa_ord, color = "Sex") + stat_ellipse() #Jaccard distances
bc_sex <- plot_ordination(covid_ps, bc_coa_ord, color = "Sex") + stat_ellipse() #Bray-Curtis distances



#PERMUTATIONAL MULTIVARIATE ANALYSIS OF VARIANCE

#The following script list the PERMANOVA results for each variable of a given sample_data table
# for a phyloseq distances object created by the ordination method of choice
pmav_all <- function(distances){ #distances is distance object created by phyloseq::distance function
  pmavlist <- vector("list", length(samcols)) #empty list to save the PERMANOVA results
  names(pmavlist) <- samcols #sample_data variables
  for(i in 1:length(samcols)){
    rhsi <- c(as.vector(unlist(samtable[,i]))) #create independent variable vector from the columns of the sample_data table
    set.seed(100) #establish a seed so the p-values are the same in each repetition
    pmavi <- adonis(distances ~ rhsi, permutations = 1000)
    pmavlist[[i]] <- pmavi
  }
  return(pmavlist)
}

#PERMANOVA RESULTS
#List with the PERMANOVA results using every variable of the sample data of the phyloseq object, covid_ps
#These variables are used as independent variables, RHS, in the PERMANOVA, to test a dissimilarity matrix, LHS
#See adonis function, ?adonis
jc_pmav <- pmav_all(jc_dist) #PERMANOVA results using Jaccard distances
bc_pmav <- pmav_all(bc_dist) #PERMANOVA results using Bray-Curtis distances



#DIFFERENTIAL ABUNDANCE

#ANALYSIS OF COMPOSITION OR ANCOM
#code modified from https://github.com/FrederickHuangLin/ANCOM-BC-Code-Archive/blob/master/scripts/figure_6_table_s1_s2_s3.Rmd
# and https://rdrr.io/github/FrederickHuangLin/ANCOMBC/f/vignettes/ANCOMBC.Rmd

#Group to a taxa rank
rankovid <- aggregate_taxa(covid_ps, "phylum")

#Analisys of composition, choosing your variable and values
ancovid <- ancombc(phyloseq = rankovid, formula = "Stage", 
                   p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                   group = "Stage", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                   max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
#Make tables with the data to plot
#Part 1 retrieve and adjust the log fold change between the pairwise comparisons
dfig_pt1 <- data.frame(ancovid$res$beta * ancovid$res$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
#Part 2 retrieve and adjust the standard error between the pairwise comparisons
dfig_pt2 <- data.frame(ancovid$res$se * ancovid$res$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(dfig_pt2)[-1] <- paste0(colnames(dfig_pt2)[-1], "SD") #Names to distinguish each value
#Paste together part 1 and 2, filter the comparison you wish to plot
dfig_c <- dfig_pt1 %>%
  left_join(dfig_pt2, by = "taxon_id") %>%
  filter(StageInfection != 0) %>%
  arrange(desc(StageInfection))
dfig_c$taxon_id <- factor(dfig_c$taxon_id, levels = dfig_c$taxon_id)

#Log fold change plot
ggplot(data = dfig_c, #Data is adjusted to plot only the top twenty changes
           aes(x = taxon_id, y = StageInfection)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = StageInfection - StageInfectionSD, 
                    ymax = StageInfection + StageInfectionSD), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  labs(x=NULL, y="Log Fold Change") + 
  coord_flip()


samp_frac <- ancovid$samp_frac
samp_frac[is.na(samp_frac)] <- 0

log_obs_abn <- log(abundances(covid_ps) + 1)
log_oabn_adj <- t(t(log_obs_abn)-samp_frac)

