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

##CREATE THE INITIAL OBJECTS

#Read the phyloseq file with the normalized bacterial OTUs
#The process to make this file is in the transform_merge.R script found in this same directory
setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq/diversity") #Set path to the merged file
covid_ps <- read_rds("normbac.rds") #Read the rds file, it is a phyloseq object
samtable <- sample_data(covid_ps) #Extract the metadata as a table
samcols <- colnames(samtable) #This vector contains all the variables from the metadata, even the ones of interpretation
ipvars <- c("Family", "Sex", "Pacient.Age", "Condition", "Stage") #Only contain the variables of interest showed in the article

##STATISTICAL TESTS ON ALPHA INDEXES 

#Richness data
richness <- plot_richness(covid_ps, measures = c("Chao1", "Shannon", "Simpson"))
richness <- richness$data

#METADATA VARIABLES

#In the following statistical test of the alpha indexes only "ipvars" will be used.
#in order to avoid repetition and to facilitate the information visualization.

#ALPHA INDEXES MEANS

#The following function retrieves the means of each index (Shannon, Simpson, Chao1)
#grouped by variable "gvar", the metadata variable to obtain the mean, wrote between ""
#richness is a table that contains the metadata variables with its respective alpha diversity values, value
#and an identifier for the alpha diversity index used, variable
indmean <- function(gvar, richness) {
  inval <- richness$value
  index <- richness$variable
  categ <- richness[ ,gvar]
  im <- tapply(X=inval, INDEX = list(index, categ), mean)
  return(im)
}

#Variable means

#TABLE S1. Alpha indexes means
#A list to contain the variable means
imlist <- vector("list", length(ipvars)) 
names(imlist) <- ipvars #to preserve the identity of each variable mean in the list
#For cycle to obtain the means of each variable of interest
for(i in 1:length(ipvars)){ 
  imi <- indmean(gvar = ipvars[i], richness = richness)
  imlist[[i]] <- imi
}


#Subset richness data by alpha diversity index
H <- richness %>% filter(variable == "Shannon") #Just Shannon values
S <- richness %>% filter(variable == "Simpson") #Just Simpson values
C <- richness %>% filter(variable == "Chao1") #Just Shannon values

#KRUSKAL-WALLIS

#Kruskal Wallis test for each index subset grouped by a chosen variable, in this case, those in "ipvars"

#A function that will display a list with the results of each variable to the Kruskal-Wallis test
#for a given set of index values, indata
#This will help explore if there is any statistical significant difference in a given variable
#the interest variable(s) are added in a character vector "ipvars"
kt_indata <- function(indata, ipvars) {
  klist <- vector("list", length(ipvars)) #list to save the results of the for cycle
  names(klist) <- ipvars #to preserve the identity of each variable tested in the list
  for(i in 1:length(ipvars)){
    groupi <- indata[,ipvars[i]]
    kti <- kruskal.test(indata$value, groupi)
    klist[[i]] <- as.vector(kti)
  } 
  return(klist)
}

#TABLE S2. Alpha indexes Kruskal-Wallis results
#Kruskal-Wallis results for each variable in the sample data of the phyloseq object, covid_ps
#The results are saved as list of lists
H_kw_res <- kt_indata(indata = H, ipvars = ipvars) #For the Shannon indexes
S_kw_res <- kt_indata(indata = S, ipvars = ipvars) #For the Simpson indexes
C_kw_res <- kt_indata(indata = C, ipvars = ipvars) #For the Chao1 indexes

#Save a list in a text file
#"lista" is the list to be saved
#"filename" is the name of the produced text file, wrote it between ""
listext <- function(lista, filename){
  sink(filename)
  print(lista)
  sink()
}

#Save the list with the results on the alpha test in text files
listext(imlist, "norm-amean.txt") #Save the means of all the alpha values grouped by each variable of interest category
listext(H_kw_res, "norm-hkw.txt") #Save the kruskal wallis results of the shannon values grouped by each variable of interest category
listext(S_kw_res, "norm-skw.txt") #Save the kruskal wallis results of the simpson values grouped by each variable of interest category
listext(C_kw_res, "norm-ckw.txt") #Save the kruskal wallis results of the chao1 values grouped by each variable of interest category

##BETA DIVERSITY

#To obtain the beta diversity two things are needed: a distance matrix and an ordination method.
#The graphic representations of the beta diversity were done in other analysis found here https://github.com/nselem/mg-covid/blob/master/diego/taxonomy/analysis-covid-240521.R

#DISTANCES 
#These can vary due to the method used to calculate them.
#Distances can be calculated directly in the ordination function, but they are required as individual object for the PERMANOVA
#Two methods were chosen for the distances:
jc_dist <- phyloseq::distance(covid_ps, method = "jaccard") #Jaccard is a presence/absence based distance
bc_dist <- phyloseq::distance(covid_ps, method = "bray") #Bray-Curtis is an abundance based distance
#Other useful distance method is Unifrac unweighted/weighted, a phylogeny based one, it requires a tree in the phyloseq object

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
#PERMANOVA results for all the metadata variables even the ones not included in the article
jc_pmav <- pmav_all(jc_dist) #for Jaccard distances
bc_pmav <- pmav_all(bc_dist) #for Bray-Curtis distances

#TABLE S3. PERMANOVA results
jcipmav <- jc_pmav[ipvars] #PERMANOVA results using Jaccard distances for the variables in the article
bcipmav <- bc_pmav[ipvars] #PERMANOVA results using Bray-Curtis distances for the variables in the article

listext(jcipmav, "pmav-jcnorm.txt") #Save the Jaccard PERMANOVA results in a text file
listext(bcipmav, "pmav-bcnorm.txt") #Save the Bray-Curtis PERMANOVA results in a text file

##PAIRED SAMPLES

#Subset the samples that belong to the same individuals sampled in different times
#Leave only those that in different times were in different stages, post-covid and infected
same <- prune_samples(samtable$Same != 0 & samtable$Same != 2, covid_ps)
pmet <- sample_data(same)

#Alpha diversity
#This will plot the richness values grouped by stage
prplot <- plot_richness(same, x = "Stage", 
                        measures = c("Chao1", "Shannon", "Simpson"),
                        color = "Stage") +
  geom_boxplot()
#Save the richness values in a table
parich <- prplot$data
#Table S4. Alpha diversity means for the paired data set
#Explore the alpha diversity means
indmean(gvar = "Stage", richness = parich)
#Separate the richness values by the type of used index
Cpar <- parich %>% filter(variable == "Chao1")
Hpar <- parich %>% filter(variable == "Shannon")
Spar <- parich %>% filter(variable == "Simpson")
#TABLE S4. ALfa diversity wilcoxon test results for the paired data set
#Compare them trough a Wilcoxon test for paired samples
wilcox.test(Cpar$value ~ Cpar$Stage, paired = TRUE)
wilcox.test(Hpar$value ~ Hpar$Stage, paired = TRUE)
wilcox.test(Spar$value ~ Spar$Stage, paired = TRUE)

#Beta diversity
#Generate the distance matrixes
par_jdist <- phyloseq::distance(same, method = "jaccard")
par_bdist <- phyloseq::distance(same, method = "bray")
#Generate the ordination method
par_bord <- phyloseq::ordinate(same, method = "PCoA", distance = "bray")
par_jord <- phyloseq::ordinate(same, method = "PCoA", distance = "jaccard")
#Plots
#FIG 4. PCoA done with Bray-Curtis distances
par_bcoa <- plot_ordination(same, par_bord, color = "Stage") + 
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Stage),
               type = "norm", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  geom_point(size=3, alpha=0.8) +
  geom_text(label = pmet$Same, size=5, vjust = "outward", hjust = "outward") +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()
#TABLE S4. PERMANOVA results for the beta diversity of the paired samples
set.seed(100)
bppmav <- adonis(par_bdist ~ pmet$Stage, permutations = 1000)
set.seed(100)
jppmav <- adonis(par_jdist ~ pmet$Stage, permutations = 1000)

#DIFFERENTIAL ABUNDANCE

#ANALYSIS OF COMPOSITION MICROBIOMES OR ANCOM
#code modified from https://github.com/FrederickHuangLin/ANCOM-BC-Code-Archive/blob/master/scripts/figure_6_table_s1_s2_s3.Rmd
# and https://rdrr.io/github/FrederickHuangLin/ANCOMBC/f/vignettes/ANCOMBC.Rmd

#Group OTUS to a taxa rank
genera_ps <- aggregate_taxa(covid_ps, "genus")

#Using ANCOM BC
#The following function sets the parameters used in this work
#It also automates the replacement of the phyloseq object used, taxrank, which is the phyloseq object grouped to a taxa rank of preference
#and of the group variable, ivar, which is could be any variable of interest present in sample data of the phyloseq object 
ancovid <- function(taxrank, ivar){
  ancombc(phyloseq = taxrank, formula = ivar, 
          p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
          group = ivar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
          max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
}
#A warning will be displayed for the variables with less than 3 categories, this is not a problem, it only means that the global analysis will not be done, but in those cases it is not needed

#Adjusted log fold change and standard error tables
#The following function produces a table with the comparisons between one of the categories with the rest of them in the chosen group variable
#Be warned, if your group variable has more than two categories and you want to compare between all of them all the steps are need to be repeated
lf_sd_tabs <- function(ancomres){
  #Part 1 retrieve and adjust the log fold change between the pairwise comparisons
  dfig_pt1 <- data.frame(ancomres$res$beta * ancomres$res$diff_abn, check.names = FALSE) %>%
    rownames_to_column("taxon_id")
  #Part 2 retrieve and adjust the standard error between the pairwise comparisons
  dfig_pt2 <- data.frame(ancomres$res$se * ancomres$res$diff_abn, check.names = FALSE) %>%
    rownames_to_column("taxon_id")
  colnames(dfig_pt2)[-1] <- paste0(colnames(dfig_pt2)[-1], "SD") #Names to distinguish each value
  #Paste together part 1 and 2, filter the comparison you wish to plot
  dfig_c <- dfig_pt1 %>%
    left_join(dfig_pt2, by = "taxon_id") 
  dfig_c$taxon_id <- factor(dfig_c$taxon_id, levels = dfig_c$taxon_id)
  return(dfig_c)
}


#Function to filter and order the log fold change and standard error table, lfsdtab
#to make it work you need a object of vector class with the column names, compar
#This is required to do the plot
filtandord <- function(lfsdtab, compar) {
  lfsdtab %>%
    filter(.data[[compar[[2]]]] != 0) %>% #here we keep only the taxa with a log fold change different to zero, for that we use the second element of compar, that refers to coefficient column
    arrange(desc(abs(.data[[compar[[2]]]]))) #arrange the rows in descending order by their absolute value
}

#Log fold change plot
lfc_plot <- function(plotab, compar){
  ggplot(data = head(plotab, n=5L), #Data is adjusted to plot only the top five taxa changes
           aes(x = Taxon, y = .data[[compar[[2]]]])) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = .data[[compar[[2]]]] - SD, 
                    ymax = .data[[compar[[2]]]] + SD), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  labs(x=NULL, y="Log Fold Change") + 
  coord_flip()
}


#ANCOM BC results

#The variables of interest have three categories
#i.e. a, b and c
#So there will be three pairwise comparisons
#i.e. a vs b, a vs c, b vs c
#ANCOMBC only compares the first category in ascending order against the rest
#i.e. a vs b, a vs c
#To make the missing comparison, a new variable must be made with the categories renamed, where a new category is the first
#i.e old <- c(a, b, c), new <- c(control, b, c)
#So now the missing comparison can be done
#i.e. b vs control, b vs c
#If you do not want to do this process you can simply make sure that your control category is the first one in ascending
#Or whichever category your are interested to compare, must be the first in ascending alfabetically order

#Condition
#This variable has three categories: "Asymptomatic" "Symptomatic"  "Control"
#Three possible comparisons
cond_comparAC <- c("Taxon", "AsymptomaticVSControl", "SD") #Asymptomatic vs Control
cond_comparAS <- c("Taxon", "AsymptomaticVSSymptomatic", "SD") #Asymptomatic vs Symptomatic
cond_comparCS <- c("Taxon", "ControlVSSymptomatic", "SD") #Control vs Symptomatic
#New "Condition" variable, "Cond"
#Categories renamed as follows: Asymptomatic = NO, Symptomatic = SI, and Control = Control, making control the new first element
Cond <- c(replace(x = samtable$Symptoms, list = which(samtable$Condition == "Control"), values = "Control"))
#Add Cond variable to the sample_data of the phyloseq object
sample_data(geneera_ps)$Cond <- Cond
#Group by genus
cond_ge_ancA <- ancovid(genera_ps, "Condition") #list with ANCOM results 
cond_ge_tabA <- lf_sd_tabs(cond_ge_ancA) #unfiltered table with fold change with two out of three comparisons 
sample_data(genera_ps)$Cond <- Cond #Add the Cond column to the phyloseq object grouped by genera
cond_ge_ancB <- ancovid(genera_ps, "Cond") #analysis with the renamed categories of Condition
cond_ge_tabB <- lf_sd_tabs(cond_ge_ancB) #unfiltered table with fold change, two comparisons, one is the  missing one
#Make the table and plots for each comparison
#Control vs Asymptomatic
cond_ge_tabAC <- cond_ge_tabB[ ,c(1,2,4)]
colnames(cond_ge_tabAC) <- cond_comparAC
cond_ge_tabAC <- filtandord(cond_ge_tabAC, cond_comparAC)
cond_ge_ploAC <- lfc_plot(cond_ge_tabAC, cond_comparAC)
#Asymptomatic vs Symptomatic
cond_ge_tabAS <- cond_ge_tabA[ ,c(1,3,5)]
colnames(cond_ge_tabAS) <- cond_comparAS
cond_ge_tabAS <- filtandord(cond_ge_tabAS, cond_comparAS)
cond_ge_ploAS <- lfc_plot(cond_ge_tabAS, cond_comparAS)
#Control vs Symptomatic
cond_ge_tabCS <- cond_ge_tabB[ ,c(1,3,5)]
colnames(cond_ge_tabCS) <- cond_comparCS
cond_ge_tabCS <- filtandord(cond_ge_tabCS, cond_comparCS)
cond_ge_ploCS <- lfc_plot(cond_ge_tabCS, cond_comparCS)
