---

source: md
title: "Analysis of taxonomic assignation data from a Covid project"

---

# Covid project. Correlation of mouth micorbiome diversity with covid infection

![image](https://user-images.githubusercontent.com/67386612/121607065-01e76d00-ca15-11eb-9a62-94f5b034646a.png)

*Words are a weapon stronger than he knows. And songs are even greater. The words wake the mind. The melody wakes the heart.
I came from a people of song and dance. I do not need him to tell me the power of words. But I smile nonetheless. Let's learn how to sing and dance*

	-Pierce Brown

I entitled myself to support a great and diverse team in searching for clues about how the microbiome affects the covid symptomatology.
One of the task assigned to me was to analyse the information after the taxonomic assignation of the reads. In order to do that, I used 
R and some of the packages that the community have developed:

| Software | Version | Manual | Description |
| -------- | ------------ | ------ | ------------- | 
| [phyloseq](https://github.com/joey711/phyloseq) | 1.36.0 | [Link](https://joey711.github.io/phyloseq/) | Explore, manipulate and analyze microbiome profiles with R |
| [ggplot2](https://cloud.r-project.org/web/packages/ggplot2/index.html) | 3.3.3 | [Link](https://ggplot2.tidyverse.org/) | System for declaratively creating graphics, based on The Grammar of Graphics |
| [vegan](https://cran.r-project.org/web/packages/vegan/index.html) | 2.5.7 | [Link](https://rdrr.io/cran/vegan/man/vegan-package.html) | Tools for descriptive community ecology | 
| [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) | 3.34.0 | [Link](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf) | Differential expression analysis of genomic data |

## Making preparations

I have traveled around different lessons of R in these six months, and something that they always mention is that organization of
your project is a top priority. I have experienced that this is true. That is why I encourage/urge anyone who wishes to follow this page,
to create a project for this, and future analysis with R.

Packages in R and other programming languages are like instruments in the shelve of a lab, you can take them(call them) and they will
provide a myriad of possibilities.
You can install the packages in the links provided. Let's begin by loading these packages into R:
~~~
library("phyloseq")
library("ggplot2")
library("edgeR")
library("RColorBrewer")
library("vegan")
~~~
{: .language-r}

If you decided to take the harsh path, let's define the directiry where our files will be allocated with the command `setwd()`.
In the case of my computer, I need to assing the folder with the next command:
~~~
setwd("C:/Users/cairo/Documents/sideprojects/covid/from-reads")
~~~
{: .language-r}

In this directory, you need to have two of the files provided in this github directory: `covid.biom` and `metadata-covid3`.
With this two files we can began the preparation of the data for the analysis.

## Loading and trimming the data 

We need to load the data into a phyloseq object with the command `import_biom`. We will save its content in a new objects called `covid`:
~~~
covid <- import_biom("non-human/covid.biom")
covid
~~~
{: .language-r}
~~~
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3851 taxa and 31 samples ]
sample_data() Sample Data:       [ 31 samples by 10 sample variables ]
tax_table()   Taxonomy Table:    [ 3851 taxa by 7 taxonomic ranks ]
~~~
{: .output}

Now we have our phyloseq object with the needed information. But if we take a look at the tax_table in it, we will see that 
needs some trimming:
~~~
head(covid@tax_table@.Data)
~~~
{: .language-r}
~~~
        Rank1         Rank2               Rank3                    Rank4                 Rank5                  
573     "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
571     "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
1463165 "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
1134687 "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
548     "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
2026240 "k__Bacteria" "p__Proteobacteria" "c__Gammaproteobacteria" "o__Enterobacterales" "f__Enterobacteriaceae"
        Rank6           Rank7               
573     "g__Klebsiella" "s__pneumoniae"     
571     "g__Klebsiella" "s__oxytoca"        
1463165 "g__Klebsiella" "s__quasipneumoniae"
1134687 "g__Klebsiella" "s__michiganensis"  
548     "g__Klebsiella" "s__aerogenes"      
2026240 "g__Klebsiella" "s__quasivariicola" 
~~~
{: .output}

We will get rid of the taxonomic rank identifier and rename the columns with the next commands:
~~~
covid@tax_table@.Data <- substring(covid@tax_table@.Data, 4)
colnames(covid@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(covid@tax_table@.Data)
~~~
{: .language-r}
~~~
~~~
{: .output}

~~~
        Kingdom    Phylum           Class                 Order              Family               Genus       
573     "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
571     "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
1463165 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
1134687 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
548     "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
2026240 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Klebsiella"
        Species          
573     "pneumoniae"     
571     "oxytoca"        
1463165 "quasipneumoniae"
1134687 "michiganensis"  
548     "aerogenes"      
2026240 "quasivariicola" 
~~~
{: .output}

That its much better for the analysis. 
Next, we need to load our metadata file:
~~~
meta <- read.csv("metadata-covid3.csv", row.names = 1)
str(meta)
~~~
{: .language-r}
~~~
'data.frame':	32 obs. of  10 variables:
 $ Age        : int  25 25 27 33 32 32 65 35 29 42 ...
 $ Sex        : chr  "F" "F" "F" "F" ...
 $ Pacient    : chr  "Positive" "Positive" "Positive" "Positive" ...
 $ Same       : int  1 1 0 0 2 2 0 0 0 3 ...
 $ Pacient.Age: chr  "20-30" "20-30" "20-30" "30-40" ...
 $ Symptoms   : chr  "SI" "SI" "SI" "SI" ...
 $ Saliva     : chr  "Positive" "Negative" "Negative" "Negative" ...
 $ Family     : chr  "Fam1" "Fam1" "Fam1" "Fam1" ...
 $ Control    : chr  "Positive" "Negative" "Negative" "Negative" ...
 $ Covid      : chr  "Positive" "P-Negative" "P-Negative" "P-Negative" ...
~~~
{: .output}

We can now take the row names from this data.frame, and use them to 
name our columns as well in the phyloseq object, this will be the names of the samples: 
~~~
colnames(covid@otu_table@.Data) <- row.names(meta)
sample_names(covid)
~~~
{: .language-r}
~~~
 [1] "SS07" "SS03" "SS02" "SS05" "SS06" "SS01" "SS04" "SS14" "SS15" "SS08" "SS16" "SS17" "SS09" "SS18" "SS12"
[16] "SS10" "SS19" "SS11" "SS20" "SS21" "SS13" "SS22" "SS23" "SS24" "SS25" "SS26" "SS27" "SS28" "SS29" "SS30"
[31] "SS31" "SS32"
~~~
{: .output}

With the metadata in an R object, we will tell R that it need to transform it in order to be used as part of a phyloseq object
with the `sample_data` command:
~~~
meta <- sample_data(meta)
~~~
{: .language-r}

And finally, we will merge this object whit the one where our tax and otu tables are located by the command `merge_phyloseq`:
~~~
covid <- merge_phyloseq(covid,meta)
covid
~~~
{: .language-r}
~~~
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3851 taxa and 32 samples ]
sample_data() Sample Data:       [ 32 samples by 10 sample variables ]
tax_table()   Taxonomy Table:    [ 3851 taxa by 7 taxonomic ranks ]
~~~
{: .output}

### Getting rid of the data that will bias/deviate our analysis

It is common to find reads that were assigned to another kingdom apart from bacteria. Also, 
it is not rare to have some mitochondrial and Chloroplast reads in metagenome
data. We will see if this is the case:
~~~
summary(covid@tax_table@.Data[,1] == "Bacteria")
summary(covid@tax_table@.Data[,5] != "mitocondria")
summary(covid@tax_table@.Data[,3] != "Chloroplast")
~~~
{: .language-r}
~~~
   Mode    TRUE 
logical    3851 

   Mode    TRUE 
logical    3851 

   Mode    TRUE 
logical    3851 
~~~
{: .output}

Luckily, we did not find any clue of misleading assigned reads. In the case that we have to face them, we could have used `subset_taxa` command to get rid of them:
~~~
covid <- subset_taxa(covid, Kingdom == "Bacteria")
covid <- subset_taxa(covid, Family != "mitocondria" & Class != "Chloroplast")
~~~
{: .language-r}

Another point to take into account is the depth of our sequencing, and how is distributed along with our samples. 
We will create a data.frame to allocate the data needed to see this information with `ggplot2`:
~~~
deepn <- data.frame(
  samples=as.character(map(strsplit(colnames(covid@otu_table@.Data), "_"),1)),
  reads= sample_sums(covid))
str(deepn)
~~~
{: .language-r}
~~~
'data.frame':	32 obs. of  2 variables:
 $ samples: chr  "SS07" "SS03" "SS02" "SS05" ...
 $ reads  : num  565701 1426441 1993276 3438792 1237777 ...
~~~
{: .output}

With the next command, we can create the plot with `ggplot2`. We will use a bar `geom` along with other specifications:
~~~
ggplot(data = deepn, aes(y= reads, x = samples))+
  geom_bar(stat="identity", fill="violetred4") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept = 4242198, col= "cyan3", size = 1.5, alpha = 0.5) +
  xlab("Samples") + ylab("Number of reads") +
  ggtitle("Deepness of covid libraries") +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))
~~~
{: .language-r}

![image](https://user-images.githubusercontent.com/67386612/121609078-c2228480-ca18-11eb-8573-b9a56f105345.png)

###### Figure 1. Depth of our samples with an arbitrary threshold line. 

~~~
~~~
{: .language-r}
~~~
~~~
{: .output}

