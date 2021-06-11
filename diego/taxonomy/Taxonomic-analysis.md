---

source: md
title: "Analysis of taxonomic assignation data from a Covid project"

---

# Covid project. Correlation of mouth micorbiome diversity with SARS-CoV-2 infection

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

By this result, we see that the sample SS09 is underrepresented so we will leave it out of the analysis. 
Let's use the next lines of code to take its information out of the metadata file
and out of our phyloseq object:
~~~
meta <- meta[-13,]

covid.otu <- covid@otu_table@.Data
covid.otu <- as.data.frame(covid.otu)
covid.otu <- select(covid.otu, -13)

colnames(covid.otu) <- row.names(meta)

covid.otu <- otu_table(covid.otu, taxa_are_rows = TRUE)

covid.tax <- covid@tax_table@.Data
covid.tax <- tax_table(covid.tax)

covid <- merge_phyloseq(covid.otu, covid.tax,meta)
~~~
{: .language-r}

You can see that we needed to subset the individual parts of our phyloseq
file (otu_table, tax_table, sam_data), reassign them to a phyloseq format, and finally merged them
into a new covid phyloseq object.

### Normalization of the data

With our trimmed and selected samples. We need to normalize them in order to compare them
and do the analysis. In order to normalize our data, we used code and
knowledge provided by [McMurdie et al., 2014](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531). We highly recommend you to
give their paper a review. But in the mean time, we used this part of their code, a function:
~~~
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
str(z)
~~~
{: .language-r}
~~~
Formal class 'DGEList' [package "edgeR"] with 1 slot
  ..@ .Data:List of 2
  .. ..$ : num [1:3851, 1:31] 147378 14 10 8 9 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:3851] "573" "571" "1463165" "1134687" ...
  .. .. .. ..$ : chr [1:31] "SS07" "SS03" "SS02" "SS05" ...
  .. ..$ :'data.frame':	31 obs. of  3 variables:
  .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..$ lib.size    : num [1:31] 569552 1430292 1997127 3442643 1241628 ...
  .. .. ..$ norm.factors: num [1:31] 1.687 1.06 0.839 0.949 1.168 ...
~~~
{: .output}

We obtained a DGEList object that we need to transform into an object that phyloseq can read,
and merge. Again, we will use the `otu_table` and `merge_phyloseq` command:
~~~
covid.otu2 <- otu_table(z@.Data[[1]], taxa_are_rows = TRUE)
covid2 <- merge_phyloseq(covid.otu2, covid.tax,meta)
~~~
{: .language-r}

We will transform the absolute counts to relative abundance so as to compare 
the samples between them. We will use a useful command in the phyloseq package`transform_sample_counts` to accomplish this:
~~~
head(covid2@otu_table@.Data)[,c(1:6)]
~~~
{: .language-r}
~~~
          SS07   SS03   SS02   SS05   SS06  SS01
573     147378 117636 106001 111841 109690 96293
571         14     20     31     51     18    14
1463165     10      9     14      6     10     9
1134687      8      5      6     12      5     7
548          9     10     20     49     15     4
2026240      5      7      5     11      1     3
~~~
{: .output}
~~~
covid2 <- transform_sample_counts(physeq = covid2, function(x) x*100/sum(x))
head(covid2@otu_table@.Data)[,c(1:6)]
~~~
{: .language-r}
~~~
                SS07         SS03         SS02         SS05         SS06         SS01
573     25.876127202 8.2246142746 5.3076744744 3.2486958421 8.834369e+00 2.173912e+01
571      0.002458072 0.0013983159 0.0015522298 0.0014814199 1.449710e-03 3.160642e-03
1463165  0.001755766 0.0006292421 0.0007010070 0.0001742847 8.053942e-04 2.031841e-03
1134687  0.001404613 0.0003495790 0.0003004316 0.0003485694 4.026971e-04 1.580321e-03
548      0.001580189 0.0006991579 0.0010014386 0.0014233250 1.208091e-03 9.030405e-04
2026240  0.000877883 0.0004894106 0.0002503596 0.0003195219 8.053942e-05 6.772804e-04
~~~
{: .output}

## Beta diversity analysis

Phyloseq have a useful catalog of functions, among them is `ordinate`, which lead you
to construct a distance matrix for beta diversity analysis, with a desired `method` 
and distance `metric`. We will use [NMDS](https://academic.oup.com/bioinformatics/article/21/6/730/199398) as our method, and [Bray-Curtis](http://www.pelagicos.net/MARS6300/readings/Bray_&_Curtis_1957.pdf) as our distance metric.
There is a sensible reason for this, but a discussion on that issue is beyond the scope of this little document. Let's just mention,
that you can display all the distance metrics phyloseq can use by typing the command `distanceMethodList`
~~~
covid.ord <- ordinate(physeq = covid2,method = "NMDS", distance = "bray")
head(covid.ord$points)
~~~
{: .language-r}
~~~
           MDS1        MDS2
SS07 -0.2373059  0.06957890
SS03 -0.2015262 -0.01158747
SS02 -0.1498378  0.02389948
SS05 -0.1049601 -0.11142705
SS06 -0.1201024 -0.15001822
SS01 -0.1104110  0.15863121
~~~
{: .output}

Now we have a two dimentional representation of the diversity among our data. We will use `ggplot2` to 
plot it:
~~~
plot_ordination(covid2, covid.ord, color = "Covid") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Covid),
               type = "norm", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  theme_bw()+
  geom_point(size=4, alpha=0.8) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis")
~~~
{: .language-r}

![image](https://user-images.githubusercontent.com/67386612/121611145-36f7bd80-ca1d-11eb-82b8-ebd2289166fe.png)

##### Figure 2. Beta diversity by NMDS with Bray-Curtis distance of the data

## Multivariate analysis of variance with vegan
In order to use vegan for the multivariate analysis, we will extract our data from the phyloseq objects, and we will use the 
command `t()` to obtain the transposed data-frame as `vegan` needs it:
~~~
d.covid <- t(covid2@otu_table@.Data)
~~~
{: .language-r}

Also, we will extract the metadata from our `meta` object, and put it in a new data-frame:
~~~
meta.covid <- data.frame(Age = as.factor(covid2@sam_data@.Data[[5]]),
                         Sex = as.factor(covid2@sam_data@.Data[[2]]),
                         Pacient = as.factor(covid2@sam_data@.Data[[3]]),
                         Symt = as.factor(covid2@sam_data@.Data[[6]]),
                         Diag = as.factor(covid2@sam_data@.Data[[7]]),
                         Covid = as.factor(covid2@sam_data@.Data[[10]]))
str(meta.covid)
~~~
{: .language-r}
~~~
'data.frame':	31 obs. of  6 variables:
 $ Age    : Factor w/ 5 levels "<18",">55","20-30",..: 3 3 3 4 4 4 2 4 3 5 ...
 $ Sex    : Factor w/ 2 levels "F","M": 1 1 1 1 2 2 2 1 2 2 ...
 $ Pacient: Factor w/ 2 levels "Negative","Positive": 2 2 2 2 2 2 2 2 2 2 ...
 $ Symt   : Factor w/ 2 levels "NO","SI": 2 2 2 2 1 1 1 2 1 1 ...
 $ Diag   : Factor w/ 2 levels "Negative","Positive": 2 1 1 1 1 1 1 1 1 2 ...
 $ Covid  : Factor w/ 3 levels "Negative","P-Negative",..: 3 2 2 2 2 2 2 2 2 3 ...
~~~
{: .output}

Finally, with the function `adonis`, we can se the relationship that has a specific variable over the variance of
the data(abundances):
~~~
adonis(d.covid ~ Diag , data = meta.covid, permutations = 999)
~~~
{: .language-r}
~~~
Call:
adonis(formula = d.covid ~ Diag, data = meta.covid, permutations = 999) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Diag       1    0.2328 0.23283  1.2921 0.04266  0.236
Residuals 29    5.2255 0.18019         0.95734       
Total     30    5.4583                 1.00000
~~~
{: .output}

Repeating the last command with different equations as models, we can construct a table where the information concerning the results can be allotted

###### Table 1. Results from the analysis with adonis

![image](https://user-images.githubusercontent.com/67386612/121611894-f39e4e80-ca1e-11eb-9d63-d6e261816c4e.png)

This took me early-mornings and nights to [grok](https://en.wikipedia.org/wiki/Grok). This is the main reason why I desired to share it 
not only by the code, but with an explanatory document to help other fellow bioinformatician to cope with the academic life.

