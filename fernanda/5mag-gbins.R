##FIND GOOD BINS

#The following script helps you to do a fast identification 
#of good bins to use as MAGS
#It uses the small new file produced by the script 4mag-tocsv.sh
#It also list these good bins in a list, goodbins.ls

#Required libraries
library(stringr)
library(dplyr)

#Path to your file
path <- 'path/to/your/checkm/small/file' #Change this to your path
setwd(path)

#Read file
smfile <- read.csv("smallfile.csv", header = FALSE) #Read the file

#Format the data frame smfile

#Remove words to leave only values
torm <- c("bin", "'Completeness': ", "'Contamination': ") #words to remove must be equal in size to the number of columns
valdf <- t(apply(smfile, 1, str_remove_all, torm)) #Apply by row the str_remove_all function and using t to re-arrange columns and rows

#Name the columns
titulos <- c("bin", "completeness", "contamination") #column names
colnames(valdf) <- titulos

#Make numeric the completeness and contamination columns
valdf[ ,2:3] <- apply(valdf[ ,2:3], 2, as.numeric) #Apply by column the as.numeric function

valdf <- as.data.frame(valdf)

#Identify the quality bins
comp <- 80 #Mininum completeness, change it as it fits for you
cont <- 10 #Maximum contamination allowed, choose a treshold
good <- valdf %>% 
  filter(completeness > comp) %>% 
  filter(contamination < cont)
goodbins <- good$bin
write(goodbins, file = "goodbins.ls")
