###Adding some interpretations to the metadata of the sample table

#READING FILES

#Set path
setwd("/home/mfcg/Descargas/covid/microbiome/taxonomia/phyloseq")

#The original table was found here https://github.com/nselem/mg-covid/blob/master/diego/taxonomy/metadata-covid3.csv

#Reading data
SAMtable <- read.csv("metadata-covid3.csv", header = TRUE, row.names = 1)

#CREATE CONDITION COLUMN
#Condition indicates if a Pacient is asymtomatic, symtomatic or control

#Required columns from the original data
Pacient <- SAMtable$Pacient 
Symptoms <- SAMtable$Symptoms

#Create index by condition
cindex <- which(Pacient == "Negative") #control index
sindex <- which(Pacient == "Positive" & Symptoms == "SI") #symtomatic index
aindex <- which(Pacient == "Positive" & Symptoms == "NO") #asymtomatic index

#Create an empty vector for condition the same length as the vectors that we are using to create it
Condition <- vector(mode = "character", length = 32) 
#Fill the empty vector with the condtion name using the index
Condition[cindex] <- rep("Control", length(cindex)) 
Condition[aindex] <- rep("Asymptomatic", length(aindex))
Condition[sindex] <- rep("Symptomatic", length(sindex))

#Add the condition vector as a column to a new metadata table
newSAMtable <- cbind(SAMtable, Condition)

#CREATE STAGE
#Indicates the patients disease stage at the moment of the microbiome sampling
#Stage can be control, infection and post-infection

#Required columns from the original data 
Pacient <- SAMtable$Pacient
PCR <- SAMtable$Saliva

#Create index by moment
cindex <- which(Pacient == "Negative") #control index
iindex <- which(Pacient == "Positive" & PCR == "Positive") #infection index
pindex <- which(Pacient == "Positive" & PCR == "Negative") #post-infection index

#Create an empty vector for stage
Stage <- vector(mode = "character", length = length(Pacient))
#Fill with the appropiate stage 
Stage[cindex] <- rep("Control", length(cindex))
Stage[iindex] <- rep("Infection", length(iindex))
Stage[pindex] <- rep("Post-infection", length(pindex))

#Add the stage vector as a column to a new metadata table
newSAMtable <- cbind(newSAMtable, Stage)

#Save new table
write_csv(newSAMtable, file = "covid_isam_table.csv")
