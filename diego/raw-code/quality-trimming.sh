#!/bin/bash

#In order to run this program, the user needs to install the next programs(links to software pages inside parenthesis):
#	fastqC(https://github.com/s-andrews/FastQC)
#	Trimmomatic(https://github.com/timflutre/trimmomatic)
#It is recommended to install this programs using conda, otherwise be sure to add them to the PATH


#If the forward and reverse files are in more than one file each, it is neccesary to concatenate the files into two unique files.
#Usually the forward files have the "R1" characters on it, as the reverse has "R2". A good option for this is to use the next command:
#	$ cat *R1* > $filename-forward.fastq 
#-----------------------------------------------------------------------------------------------------------------------------------------

#This program requies that the user specify 2 things in order. The first one is the name of the reads file with the forward sequences
#The sencond one is the name of the reads file with the reverse sequences.

FILE1=$1 #Reads file with the forward sequences 
FILE2=$2 #Reads file with the reverse sequences 

root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
sign='$'

mkdir QUALITY
mkdir QUALITY/SCRIPTS

mkdir TRIM
mkdir TRIM/SCRIPTS

cat > quality-fastqc.sh <<EOF
#!/bin/bash

fastqc $FILE1 $FILE2 -o QUALITY/

EOF

#The parameters in the next program (LEADING, TRAILING, SLIDINGWINDOW, MINLEN) are especified with values that have been useful in 
# plant-microbiome analysis. It is important to read the Trimmomatic manual to understand these parameters
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

cat > trimming.sh << EOF1

trimmomatic PE $FILE1 $FILE2 p$FILE1 u$FILE1 p$FILE2 u$FILE2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

EOF1

#moves the quality-fastqc script to the corresponding sub-folder
mv quality* QUALITY/SCRIPTS

#Moves the trimming script to the corresponding folder
mv trimming* TRIM/SCRIPTS