#!/bin/bash

#In order to run this program, the user needs to install the next programs(links to software pages inside parenthesis):
#	Spades(https://github.com/ablab/spades)
#	IDBA(https://github.com/loneknightpy/idba)
#	MaxBin(https://sourceforge.net/projects/maxbin/)
#	bowtie2(https://github.com/BenLangmead/bowtie2)
#	bwa(https://github.com/lh3/bwa)
#	bamtools(https://github.com/pezmaster31/bamtools)
#	samtools(https://github.com/samtools/samtools)
#	Kraken(https://github.com/DerrickWood/kraken2)
#	Bracken(https://github.com/jenniferlu717/Bracken)
#It is recommended to install this programs using conda, otherwise be sure to add them to the PATH

#------------------------------------------------------------------------------------------------------------------------------------------
#So as to run kraken-bracken, a database must be made, the next command can be used to create the standar kraken database.:
#	$ kraken-build --standard --db $DBNAME
#(Replace "$DBNAME" above with your preferred database name/location.)
#This will download NCBI taxonomic information, as well as the complete genomes in RefSeq for the bacterial, archaeal, and viral domains. 

#After downloading all this data, the build process begins; this is the most time-consuming step. 
#If you have multiple processing cores, you can run this process with multiple threads, e.g.:
#	$ kraken-build --standard --threads 24 --db $DBNAME

#After building the database, to remove any unnecessary files (including the library files no longer needed), run the following:
#	$kraken-build --db $DBNAME --clean
#------------------------------------------------------------------------------------------------------------------------------------------

#This program requies that the user specify 4 things in order. The first one is the name of the reads file with the forward sequences
#The sencond one is the name of the reads file with the reverse sequences. 
#The third one is a the prefix that will be added at the begenning of each output file.
#The forth one is the directory were the kraken database has been build.
#If the reads files are not in the same directory as this script, please provide it   

#As an example: sh metaTOmag.sh SSR3255-forward.fastq SSR3255-reverse.fastq coralloid-root

#It is important to mention that the binning-mapping-MAGs_assembly processes(mapMAXBIN.sh, mapping.sh) can run only if 
# main metagenome assembled has been done.

FILE1=$1 #Reads file with the forward sequences 
FILE2=$2 #Reads file with the reverse sequences 
prefix=$3 #Suffix to be appended at the begenning of the files 
krakendb=$4 #Path/name to krakendatabase

root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
sign='$'    

#Creates the outputs directories 
mkdir ASSEMBLIES
mkdir ASSEMBLIES/ASSEMBLY_LOGS
mkdir ASSEMBLIES/ASSEMBLY_LOGS/SCRIPTS

mkdir TAXONOMY 
mkdir TAXONOMY/KRAKEN
mkdir TAXONOMY/TAXONOMY_LOGS
mkdir TAXONOMY/TAXONOMY_LOGS/SCRIPTS

mkdir BINS
mkdir BINS/MaxBin
mkdir BINS/TAXONOMY 
mkdir BINS/TAXONOMY/KRAKEN
mkdir BINS/MAPPINGS
mkdir BINS/MAPPINGS/INDEXES
mkdir BINS/MAPPINGS/FILE
mkdir BINS/MAPPINGS/FILE/MAP
mkdir BINS/SPADES
mkdir BINS/BINS_LOGS
mkdir BINS/BINS_LOGS/SCRIPTS

mkdir MAPPINGS
mkdir MAPPINGS/OUTPUTS
mkdir MAPPINGS/SCRIPTS
mkdir MAPPINGS/INDEXES
mkdir MAPPINGS/FILE
mkdir MAPPINGS/FILE/MAP



#Metaspades assembler uses the best k-mer found and two more values, one smaller and one bigger. 
#For example, if the best k-mer value found is 51, the values that the assembler will use will be 49, 51, and 53.  
#Creates the file that will run the Metaspades:
cat > runMETASPADES.sh <<EOF
#!/bin/bash

fq2fa --merge --filter $FILE1 $FILE2 $root/MERGED_READS.fa

metaspades.py --pe1-1 $FILE1 --pe1-2 $FILE2 -o METASPADES

cp METASPADES/scaffolds.fasta ASSEMBLIES/${prefix}_metaspades_scaffolds.fasta 2>>/dev/null || cp METASPADES/contigs.fasta ASSEMBLIES/${prefix}_metaspades_contigs.fasta
EOF

#Creates the file that will run the the taxonomic assingment program kraken-bracken:
cat  > txnmKRAKEN_BRACKEN.sh <<EOF1
#!/bin/bash

kraken2 --db $krakendb --threads 12 --paired --fastq-input $FILE1 $FILE2 --output TAXONOMY/KRAKEN/${prefix}_kraken.kraken --report TAXONOMY/KRAKEN/${prefix}_kraken.report
bracken -d kraken-db -i TAXONOMY/KRAKEN/${prefix}_kraken.report -o TAXONOMY/KRAKEN/${prefix}.bracken

EOF1

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Creates the file that will run the the binning process and the mapping of these bins in order to assemble the MAGs. 
#Note that this file can only be run once the metagenome assembly is finished :
cat  > mapMAXBIN.sh <<EOF2
#!/bin/bash

ls ASSEMBLIES/*fasta | while read line; do assembler=${sign}(echo ${sign}line | cut -d"/" -f2 |cut -d"_" -f2); 
file=${sign}(echo ${sign}line | cut -d"/" -f2) ;
run_MaxBin.pl -thread 12 -contig ${sign}line -reads MERGED_READS.fa -out BINS/MaxBin/$prefix-${sign}assembler-MaxBin; done

ls BINS/MaxBin/*.fasta | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f3);
kraken2 --db kraken-db --threads 12 -input ${sign}line --output BINS/TAXONOMY/KRAKEN/${sign}file-kraken.kraken --report BINS/TAXONOMY/KRAKEN/${sign}file-kraken.report;
bracken -d kraken-db -i BINS/TAXONOMY/KRAKEN/${sign}file-kraken.report -o BINS/TAXONOMY/KRAKEN/${sign}file.bracken; 
bowtie2-build ${sign}line BINS/MAPPINGS/INDEXES/${sign}file; 
bowtie2 --threads 12 --sensitive-local -x BINS/MAPPINGS/INDEXES/${sign}file -1 $FILE1 -2 $FILE2 -S BINS/MAPPINGS/FILE/${sign}file.sam;
samtools view -F 4 -bS BINS/MAPPINGS/FILE/${sign}file.sam > BINS/MAPPINGS/FILE/${sign}file.bam; 
bamtools convert -in BINS/MAPPINGS/FILE/${sign}file.bam -format fastq > BINS/MAPPINGS/FILE/MAP/${sign}file.fastq;
sed -n '1~4s/^@/>/p;2~4p' BINS/MAPPINGS/FILE/MAP/${sign}file.fastq > BINS/MAPPINGS/FILE/MAP/${sign}file.fasta; 
mkdir BINS/SPADES/${sign}file;
spades.py -s BINS/MAPPINGS/FILE/MAP/${sign}file.fastq -o BINS/SPADES/${sign}file ;
cp BINS/SPADES/${sign}file/scaffolds.fasta BINS/SPADES/${sign}file.scaffolds.fasta 2>>/dev/null; done

ls BINS/MaxBin/*.tooshort | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f3);
kraken2 --db kraken-db --threads 12 -input ${sign}line --output BINS/TAXONOMY/KRAKEN/${sign}file-kraken.kraken --report BINS/TAXONOMY/KRAKEN/${sign}file-kraken.report;
bracken -d kraken-db -i BINS/TAXONOMY/KRAKEN/${sign}file-kraken.report -o BINS/TAXONOMY/KRAKEN/${sign}file.bracken; 
bowtie2-build ${sign}line BINS/MAPPINGS/INDEXES/${sign}file; 
bowtie2 --threads 12 --sensitive-local -x BINS/MAPPINGS/INDEXES/${sign}file -1 $FILE1 -2 $FILE2 -S BINS/MAPPINGS/FILE/${sign}file.sam;
samtools view -F 4 -bS BINS/MAPPINGS/FILE/${sign}file.sam > BINS/MAPPINGS/FILE/${sign}file.bam; 
bamtools convert -in BINS/MAPPINGS/FILE/${sign}file.bam -format fastq > BINS/MAPPINGS/FILE/MAP/${sign}file.fastq;
sed -n '1~4s/^@/>/p;2~4p' BINS/MAPPINGS/FILE/MAP/${sign}file.fastq > BINS/MAPPINGS/FILE/MAP/${sign}file.fasta; 
mkdir BINS/SPADES/${sign}file;
spades.py -s BINS/MAPPINGS/FILE/MAP/${sign}file.fastq -o BINS/SPADES/${sign}file ;
cp BINS/SPADES/${sign}file/scaffolds.fasta BINS/SPADES/${sign}file.scaffolds.fasta 2>>/dev/null; done


EOF2

#Creates the file that will run the the mapping process of the contigs-scaffolds from the assemble to the original reads. 
#Note that this file can only be run once the metagenome assembly is finished :
cat  > mapping.sh <<EOF3

ls ASSEMBLIES/*fasta | while read line; do assembler=${sign}(echo ${sign}line | cut -d"/" -f2 | cut -d'_' -f2); 
file=${sign}(echo ${sign}line | cut -d"/" -f2)
bowtie2-build ${sign}line MAPPINGS/INDEXES/${sign}assembler; 
bowtie2 --threads 12 --sensitive-local -x MAPPINGS/INDEXES/${sign}assembler -1 $FILE1 -2 $FILE2 -S MAPPINGS/FILE/${prefix}-${sign}assembler.sam;
samtools view -F 4 -bS MAPPINGS/FILE/${prefix}-${sign}assembler.sam > MAPPINGS/FILE/${prefix}-${sign}assembler.bam; 
bamtools convert -in MAPPINGS/FILE/${prefix}-${sign}assembler.bam -format fastq > MAPPINGS/FILE/MAP/${prefix}-${sign}assembler.fastq;
sed -n '1~4s/^@/>/p;2~4p' MAPPINGS/FILE/MAP/${prefix}-${sign}assembler.fastq > MAPPINGS/FILE/MAP/${prefix}-${sign}assembler.fasta; done

EOF3
 

#moves the TAXONOMY scripts to the corresponding sub-folder
mv txnm* TAXONOMY/TAXONOMY_LOGS/SCRIPTS

#Moves the assembly scripts to the corresponding folder
mv run* ASSEMBLIES/ASSEMBLY_LOGS/SCRIPTS

#Moves the binning scripts to the corresponding folder
mv binMAXBIN.sh BINS/BINS_LOGS/SCRIPTS

#Moves the mapping scripts to the corresponding folder
mv mapping.sh MAPPINGS/SCRIPTS