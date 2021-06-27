#PBS -N kraken2
#PBS -q default
#PBS -l nodes=1:ppn=12,mem=48g,vmem=48g,walltime=99:99:99
#PBS -V
#PBS -e 7kraken2.err
#PBS -o 7kraken2.out

cd /PATH/TO/DIRECTORIES/THAT/CONTAIN/KRAKEN/OUTPUT/AND/METABAT/BINS #Change this before executing, it is the directory that contain all your sample files

module load kraken/2.0.7

mkdir krakenmag/

#You need the goodbins.ls file produced by the 5mag-gbins.R script
cat goodbins.ls | while read line; do kraken2 --db /PATH/TO/YOUR/KRAKEN/DATABASE --threads 12 --output krakenmag/${line}.kraken --report krakenmag/${line}.report metabat/${line}.fa; done
