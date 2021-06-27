#!/bin/bash

READ1=$1 #Contains the forward reads in fastq.gz format
READ2=$2 #Contains the reverse reads in fastq.gz format
prefix=$3 #An assigned name to identify your sample
path=$(pwd)
sign='$'

cat > 1bowsam-${prefix}.sh << EOF
#PBS -N bowtie${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=24g,vmem=24g,walltime=99:00:00
#PBS -e 1bowsam${prefix}.err
#PBS -o 1bowsam${prefix}.out
#PBS -V

cd $path

module load bowtie2/2.3.5.1 samtools/1.9

mkdir bowtie/

bowtie2-build contig.fa bowtie/${prefix}.ind

bowtie2 -p 8 -x bowtie/${prefix}.ind -1 $READ1 -2 $READ2 -S bowtie/${prefix}.map.sam

samtools sort -o bowtie/${prefix}.map.sorted.bam -O bam bowtie/${prefix}.map.sam

EOF

cat > 2metabat-${prefix}.sh << EOF
#PBS -N MetaBat${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=99:00:00
#PBS -e 2metabat${prefix}.err
#PBS -o 2metabat${prefix}.out
#PBS -V

cd $path

module load metabat/2.13

mkdir metabat/

jgi_summarize_bam_contig_depths --outputDepth metabat/${prefix}depth.txt bowtie/${prefix}.map.sorted.bam

metabat2 --saveCls -i contig.fa -a metabat/${prefix}depth.txt -o metabat/${prefix}bin -m 2500 --maxEdges 600 --noAdd --seed 4

EOF

cat > 3checkm-${prefix}.sh << EOF
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=9:59:59,vmem=32gb,mem=32gb
#PBS -N check${prefix}
#PBS -V
#PBS -o 3check${prefix}.out
#PBS -e 3check${prefix}.err

cd $path

module load CheckM/1.1.3 Prodigal/2.6.2 hmmer/3.1b2

mkdir checkm/

checkm taxonomy_wf domain Bacteria -t 8 -x fa metabat/ checkm/

EOF
