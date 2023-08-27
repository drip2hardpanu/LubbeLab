#!/bin/bash 
#SBATCH -A b1042
#SBATCH -p genomicslong 
#SBATCH -t 167:00:00 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH --mem 50000 

/projects/b1049/pranav5/moad/stratification/CADD-scripts/CADD.sh -g GRCh38 -o sample_out.tsv.gz forCADD.vcf

awk '{if ($6>12.37) print $1"_"$2"_"$3"_"$4}' sample_out.tsv.gz > 12.37_CADD.txt

awk '{if ($6>20) print $1"_"$2"_"$3"_"$4}' sample_out.tsv.gz > 20_CADD.txt
