#!/bin/bash 
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 23:59:00 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH --mem 50000 


module load bcftools 
prefix=sample
bcftools annotate --rename-chrs chromRename.txt -Oz -o $prefix.1.vcf.gz $prefix.vcf.gz
bcftools index -t $prefix.1.vcf.gz
bcftools view --no-header $prefix.1.vcf.gz > $prefix-noheader.vcf.gz
awk '{print $1"	"$2"	"".""	"$4"	"$5}' $prefix-noheader.vcf.gz > forCADD.vcf