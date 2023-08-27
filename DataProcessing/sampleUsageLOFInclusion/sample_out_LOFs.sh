#!/bin/bash 
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH --mem 50000 

module load bcftools
bcftools view -i 'ID=@LOFs.txt' --threads 10 -Oz -o LOFs.vcf.gz sample.vcf.gz
bcftools index -t LOFs.vcf.gz 

