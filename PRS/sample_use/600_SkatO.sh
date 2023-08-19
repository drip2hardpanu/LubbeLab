#!/bin/bash 
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH --mem 50000 

/projects/b1049/genetics_programs/Rvtest/rvtests/executable/rvtest --inVcf input.vcf.gz --pheno 600_PRSids.txt --pheno-name PHENO --covar 600_PRSids.txt --out 600_rareBurden.sh --covar-name covariate1,covariate2,covariate3 --kernel skato