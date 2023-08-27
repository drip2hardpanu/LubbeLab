#!/bin/bash 
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH --mem 50000 

module load bcftools

awk '{print $1"	"$2"	"$2"	"$4"	"$5}' forCADD.vcf > forAnno.txt

perl /projects/b1049/pranav5/moad/stratification/annovar/table_annovar.pl forAnno.txt /projects/b1049/pranav5/moad/stratification/annovar/humandb/ -buildver hg38 -protocol refGene -remove -otherinfo -operation g -out sample_out_annotated

awk '{if ($6 == "splicing" || $8 =="frameshift substitution" || $8 =="stopgain" || $8 =="stoploss") print $1"_"$2"_"$4"_"$5}' sample_out_annotated.hg38_multianno.txt > LOFs.txt 

