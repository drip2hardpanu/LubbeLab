## script generator + executor for stratifying by MAF
##auxiliary functions on stratification_functions.py

##Created by Pranav Patel on June 10th, 2023
##Dr. Steven Lubbe Lab
##Northwestern FSM 

import os
import sys

n = len(sys.argv)
listOfOptions = ["-p", "--prefix", "-i", "--input"]
listofMAFs = []
x=0

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print("\n -h help")
    print(
        "\n stratifying.py takes command line input:"
        '\n example usage: run stratification.py -p "LOFs" -i "input" .1 .2 .4 "singletons" '
        "\n -i, --input : input gzipped vcf PREFIX, use quotation marks"
        "\n -p, --prefix : prefix for output files, use quotation marks"
    )
    sys.exit()

##input flag

if "-i" in sys.argv[1:] or "--input" in sys.argv[1:]:
    for i in range(1,n):
        if sys.argv[i] == "-i" or sys.argv[i] == "--input":
            inputVCF = sys.argv[i+1]
            print("inputVCF == " + inputVCF + ".vcf.gz")
else:
    print("stratifying.py requires input flag [-i] or [--input]")
    sys.exit()

##prefix flag

if "-p" in sys.argv[1:] or "--prefix" in sys.argv[1:]:
    for i in range(1,n):
        if sys.argv[i] == "-p" or sys.argv[i] == "--prefix":
            prefix = sys.argv[i+1]
            print("prefix == " + prefix)
else:
    print("stratifying.py requires prefix flag [-p] or [--prefix]")
    sys.exit()

##maf frequencies
z = 1
while z < n:
    if sys.argv[z] in listOfOptions:
        z += 2
    elif sys.argv[z] == "singletons":
        z += 1
    else:
        try:
            float((sys.argv[z]))
        except:
            print(sys.argv[z] + " is not a valid MAF, try again")
            sys.exit()
        listofMAFs.append(float(sys.argv[z]))
        z += 1

lenOfLOM=len(listofMAFs)
print("\n MAF frequencies calculated:", end = " ")
if "singletons" in sys.argv:
    print("singletons", end=" ")
for maf in range(lenOfLOM):
    print(listofMAFs[maf], end = " ")
 
#manipulate parameters as needed
p = open("stratification.sh", 'w')
p.write("#!/bin/bash \n"
        "#SBATCH -A b1042 \n"
        "#SBATCH -p genomics \n"
        "#SBATCH -t 2:00:00 \n"
        "#SBATCH -N 1 \n"
        "#SBATCH -n=10 \n"
        "#SBATCH --mem=16000 \n"
        "\n"
        )

#####################
p.write("module load bcftools \n"
        "module load plink \n"
        "\n"
        )
#################
#handling everything but singletons
listofMAFs.sort(reverse=True)
LenSort= len(listofMAFs)

if LenSort > 0:
    q = 0
    while q != LenSort:
        #taking from the input file
        if q == 0:
            p.write("bcftools view -Q "+ str(listofMAFs[q]) + ":minor --threads 10 -Oz -o "+str(listofMAFs[q])+"_"+prefix+".vcf.gz " +inputVCF+".vcf.gz \n"
                    "bcftools index -t "+str(listofMAFs[q])+"_"+prefix+".vcf.gz \n"
                    "\n"                
                    )
            q += 1
        else:
            p.write("bcftools view -Q "+ str(listofMAFs[q]) + ":minor --threads 10 -Oz -o "+str(listofMAFs[q])+"_"+prefix+".vcf.gz " +str(listofMAFs[q-1])+"_"+prefix+".vcf.gz \n"
                    "bcftools index -t "+str(listofMAFs[q])+"_"+prefix+".vcf.gz \n"
                    "\n"                
                    )
            q += 1
    
    if "singletons" in sys.argv:
        p.write(
            "plink2 --vcf "+str(listofMAFs[LenSort-1])+"_"+prefix+".vcf.gz --make-bed --out temp \n"
            "plink2 --bfile temp --freq counts \n"
            "awk '{if ($5 == 1) print $2}' plink2.acounts > inclusion.txt \n"
            "bcftools view -i 'ID=@inclusion.txt' --threads 10 -Oz -o singletons"+prefix+".vcf.gz " +str(listofMAFs[LenSort-1])+"_"+prefix+".vcf.gz"
            )
        
#singletons only
elif LenSort == 0 and "singletons" in sys.argv:
    p.write(
    "plink2 --vcf "+inputVCF+".vcf.gz --make-bed --out temp \n"
    "plink2 --bfile temp --freq counts \n"
    "awk '{if ($5 == 1) print $2}' plink2.acounts > inclusion.txt \n"
    "bcftools view -i 'ID=@inclusion.txt' --threads 10 -Oz -o singletons"+prefix+".vcf.gz " +inputVCF+".vcf.gz"
    )

else:
    print("no mafs were given!")
    sys.exit()

p.close()
os.system("sbatch stratification.sh")
