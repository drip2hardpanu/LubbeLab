import os
os.system("module load python/anaconda3.6 \n")

import sys
import pandas
import fileinput

#sbatch parameters
def sbatchWritingBlank(self):
    """function to write common sbatch parameters, manipulations will show up in all script-generating functions"""
    self.write(
        "#!/bin/bash \n"
        "#SBATCH -A  \n"
        "#SBATCH -p  \n"
        "#SBATCH -t  \n"
        "#SBATCH -N  \n"
        "#SBATCH -n  \n"
        "#SBATCH --mem  \n"
        "\n"
    )

def sbatchWritingLong(self):
    """function to write common sbatch parameters, manipulations will show up in all script-generating functions"""
    self.write(
        "#!/bin/bash \n"
        "#SBATCH -A b1042\n"
        "#SBATCH -p genomicslong \n"
        "#SBATCH -t 167:00:00 \n"
        "#SBATCH -N 1 \n"
        "#SBATCH -n 10 \n"
        "#SBATCH --mem 50000 \n"
        "\n"
    )

def sbatchWritingShort(self):
    """function to write common sbatch parameters, manipulations will show up in all script-generating functions"""
    self.write(
        "#!/bin/bash \n"
        "#SBATCH -A b1042\n"
        "#SBATCH -p genomics\n"
        "#SBATCH -t 4:00:00 \n"
        "#SBATCH -N 1 \n"
        "#SBATCH -n 10 \n"
        "#SBATCH --mem 50000 \n"
        "\n"
    )

def sbatchWritingMid(self):
        self.write(
        "#!/bin/bash \n"
        "#SBATCH -A b1042\n"
        "#SBATCH -p genomics\n"
        "#SBATCH -t 10:00:00 \n"
        "#SBATCH -N 1 \n"
        "#SBATCH -n 10 \n"
        "#SBATCH --mem 50000 \n"
        "\n"
    )
        
def sbatchWritingMidLong(self):
        self.write(
        "#!/bin/bash \n"
        "#SBATCH -A b1042\n"
        "#SBATCH -p genomics\n"
        "#SBATCH -t 23:59:00 \n"
        "#SBATCH -N 1 \n"
        "#SBATCH -n 10 \n"
        "#SBATCH --mem 50000 \n"
        "\n"
    )
####LiftOver preperation functions ###
def readyForLift(input_file):
    import sys
    import os

    file = open(input_file, 'r')
    newFile = open(input_file[:-4]+"_1.txt", 'w')

    for line in file.readlines():
        i = str(line.rstrip())
        n = len(i)
        for char in range(n):
            if ":" in i:
                if i[char] == ":":
                    newFile.write(
                        str(i[:char])+"\t"+str(i[char+1:])+"\t"+str(i[char+1:])+"\n")
                    
    file.close()
    newFile.close()


def readyForLift_offline(input_file):
    import sys
    import os

    file = open(input_file, 'r')
    newFile = open(input_file[:-4]+"_1.bed", 'w')

    for line in file.readlines():
        i = str(line.rstrip())
        n = len(i)
        for char in range(n):
            if ":" in i:
                if i[char] == ":":
                    newFile.write(
                        str(i[:char])+"\t" +str(i[char+1:]) + "\t"+str(i[char+1:]) +"\n")
                    
    file.close()
    newFile.close()

##VCF to PLINK converters
def commandLineVCFtoPlink(vcfFile, CovariateFile):
    """requires IID, FID ... covariate file"""
    prefix = vcfFile[:-7]
    os.system("module load plink")
    #turning vcf into bed
    os.system("plink2 --vcf "+vcfFile+"--make-bed --out "+prefix+"_unUpdated")

    covariatedf = pandas.read_table(CovariateFile)
    covariatedf_head = list(covariatedf.columns)
    
    #creating phenotype file
    index_pheno = covariatedf_head.index("PHENO")
    os.system('awk '"'{print "'$1"\t"$2"\t"$'+str(index_pheno+1)+"}' "+CovariateFile+" > "+prefix+"_pheno.txt")

    #creating sex file
    index_sex = covariatedf_head.index("SEX")
    os.system('awk '"'{print "'$1"\t"$2"\t"$'+str(index_sex+1)+"}' "+CovariateFile+"> "+prefix+"_sex.txt")

    os.system("plink2 --bfile "+prefix+"_unUpdated --update-sex "+prefix+"_sex.txt --pheno "+prefix+"_pheno.txt --out "+prefix)
    os.system("rm "+prefix+"_unUpdated*")

def scriptVCFtoPlink(vcfFile, CovariateFile):
    prefix = vcfFile[:-7]
    script = open(prefix+"vcf_plink.sh", 'w')
    sbatchWritingMidLong(script)

    script.write(
         "module load plink"
         "\n"
    )

    #covariate dataframe for columns
    covariatedf = pandas.read_table(CovariateFile)
    covariatedf_head = list(covariatedf.columns)

    #creating phenotype file
    index_pheno = covariatedf_head.index("PHENO")
    os.system('awk '"'{print "'$1"\t"$2"\t"$'+str(index_pheno+1)+"}' "+CovariateFile+" > "+prefix+"_pheno.txt")

    #creating sex file
    index_sex = covariatedf_head.index("SEX")
    os.system('awk '"'{print "'$1"\t"$2"\t"$'+str(index_sex+1)+"}' "+CovariateFile+"> "+prefix+"_sex.txt")

    script.write(
        "plink2 --vcf "+vcfFile+" --make-bed --out "+prefix+"_unUpdated"
        "\n"
        "plink2 --bfile "+prefix+"_unUpdated --update-sex "+prefix+"_sex.txt --pheno "+prefix+"_pheno.txt --out "+prefix+
        "\n"
        "rm "+prefix+"_unUpdated*"
            )
    
###Data processing function###
def cadAnno_VCFfiltering(vcfFile, build, prefix, thresholdsList):
    annotation = open(prefix+"_cadAnno.sh", 'w')
    wrapper = open(prefix+"_CADDwrapper.sh", 'w')


    ##annotation file
    sbatchWritingLong(annotation)
    annotation.write("/projects/b1049/pranav5/moad/stratification/CADD-scripts/CADD.sh -g "+build+" -o "+prefix+".tsv.gz "+vcfFile)

    ##thresholds
    n = len(thresholdsList)
    for i in range(n):
        annotation.write("\nawk '{if ($6>"+str(thresholdsList[i])+") print $1"'"_"'"$2"'"_"'"$4"'"_"'"$3}' "+ prefix+".tsv.gz > "+str(thresholdsList[i])+"_CADD.txt")
    annotation.close()
    wrapper.write(
        "jid0=($(sbatch "+prefix+"_cadAnno.sh))\n"
        "echo "'"jid0 ${jid0[-1]}"'" >> slurm_ids\n"
         )

    for i in range(n):
        script= open(str(thresholdsList[i])+"_filtering.sh", 'w')
        sbatchWritingShort(script)
        
        script.write(
            "module load bcftools\n"
            "bcftools view -i 'ID=@"+str(thresholdsList[i])+"_CADD.txt' --threads 10 -Oz -o "+ prefix+"_"+str(thresholdsList[i]) +"_CADD.vcf.gz "+ vcfFile +"\n"
            "bcftools index -t "+ prefix+"_"+str(thresholdsList[i]) +"_CADD.vcf.gz\n"
            "\n"  
        )

        ##bcftools doesn't manipulate the file directly, meaning it is okay for two jobs to be running on the same file as each job is reading and outputting to another file
        script.close()
        
        wrapper.write(
            "jid"+str(i+1)+"=($(sbatch --dependency=afterok:${jid0[-1]} "+str(thresholdsList[i])+"_filtering.sh))\n"

        )    

def allLOFinclusion_pipeline(vcfFile, build, prefix, thresholdsList):
    """dependency shell script pipeline generator for damaging variant filtration + LOF inclusion"""
    cadannotation = open(prefix+"_cadAnno.sh", 'w')
    annovarannotation = open(prefix+"_annovarAnno.sh", 'w')
    LOFextraction = open(prefix+"_LOFs.sh", 'w')
    wrapper = open(prefix+"_CADDwrapper.sh", 'w')
    wrapper.write("!/bin/bash -e \n")

    #creates a file + script for Chromosome Rename
    chr = open('rename_chromosomes.txt', 'w')
    for x in range(1,23):
         chr.write("chr"+str(x)+" "+str(x)+"\n")
        
    chr.close()

    rename_script = open("chromosomeRename.sh", 'w')
    sbatchWritingMid(rename_script)

    rename_script.write(
         "\n"
         "module load bcftools \n"
         "prefix="+vcfFile[:-7]+
         "bcftools annotate --rename-chrs chromRename.txt -Oz -o $prefix.1.vcf.gz $prefix.vcf.gz\n"
         "bcftools index -t $prefix.1.vcf.gz\n"
    )
    newvcf = vcfFile[:-7]+".1.vcf.gz"

    ##creates an annotation file for cad (to field damaging variant scores)
    sbatchWritingLong(cadannotation)
    cadannotation.write(
        "/projects/b1049/pranav5/moad/stratification/CADD-scripts/CADD.sh -g "+build+" -o "+prefix+".tsv.gz "+vcfFile +"\n"
        )

    #creates an annotation file for annovar (for functional annotation)
    sbatchWritingMid(annovarannotation)
    annovarannotation.write(
        "module load bcftools\n"
        "bcftools view --no-header "+vcfFile+" > "+ prefix+"_noHeader.vcf.gz \n"
        "\n"
        "awk '{print $1"'"\t"'"$2"'"\t"'"$2"'"\t"'"$4"'"\t"'"$5}' "+ prefix+"_noHeader.vcf.gz > forAnno.txt\n"
        "\n"
        "perl /projects/b1049/pranav5/moad/stratification/annovar/table_annovar.pl forAnno.txt /projects/b1049/pranav5/moad/stratification/annovar/humandb/ -buildver hg38 -protocol refGene -remove -otherinfo -operation g -out " + prefix+"_annotated\n"
        "\n"
        "awk '{if ($6 == "'"splicing"'" || $8 =="'"frameshift substitution"'" || $8 =="'"stopgain"'" || $8 =="'"stoploss"'") print $1"'"_"'"$2"'"_"'"$4"'"_"'"$5}' "+prefix+"_annotated.hg38_multianno.txt > LOFs.txt \n"
        "\n"    
    )

    ##thresholds; filter by cadd thresholds
    n = len(thresholdsList)
    for i in range(n):
        cadannotation.write("\nawk '{if ($6>"+str(thresholdsList[i])+") print $1"'"_"'"$2"'"_"'"$4"'"_"'"$3}' "+ prefix+".tsv.gz > "+str(thresholdsList[i])+"_CADD.txt\n"
        )

    #creates a wrapper
    cadannotation.close()
    wrapper.write(
        "jid0=($(sbatch "+prefix+"_cadAnno.sh))\n"
        "echo "'"jid0 ${jid0[-1]}"'" >> slurm_ids\n"
        "\n"
        "jid1=($(sbatch "+prefix+"_annovarAnno.sh))\n"
        "echo "'"jid1 ${jid1[-1]}"'" >> slurm_ids\n"
         )

    sbatchWritingShort(LOFextraction)
    
    #pulls LOFs from the vcf file
    LOFextraction.write(
            "module load bcftools\n"
            "bcftools view -i 'ID=@LOFs.txt' --threads 10 -Oz -o LOFs.vcf.gz "+ vcfFile +"\n"
            "bcftools index -t LOFs.vcf.gz \n"
            "\n"  
        )
    wrapper.write(
            "jid"+str(2)+"=($(sbatch --dependency=afterok:${jid0[-1]} "+prefix+"_LOFs.sh))\n"
            "echo "'"jid2 ${jid2[-1]}"'" >> slurm_ids\n"
        )

    for i in range(n):
        script= open(str(thresholdsList[i])+"_filtering.sh", 'w')
        sbatchWritingShort(script)
        
        script.write(
            "module load bcftools\n"
            "bcftools view -i 'ID=@"+str(thresholdsList[i])+"_CADD.txt' --threads 10 -Oz -o "+ prefix+"_"+str(thresholdsList[i]) +"_CADD.vcf.gz "+ vcfFile +"\n"
            "bcftools index -t "+ prefix+"_"+str(thresholdsList[i]) +"_CADD.vcf.gz\n"
            "\n"
            "bcftools merge --force-samples "+prefix+"_"+str(thresholdsList[i]) +"_CADD.vcf.gz LOFs.vcf.gz > "+prefix+"_LOFs_wDuplicates.vcf.gz \n"
            "bcftools norm -D "+prefix+"_LOFs_wDuplicates.vcf.gz --threads 5 -Oz -o " + prefix+"_LOFs.vcf.gz\n"
            "bcftools index -t "+prefix+"_LOFs.vcf.gz"   
        )    
        ##bcftools doesn't manipulate the file directly, meaning it is okay for two jobs to be running on the same file as each job is reading and outputting to another file
        script.close()
        
        wrapper.write(
            "jid"+str(i+3)+"=($(sbatch --dependency=afterok:${jid2[-1]} "+str(thresholdsList[i])+"_filtering.sh))\n"
        )