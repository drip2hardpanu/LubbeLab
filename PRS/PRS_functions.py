##auxiliary functions for PRS analysis

##By Pranav Vijaykumar Patel
##Dr. Steven Lubbe Lab
#Northwestern FSM


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

##Data Processing
     
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
#chr1.OnlyPASS_DP10_GQ20_CR95.MultiSplit.SNPsInDels.RawID.bcftools.August2021.QCedPDCC_MAF01.vcf.gz
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
         "prefix="+vcfFile[:-7]
         "bcftools annotate --rename-chrs chromRename.txt -Oz -o $prefix.1.vcf.gz $prefix.vcf.gz\n"
         "bcftools index -t $prefix.1.vcf.gz\n"
    )
    newvcf = vcfFile[:-7]+".1.vcf.gz"
##stopping here for right now, because i still need to wait for header thing to be over to figure out if the pipeline will be faster doing the header step 

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

##gene burden functions

def postrun_geneburdenpipeline(inputFile, RVoutputFile):
    RVdf = pandas.read_table(RVoutputFile)
    countFiltered_df = RVdf[(RVdf['NumPolyVar'] > 1)]
    listOfGenes = list(countFiltered_df['Gene'])

    genes = ",".join(map(str, listOfGenes))

    with open(inputFile, 'r') as inputscript:
        withGene = [''.join([line.strip(), " --gene "+ genes, '\n']) for line in inputscript.readlines() if line[0:66] == "/projects/b1049/genetics_programs/Rvtest/rvtests/executable/rvtest"]
        print(withGene)
#postrun_geneburdenpipeline("300_SkatO.sh","group2.SkatO.assoc")
    sys.exit()
   # inputscript = open(inputFile)
    #for line in fileinput.input(files, inplace = 1): 
     # print line.replace("foo", "bar"),

      #  if line[0:66] == "/projects/b1049/genetics_programs/Rvtest/rvtests/executable/rvtest":
       #     newWithGene = line.replace(withGene) 

#run PRS.py -iPlink "naenae" -iPRS "amp-pdPRS.all_score" -p "naenae" -c "AMP-PD_2021_PDCC_Clean_Pheno_COVs_bbustos.wSITE.txt" -iVCF "naenae.vcf.gz" --covar-name "naenae,booboo" 300 600

def zPlotting(input_prs,covariateFile):
    """takes .best (PRSice) input and covariate file, ensure phenotype column is labeled PHENO, creates Z-distribution plot"""
    PRSdf = pandas.read_table(input_prs)
    covariatedf_preprocess = pandas.read_table(covariateFile)

    PRSdf.columns = PRSdf.columns.str.strip()

    covardf1 = covariatedf_preprocess.sort_values('FID')
    #sorted, ready to be pasted onto 
    sorted_covardf = covardf1.reset_index(drop = True)

    #for fixing data as pandas reads prsice data incorrectly
    PRSdf[["FID","IID","toDrop","PRS"]] = PRSdf.iloc[:,0].str.split(' ', expand=True)
    PRSdf = PRSdf.drop(columns = PRSdf.columns[0])

    PRSdf.columns.values[0] = "IID"
    PRSdf.columns.values[1] = "FID"
    PRSdf.columns.values[3] = "PRS"

    columns = ["FID","IID","PHENO","PRS"]

    #creates df -> FID IID PHENO PRS for conversion to different lines
    new_df = pandas.DataFrame({col: PRSdf[col] if col in PRSdf.columns else sorted_covardf[col] for col in columns})
    
    #cases

    df_cases = new_df[new_df.PHENO == 2]
    df_cases.to_csv('9_21_cases_covars_all.txt', sep='\t', index=False)



    #controls
    df_controls = new_df[new_df.PHENO == 1]
    df_controls.to_csv('9_21_controls_covars_all.txt', sep='\t', index=False)

    #plotting

    all_prs_scores = df_cases['PRS'].tolist() + df_controls['PRS'].tolist()
    new_all = [float(x) for x in all_prs_scores]
    num_ticks = 10  # Set the number of desired ticks
    tick_positions = [min(new_all) + (i * (max(new_all) - min(new_all)) / num_ticks) for i in range(1, num_ticks + 1)]
    print(tick_positions)
    #tick_positions = np.linspace(min(new_all), max(new_all), num_ticks + 1)
    #print(tick_positions)

    tick_positions = [-0.00013,-.006]
    plt.figure()
    plt.hist([df_cases['PRS'], df_controls['PRS']], bins=20, alpha=0.5, label=['Cases', 'Controls'], color=['blue', 'green'])
    plt.xticks([-.05, 0,.05],rotation='vertical')
    plt.xlabel('PRS Score')
    plt.ylabel('Frequency')
    plt.title('Distribution of PRS Scores')
    plt.legend()



    plt.show()

