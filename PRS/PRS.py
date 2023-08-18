##script generator + executor for PRS analysis 
##auxiliary functions on PRS_functions.py

##Created by Pranav Patel on July 4th, 2023
##Dr. Steven Lubbe Lab 
##Northwestern FSM 

##maybe will need to add covariate and phenotype file

#read file
#extract the threshold and the individuals
#put them in their own files
#possible run gbba

#only works with python/anaconda3.6 on FSM quest

import os
import sys
import pandas

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

n = len(sys.argv)
listOfOptions = ["-p", "--prefix", "-inputPRS", "-iPRS", "-iPlink", "--inputPlink", "--covar", "-c", "-iVCF", "--inputVCF", "--covar-name"]
individualThresholds = []
x=0

bottomIDs = []
topIDs = []

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
  print("\n -h help")
  print(
      "\n PRS.py takes command line input:"
      '\n example usage: run stratification.py -p "LOFs" -inputPlink "input" -iPRS "input.all_score" 300 600'
      "\n -iPlink, --inputPlink : input gzipped vcf.gz, use quotation marks"
      "\n -iPRS, --inputPRS : input file of PRS scores (PRSice --all-scores output)"
      "\n -p, --prefix : prefix for output files, use quotation marks"
      "\n -c, --covar : covariate file with phenotype denoted as PHENO, use quotation marks"
  )
  sys.exit()

##input plink flag

if "-iPlink" in sys.argv[1:] or "--inputPlink" in sys.argv[1:]:
  for i in range(1,n):
      if sys.argv[i] == "-iPlink" or sys.argv[i] == "--inputPlink":
          plinkInput = sys.argv[i+1]
          print("input plink == " + plinkInput + "\n")
else:
  print("stratifying.py requires input plink flag [-iPlink] or [--inputPlink]")
  sys.exit()


if "-iPRS" in sys.argv[1:] or "--inputPRS" in sys.argv[1:]:
  for i in range(1,n):
      if sys.argv[i] == "-iPRS" or sys.argv[i] == "--inputPRS":
          PRSinput = sys.argv[i+1]
          print("inputPRS == " + PRSinput+ "\n")
else:
  print("stratifying.py requires input PRS flag [-iPRS] or [--inputPRS]")
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

if "-c" in sys.argv[1:] or "--covar" in sys.argv[1:]:
  for i in range(1,n):
      if sys.argv[i] == "-c" or sys.argv[i] == "--covar":
          covariate = sys.argv[i+1]
          print("covar == " + covariate)
      
          COVdf = pandas.read_table(covariate)
          COVdf.columns = COVdf.columns.str.strip()
else:
  print("stratifying.py requires covariate flag [-c] or [--covar]")
  sys.exit()

if "--covar-name" in sys.argv[1:]:
  for i in range(1,n):
      if sys.argv[i] == "--covar-name":
          covariateNames  = sys.argv[i+1]
          print("covariate names == " + covariateNames + "\n")
else:
  print("no covariate names were provided, all information columns (headers) will be treated as covariates")

  notCovarNames = ["FID","IID","MATID","FATID","PHENO"]
  covarnames_prestring = []

  covarnames_temp = list(COVdf.columns)

  for covar in covarnames_temp:
      if covar not in notCovarNames:
          covarnames_prestring.append(covar)

  covariateNames = ",".join(map(str, covarnames_prestring))
  print("covariate names == " + covariateNames + "\n")

if "-iVCF" in sys.argv[1:] or "--inputVCF" in sys.argv[1:]:
  for i in range(1,n):
      if sys.argv[i] == "-iVCF" or sys.argv[i] == "--inputVCF":
          iVCF  = sys.argv[i+1]
          print("input VCF == " + iVCF + "\n")
else:
  print("stratifying.py requires input vcf flag [-iVCF] or [--inputVCF]")
  sys.exit()


##individual thresholds
z = 1
while z < n:
  if sys.argv[z] in listOfOptions:
      z += 2
  else:
      try:
          int((sys.argv[z]))
      except:
          print(sys.argv[z] + " is not a valid threshold, try again")
          sys.exit()
      individualThresholds.append(int(sys.argv[z]))
      z += 1

lenOfThresholds=len(individualThresholds)
print("\n Thresholds calculated:", end = " ")
for ind in range(lenOfThresholds):
  print(individualThresholds[ind], end = " ")

##dataframing for sorting
os.system("awk '{print $1"'"\t"'"$2"'"\t"'"$3}' "+PRSinput+" > " + PRSinput+"_updated.txt")
updatedPRS = PRSinput+"_updated.txt"
PRSdf = pandas.read_table(updatedPRS)

#cleaning cols
PRSdf.columns.values[0] = "IID"
PRSdf.columns.values[1] = "FID"
PRSdf.columns.values[2] = "PRS"
PRSdf.columns = PRSdf.columns.str.strip()

#bottom to top
PRSdf1 = PRSdf.sort_values(by=['PRS'])
sortedPRSdf = PRSdf1.reset_index(drop = True)

#top to bottom
PRSdf2 = PRSdf.sort_values(by=['PRS'], ascending=False)
sortedPRSdf_reverse = PRSdf2.reset_index(drop = True)

#loop to get bottom and top files, while having a 1:1 balance between cases and controls, if the threshold is odd, it will remove one sample to balance the c:c
for ind in range(lenOfThresholds):

  id_bottom = 0 #iterates through the samples
  id_counter_bottom = 0 # tells how many samples are actually written into the file

  case_counter_bottom = 0 #how many cases are written into the ile
  control_counter_bottom = 0 #how many controls are written into the file

  file = open(str(individualThresholds[ind])+"_PRSids.txt", 'w')

  #add header to file
  list(COVdf.columns)
  header_string = "\t".join(map(str, COVdf.columns)) + "\n"
  file.write(header_string)

  #while amount of bottom samples in the file is less than the threshold
  while id_counter_bottom < individualThresholds[ind]:
      
      #need to optimize this, long for no reason; this should iterate through the IDs; for top
      IDbottom = sortedPRSdf['IID'][id_bottom]
      row = COVdf[COVdf['IID'] == IDbottom].index
      stringd_row = str(row)
      numRow = int(stringd_row[7:-17])

      try:
          #checker for phenotype 
          float(COVdf.loc[numRow]['PHENO'])

          #control
          if float(COVdf.loc[numRow]['PHENO']) == 1.0:
              if control_counter_bottom <= individualThresholds[ind]//2: ##check to make sure control to case ratio isn't overshot
                  row_string = "\t".join(map(str, COVdf.loc[numRow].values)) + "\n"
                  
                  file.write(str(row_string))
                  
                  control_counter_bottom += 1 #indicate you wrote in a control
                  id_counter_bottom += 1 #indicate you wrote in an id
                  id_bottom += 1 #iterate
              
              else:
                  id_bottom += 1 #iterate for if the counter is overshot

  
          #case
          elif float(COVdf.loc[numRow]['PHENO']) == 2.0:
              if case_counter_bottom <= individualThresholds[ind]//2:
                  row_string = "\t".join(map(str, COVdf.loc[numRow].values)) + "\n"
                  
                  file.write(str(row_string))
                  
                  case_counter_bottom += 1
                  id_counter_bottom += 1
                  id_bottom += 1

              else:
                  id_bottom += 1

      except:
          id_bottom += 1

      ##to get top PRS

  id_top = 0
  id_counter_top = 0

  case_counter_top = 0
  control_counter_top = 0

  while id_counter_top < individualThresholds[ind]:
      
      IDtop = sortedPRSdf_reverse['IID'][id_top]
      row1 = COVdf[COVdf['IID'] == IDtop].index
      stringd_row1 = str(row1)
      numRow1 = int(stringd_row1[7:-17])

      try:
          #checker for phenotype 
          float(COVdf.loc[numRow1]['PHENO'])

          #control
          if float(COVdf.loc[numRow1]['PHENO']) == 1.0:
              if control_counter_top <= individualThresholds[ind]//2:
                  row_string1 = "\t".join(map(str, COVdf.loc[numRow1].values)) + "\n"
                  
                  file.write(str(row_string1))
                  
                  control_counter_top += 1

                  id_counter_top += 1
                  id_top += 1
              
              else:
                  id_top += 1

              
          #case
          elif float(COVdf.loc[numRow1]['PHENO']) == 2.0:
              if case_counter_top <= individualThresholds[ind]//2:
                  row_string = "\t".join(map(str, COVdf.loc[numRow1].values)) + "\n"
                  
                  file.write(str(row_string))
                  
                  case_counter_top += 1

                  id_counter_top += 1
                  id_top += 1

              else:
                  id_top += 1

      except:
          id_top += 1

file.close()

filenames= []
#at this point, you have a file of all of the covariates that are at the top or bottom thresholds in threshold_PRS.txt
#there isn't really an extraction process, as skat and cmc will only analyze variants with an assigned phenotype
script_log = open("PRS_log.sh", 'w')
num_scripts = 0
for length in range(lenOfThresholds): #index for threshold list
  #you want to create a skat file for each
  skat = open(str(individualThresholds[length])+"_SkatO.sh", 'w')
  sbatchWritingMid(skat)

  skat.write(
      "/projects/b1049/genetics_programs/Rvtest/rvtests/executable/rvtest --inVcf " + iVCF + " --pheno " + str(individualThresholds[length])+"_PRSids.txt --pheno-name PHENO --covar "+ str(individualThresholds[length])+"_PRSids.txt --out "+ str(individualThresholds[length])+"_rareBurden.sh --covar-name " + covariateNames + " --kernel skato" 
  )
  skat.close()

  script_log.write(
      "PRS_"+str(individualThresholds[length])+"=($(sbatch "+str(individualThresholds[length])+"_rareBurden.sh))\n"
      "echo "'"PRS_'+str(individualThresholds[length])+' ${PRS_' +str(individualThresholds[length])+'[-1]}"'" >> PRS_ids\n"
      "\n"
  )

script_log.close()
sys.exit()
