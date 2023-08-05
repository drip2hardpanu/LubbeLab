#workflow
#print all genes from 05 maf into $1, print p-values from 05 maf into $2
import pandas
 
L = []

 #file in this would be the 05 file
file = pandas.read_table("05/group3.SkatO.assoc")
f = open("final.txt", "a+")
f.write("Gene Name")
for x in file['Gene']:
    L.append(x)    
    f.write(x)
    f.write("\n")
f.close()

file05 = pandas.read_table("05/group3.SkatO.assoc")
f05 = open("05.txt", "a+")
#this is where i make a new file with just one column of p-values, if there is a p value it will show up in this file and if there isn't it'll just show up as NA

f05.write("MAF 05\n")

f05geneList = []

for geneName in file05['Gene']:
    f05geneList.append(geneName)

#this will look through the gene names in L, see if they're in the 01 gene list, and report their p values accordingly
for gene in L:
    if gene in f05geneList:
        phrase = ("Gene == "+"'"+ gene +"'" +"")
        pLocate = file05.query(phrase)["Pvalue"]
        P = pLocate.iloc[0]
        Q = str(P)
        f05.write(Q)
        f05.write("\n")
    else:
        f05.write("NA")
        f05.write("\n")
f05.close()

#file in this would be the 001 file
file01 = pandas.read_table("01/group3.SkatO.assoc")
f01 = open("01.txt", "a+")
#this is where i make a new file with just one column of p-values, if there is a p value it will show up in this file and if there isn't it'll just show up as NA

f01.write("MAF 01 P-Val\n")

f01geneList = []

for geneName in file01['Gene']:
    f01geneList.append(geneName)

#this will look through the gene names in L, see if they're in the 01 gene list, and report their p values accordingly
for gene in L:
    if gene in f01geneList:
        phrase = ("Gene == "+"'"+ gene +"'" +"")
        pLocate = file01.query(phrase)["Pvalue"]
        P = pLocate.iloc[0]
        Q = str(P)
        f01.write(Q)
        f01.write("\n")
    else:
        f01.write("NA")
        f01.write("\n")
f01.close()

#file in this would be the 001 file
file001 = pandas.read_table("001/group3.SkatO.assoc")
f001 = open("001.txt", "a+")
#this is where i make a new file with just one column of p-values, if there is a p value it will show up in this file and if there isn't it'll just show up as NA

f001.write("MAF 01 P-Val\n")

f001geneList = []

for geneName in file001['Gene']:
    f001geneList.append(geneName)

#this will look through the gene names in L, see if they're in the 01 gene list, and report their p values accordingly
for gene in L:
    if gene in f001geneList:
        phrase = ("Gene == "+"'"+ gene +"'" +"")
        pLocate = file001.query(phrase)["Pvalue"]
        P = pLocate.iloc[0]
        Q = str(P)
        f001.write(Q)
        f001.write("\n")
    else:
        f001.write("NA")
        f001.write("\n")
f001.close()

#look for those gene names in 01 maf, if there, print p-value into $1 of a new file
#^do that for all of the maf threshholds
#paste each of those files onto the original file, going from 01 - 0001python -m IPython

file0001 = pandas.read_table("0001/group3.SkatO.assoc")
f0001 = open("0001.txt", "a+")
#this is where i make a new file with just one column of p-values, if there is a p value it will show up in this file and if there isn't it'll just show up as NA

f0001.write("MAF 01 P-Val\n")

f0001geneList = []

for geneName in file0001['Gene']:
    f0001geneList.append(geneName)

#this will look through the gene names in L, see if they're in the 01 gene list, and report their p values accordingly
for gene in L:
    if gene in f0001geneList:
        phrase = ("Gene == "+"'"+ gene +"'" +"")
        pLocate = file0001.query(phrase)["Pvalue"]
        P = pLocate.iloc[0]
        Q = str(P)
        f0001.write(Q)
        f0001.write("\n")
    else:
        f0001.write("NA")
        f0001.write("\n")
f0001.close()


fileSingletons = pandas.read_table("singletons/group3.SkatO.assoc")
fSingletons = open("singletons.txt", "a+")
#this is where i make a new file with just one column of p-values, if there is a p value it will show up in this file and if there isn't it'll just show up as NA

fSingletons.write("Singletons\n")

fSingletonsgeneList = []

for geneName in fileSingletons['Gene']:
    fSingletonsgeneList.append(geneName)

#this will look through the gene names in L, see if they're in the 01 gene list, and report their p values accordingly
for gene in L:
    if gene in fSingletonsgeneList:
        phrase = ("Gene == "+"'"+ gene +"'" +"")
        pLocate = fileSingletons.query(phrase)["Pvalue"]
        P = pLocate.iloc[0]
        Q = str(P)
        fSingletons.write(Q)
        fSingletons.write("\n")
    else:
        fSingletons.write("NA")
        fSingletons.write("\n")
fSingletons.close()