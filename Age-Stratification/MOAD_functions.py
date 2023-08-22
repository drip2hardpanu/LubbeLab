import pandas
import os
#notes: ensure age is titled AGE in covariate file
#eoad < 60, moad < 70, LOAD > 70

def age_strattingPANDAS(prefix, inputCovar):
    mo = open(prefix+"_mid_onset.txt", 'w')
    lo = open(prefix+"_late_onset.txt", 'w')

    covdf = pandas.read_table(inputCovar)
    covariatedf_head = list(covdf.columns)

    #early onset
    eo_df = covdf[covdf["AGE"] < 60]
    eo_df.to_csv(prefix+"_early_onset.txt")

    #mid onset
    mo_df = covdf[covdf["AGE"] < 70]
    mo_df.to_csv(prefix+"_mid_onset.txt")

    #late onset
    lo_df = covdf[covdf["AGE"] >= 70]
    lo_df.to_csv(prefix+"_late_onset.txt")
