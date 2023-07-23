#####This will be the iterative step, run this for a given eQTL, you pick the best GWAS and run it for H0
# The input eQTL is under /users/michaelbinkley/desktop/RTCstuffs/besteQTLs.txt

import os
import load_params as LP
import prepstepsCopy1_null as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
#param_file = "param_files/ryo_local.csv"
#param_file = "param_files/binkley_local.csv" 
param_file = "param_files/trial_local.csv" 

paramDict = LP.load_param_file(param_file)



#geneName = "ENSG00000172404.4"
#chrNum = "22"
#pos = "40407887.0"
#tissue = "Breast_Mammary_Tissue"
#this step works and takes ~2 min to run
def prep_step(paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr, gwasfilename):
    #This step takes the given input eQTL from above and outputs the cold spot coordinates and the best GWAS
    #ajay added snp and pvalue to all functions
    PS.find_variant_coldspot('eQTLcoldspot' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    PS.output_selected_coldspots('selectedcoldspots' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr)
    #PS.find_variant_coldspot_gwas("gwasColdSpot.txt", paramDict, geneName, chrNum, pos, tissue, snp, pvalue) #dont have to run this because it is called within sortedGWAS.txt
    PS.sort_combined_gwas('sortedGWAS' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr, gwasfilename)
    if os.stat("/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/sortedGWAS" + prep_uuidstr + ".txt").st_size != 0:
        PS.output_header("bestGWAS" + prep_uuidstr + ".txt", prep_uuidstr)    
        PS.output_selected_coldspots('preVCF' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr)
        PS.filter_vcf('filteredgenotype' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr)
#prep_step(paramDict, geneName, chrNum, pos, tissue)
