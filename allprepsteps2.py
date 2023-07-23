#####This will be the iterative step, run this for a given eQTL, you pick the best GWAS and run it for H0
# The input eQTL is under /users/michaelbinkley/desktop/RTCstuffs/besteQTLs.txt

import os
import load_params as LP
import prepstepsCopy2 as PS
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
def prep_step(paramDict, geneName, chrNum, pos, tissue, snp, pvalue):
    #This step takes the given input eQTL from above and outputs the cold spot coordinates and the best GWAS
    #ajay added snp and pvalue to all functions
    PS.find_variant_coldspot('eQTLcoldspot.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    PS.output_selected_coldspots('selectedcoldspots.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    #PS.find_variant_coldspot_gwas("gwasColdSpot.txt", paramDict, geneName, chrNum, pos, tissue, snp, pvalue) #dont have to run this because it is called within sortedGWAS.txt
    PS.sort_combined_gwas('sortedGWAS.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    if os.stat("/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/sortedGWAS.txt").st_size != 0:
        PS.output_header("bestGWAS.txt")    
        PS.output_selected_coldspots('preVCF.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
        PS.filter_vcf('filteredgenotype.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
#prep_step(paramDict, geneName, chrNum, pos, tissue)
