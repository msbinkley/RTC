
import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
#param_file = "param_files/ryo_local.csv"
#param_file = "param_files/binkley_local.csv" 

param_file = "param_files/trial_local.csv"
paramDict = LP.load_param_file(param_file)


#ajay added --snps-only to all plink commands in allstepslist.py

#This step is incredibly fast. No need to optimize.
def run_all_steps_H0(paramDict, chrNum, uuidvalstr, prep_uuidstr):  
    #This selects a 'random' eQTL and GWAS within the selected coldspot
    AS.select_rand_variant(paramDict["gtexDir"] + "eQTLcausal" + uuidvalstr + ".txt", paramDict, chrNum,uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant(paramDict["gtexDir"] + "eQTLcausalpos" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)  
    AS.select_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausal" + uuidvalstr + ".txt", paramDict, chrNum,uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausalpos" + uuidvalstr + ".txt", paramDict,uuidvalstr, prep_uuidstr)
    
    #This identifies the linked eQTL
    #os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/eQTLcausal.txt --tag-r2 0.5")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex_copy2.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --tag-r2 0.5")  
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --tag-r2 0.5") 
    #now we will subset it
    cmd_plink="plink --vcf  " + paramDict["gtexDir"]+ "chr"+ str(chrNum) + "_GTEx_Analysis_v7_WGS.vcf.gz --show-tags " + paramDict["gtexDir"] + "eQTLcausal" + uuidvalstr + ".txt --snps-only --tag-r2 0.5 --out " + paramDict["gtexDir"]+ "plink"+ uuidvalstr
    os.system(cmd_plink)
    #import shutil
    #shutil.copy(paramDict["gtexDir"]+ "plink.tags", paramDict["gtexDir"]+ "plink"+ uuidvalstr + ".tags")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GTEx_Analysis_v7_WGS.vcf.gz --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --snps-only --tag-r2 0.5")
    AS.filter_genov_output(paramDict["gtexDir"] + "filteredSNPs" + uuidvalstr + ".txt", paramDict, prep_uuidstr)
    AS.filter_plink_output(paramDict["gtexDir"] + "filteredplink" + uuidvalstr + ".txt", paramDict,uuidvalstr, prep_uuidstr) 
    AS.select_rand_variant2(paramDict["gtexDir"] + "linkedeQTL" + uuidvalstr + ".txt", paramDict,uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos" + uuidvalstr + ".txt", paramDict,uuidvalstr, prep_uuidstr)   
    
    #commented below line out and doing it in RTC
    #AS.prepare_perm_file( "/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations.txt")
    
    #This now identifies the linked GWAS variant
    #os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/GWAScausal.txt --tag-r2 0.5")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex_copy2.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GWAScausal.txt --tag-r2 0.5")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GWAScausal.txt --tag-r2 0.5")
    #now for the subset
    cmd_plink2="plink --vcf "+ paramDict["gtexDir"]+ "chr"+ str(chrNum) + "_GTEx_Analysis_v7_WGS.vcf.gz --show-tags " + paramDict["gtexDir"] + "GWAScausal" + uuidvalstr + ".txt --snps-only --tag-r2 0.5 --out " + paramDict["gtexDir"]+ "plink"+ uuidvalstr
    print("command for plink : cmd_plink2 is : ", cmd_plink2)
    os.system(cmd_plink2)
    #import shutil
    #shutil.copy(paramDict["gtexDir"]+ "plink.tags", paramDict["gtexDir"]+ "plink"+ uuidvalstr + ".tags")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GTEx_Analysis_v7_WGS.vcf.gz --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GWAScausal.txt --snps-only --tag-r2 0.5")
    AS.filter_plink_output_gwas(paramDict["gtexDir"] + "filteredplink" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    AS.select_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWAS" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWASpos" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
#run_all_steps_H0(paramDict)

def run_all_steps_H1(paramDict, chrNum, uuidvalstr, prep_uuidstr):  
    #This selects a 'random' eQTL and GWAS within the selected coldspot
    AS.select_rand_variant(paramDict["gtexDir"] + "eQTLcausal" + uuidvalstr + ".txt", paramDict, chrNum, uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant(paramDict["gtexDir"] + "eQTLcausalpos" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)  
    AS.select_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausal" + uuidvalstr + ".txt", paramDict,chrNum, uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausalpos" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    #This identifies the linked eQTL
    #os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/eQTLcausal.txt --tag-r2 0.5")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex_copy2.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --tag-r2 0.5")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex.vcf --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --tag-r2 0.5")
    
    #now we will subset it
    cmd_plink="plink --vcf  "+ paramDict["gtexDir"] + "chr"+ str(chrNum) + "_GTEx_Analysis_v7_WGS.vcf.gz --show-tags " + paramDict["gtexDir"] + "eQTLcausal" + uuidvalstr + ".txt --snps-only --tag-r2 0.5 --out " + paramDict["gtexDir"]+ "plink"+ uuidvalstr
    os.system(cmd_plink)
    #import shutil
    #shutil.copy(paramDict["gtexDir"]+ "plink.tags", paramDict["gtexDir"]+ "plink"+ uuidvalstr + ".tags")
    #os.system("plink --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GTEx_Analysis_v7_WGS.vcf.gz --show-tags /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/eQTLcausal.txt --snps-only --tag-r2 0.5")
    AS.filter_genov_output(paramDict["gtexDir"] + "filteredSNPs" + uuidvalstr + ".txt", paramDict, prep_uuidstr)
    AS.filter_plink_output(paramDict["gtexDir"] + "filteredplink" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr) 
    AS.select_rand_variant2(paramDict["gtexDir"] + "linkedeQTL" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos" + uuidvalstr + ".txt", paramDict,uuidvalstr, prep_uuidstr)   
    
    #commented below line out and doing it in RTC
    #AS.prepare_perm_file( "/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations.txt")
    
    
    #This now identifies the linked GWAS variant
    #os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy.vcf --show-tags /users/michaelbinkley/desktop/GTEx/GWAScausal.txt --tag-r2 0.5")
    AS.filter_plink_output_gwas(paramDict["gtexDir"] + "filteredplink" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    AS.select_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWAS" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
    AS.output_pos_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWASpos" + uuidvalstr + ".txt", paramDict, uuidvalstr, prep_uuidstr)
#run_all_steps_H1(paramDict)

