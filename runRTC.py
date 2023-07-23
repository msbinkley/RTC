import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy
from scipy import stats
import sys
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import sys
import codecs
 
import random
random.seed(0)
from random import sample

import random
random.seed(0)
from random import sample

import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import allstepslist as ASL


def run_RTC(paramDict, geneName, chrNum, pos, tissue, uuidvalstr, prep_uuidstr):
    ### RUn this for the real eQTL
    import load_params as LP
    import eqtl_funcs as EF
    import allstepsH as AS
    import numpy as np
    
    #geneName = "ENSG00000172404.4"
    #chrNum = "22"
    #pos = "40407887"
    #tissue = "Breast_Mammary_Tissue"

    print("Start") #SLOW
    slope, intercept,  residuals, indivVector, indextodropforexp = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  
    print("RTC1")

    
    def calculate_pp(slope_r, intercept, residuals, chrNum, pos, indivVector, paramDict, indextodropfromindiv):
        import numpy as np
        permutedResiduals = random.sample(list(residuals), len(residuals))
        genoVector, indivGenoVector,indextodrop,namestodelete = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
        #need to do reversed because if you delete then it will move up index 
        for i in reversed(indextodropfromindiv):
            del genoVector[i]
        print("len of indivVector ", len(indivVector))
        print("index to drop", indextodrop)
        print("index to drop from indiv", indextodropfromindiv)
        print("names to delete ", namestodelete)
        for i in reversed(indextodropfromindiv):
            del indivGenoVector[i]
        #for i in reversed(indextodrop):
        #    del indivVector[i]
        print("inside pp")
        print("this is slope_r", slope_r)
        print("len of genoVector ", len(genoVector))
        print("len of indivGenoVector ", len(indivGenoVector))
        print("len of indivVector ", len(indivVector))
        filteredGenoVector, indivVector, indextodropforexp, indextokeepforexp = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivVector, indextodrop, namestodelete)
        permutedResiduals2=[permutedResiduals[x] for x in indextokeepforexp]
        for i in reversed(indextodropforexp):
            del permutedResiduals[i]
        print("permutedResiduals ",permutedResiduals)
        print("permutedResiduals len ", len(permutedResiduals))
        print("permutedResiduals2 ", permutedResiduals2)
        print("permutedResiduals2 ", len(permutedResiduals2))
        pp = float(slope_r)*np.array(filteredGenoVector) + permutedResiduals2 + intercept   
        return pp



   # random.seed(0)

    #param_file = "param_files/binkley_local.csv" 
    param_file = "param_files/trial_local.csv" 
    
    paramDict = LP.load_param_file(param_file)
    print("RTC2")
    #geneName = "ENSG00000172404.4"
    #chrNum = "22"
    #pos = "40407887"
    #tissue = "Breast_Mammary_Tissue"

    #print(pos)
    slope_r, intercept,  residuals, indivVector, indextodropfromindiv = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  
    print("RTC3")

    chrNumLinked = chrNum
    ### THIS is the position of the randomly selected
    #posLinked = output_pos_rand_variant(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)
    posLinked = AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict, uuidvalstr, prep_uuidstr)
    print("RTC4") #SLOW
    print("pos is ", pos)
    print("position linked", posLinked, posLinked[0])
    posLinked = posLinked[0]
    #old
    #pp = calculate_pp(slope_r, intercept, residuals, chrNumLinked, posLinked, indivVector, paramDict)
    #ajay 
    pp = calculate_pp(slope_r, intercept, residuals, chrNumLinked, posLinked, indivVector, paramDict, indextodropfromindiv)
    print("RTC5")

    def output_bed(indivVector, expVector, filepath, geneName, paramDict):
        selected = output_gene_info(geneName, paramDict)
        header = output_header2(indivVector)
        with open(filepath, 'w') as outfile:
            outfile.write(header + "\n")
            geneLocation = selected
            outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expVector]) + "\n")
    #print(geneName)        
    output_bed(indivVector, pp, paramDict["tmpSmall"] + "/GTExExpfile" + uuidvalstr + ".bed", geneName, paramDict)
    
    #added this below, need the file above for permutations.txt, I am uncommenting it out from ASL or allstepslist
    AS.prepare_perm_file( "/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations" + uuidvalstr + ".txt", uuidvalstr, prep_uuidstr)
    print("RTC6") #SLOW

    #This prepares the index files for the pseudophenotype
    #print("RTC6A")
    #os.system("bgzip /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed")
    #print("RTC6B")
    #os.system("tabix -p bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz")
    #print("RTC6C")
    #This calculates the RTC score
    #os.system("QTLtools rtc --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex.vcf.gz --bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz --hotspot /users/michaelbinkley/desktop/qtltools_test/hotspots.bed --gwas-cis /users/michaelbinkley/desktop/GTEx/linkedGWAS.txt /users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations.txt --normal --out rtc_results2.txt")
    
    
    #copied below for sherlock
    #This prepares the index files for the pseudophenotype
    print("RTC6A")
    os.system("bgzip -f -c /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile" + uuidvalstr + ".bed > /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile" + uuidvalstr + ".bed.gz")
    print("RTC6B")
    os.system("tabix -p bed /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile" + uuidvalstr + ".bed.gz")
    print("RTC6C")
    #This calculates the RTC score
    #os.system("QTLtools_1.2_CentOS7.8_x86_64 rtc --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/GTEx_Analysis_v7_WGS.vcf.gz --bed /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile.bed.gz --hotspot /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/hotspots.bed --gwas-cis /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/linkedGWAS.txt /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations.txt --normal --out rtc_results2.txt")
    vcf_file_path = "/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/chr"+str(chrNum)+"_GTEx_Analysis_v7_WGS.vcf.gz"
    cmd_forqtl = "QTLtools_1.2_CentOS7.8_x86_64 rtc --vcf %s" %(vcf_file_path) + " --bed /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile" + uuidvalstr + ".bed.gz --hotspot /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/hotspots.bed --gwas-cis /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/linkedGWAS" + uuidvalstr + ".txt /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations" + uuidvalstr + ".txt --normal --out rtc_results2" + uuidvalstr + ".txt" 
    os.system(cmd_forqtl)
   
    
    #print("RTC6A")
    #os.system("bgzip /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile_realeQTL.bed")
    #print("RTC6B")
    #os.system("tabix -p bed /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile_realeQTL.bed.gz")
    #os.system("QTLtools rtc --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex.vcf.gz --bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile_realeQTL.bed.gz --hotspot /users/michaelbinkley/desktop/qtltools_test/hotspots.bed --gwas-cis /users/michaelbinkley/desktop/RTCstuffs/bestGWAS_var.txt /users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations_real.txt --normal --out rtc_results_realeQTL.txt") 
    #print("RTC6C")
    #os.system("QTLtools_1.2_CentOS7.8_x86_64 rtc --vcf /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/chr22_subset_gtex.vcf.gz --bed /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile_realeQTL.bed.gz --hotspot /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/qtltools_test/hotspots.bed --gwas-cis /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/bestgwas_var.txt /oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/permutations_real.txt --normal --out rtc_results_realeQTL.txt")

    print("RTC7")
#run_RTC(paramDict)

def output_gene_info(geneName, paramDict):
    
    fileIN = open(paramDict["tmpDir"] + '/chrname.txt')
    selected = list()
    for i in fileIN: 
        i = i.rstrip().split('\t')
        
        if i[3]== geneName:
            selected = i
            break
            
    return selected

  
def output_header2(indVector): 

            
            
    header = "\t".join(['#chr',  'start', 'end' , 'gene' , 'length' , 'strand'] + [str(x) for x in indVector] )
    return header

    fileIN.close()

