#!/usr/bin/env python3

import gzip, os, subprocess 


def get_exp_vector(geneName, tissue, gtexDir):
    '''
    Output a vector of expression levels for a given gene and tissue 

    Note:  Warning, the order of this expression vector may not be the same as the genotype vector.
    '''

    #filePath = gtexDir + "/" + tissue + ".v7.normalized_expression.bed.gz"
    #edited to this 
    filePath = gtexDir + "egene_and_expbed_files/GTEx_Analysis_v7_eQTL_expression_matrices/" + tissue + ".v7.normalized_expression.bed.gz"
    fileIN = gzip.open(filePath)
    counter = 0
    matchingLine = ""
    header = fileIN.readline().decode().rstrip().split('\t')
    individualNames = header[4:]
    for i in fileIN:
        i =i.decode().rstrip().split("\t")
        if i[3]==geneName:
            matchingLine = i
            break
    fileIN.close()
    expVector = matchingLine[4:]
    expVector = [float(x) for x in expVector]
    return expVector, individualNames


def get_dos_vector(chrNum, pos, vcfDir):
    '''
    Outputs the dosage vector for a particular chromoosome and position. 
    
    #Make sure file has a .tbi. 
    If it doesn't, then do the following: 
    DN527o9v:vcfDir ryosukekita$ gunzip chr22_subset_gtex.vcf.gz
    DN527o9v:vcfDir ryosukekita$ bgzip chr22_subset_gtex.vcf
    DN527o9v:vcfDir ryosukekita$ tabix -p vcf chr22_subset_gtex.vcf.gz
    DN527o9v:vcfDir ryosukekita$ tabix chr22_subset_gtex.vcf.gz

    TEST SNP:  40051275        22_40051275_G_GC_b37  
    '''
    print("\tGetting dosage vector for :", chrNum, pos, vcfDir)
    #vcfFilePath= vcfDir + "/chr" + chrNum + "_subset_gtex.vcf.gz"
    #changed to this for chr22
    #vcfFilePath= vcfDir + "chr" + chrNum + "_subset_gtex.vcf.gz"
    #for full vcf 
    #vcfFilePath= vcfDir + "GTEx_Analysis_v7_WGS.vcf.gz"
    #for subset 
    vcfFilePath= vcfDir + "chr"+ str(chrNum)+ "_GTEx_Analysis_v7_WGS.vcf.gz"
    print('\n this is vcf path', vcfFilePath)
    command = ["tabix", vcfFilePath ,  str(chrNum) + ":" + str(pos) + "-" + str(pos)]
    print('command is ', command)
    proc = subprocess.Popen(command, stdout = subprocess.PIPE)
    output, err = proc.communicate()
    if len(output)==0:
        print("Error: There is a problem with the tabix output. Check to make sure the vcf file is tabix-indexed. May need to reindex. Exiting. ")
        #sys.exit()
    #commented below out for lines after to remove non SNPs and concatenate multiple potential lists together if there are two snps in same location
    #genotypes = output.decode().split("\t")[9:]
    
    genotypesbyn = output.decode().split("\n")
    countg=0
    genotypes=[]
    for i in genotypesbyn:
        print(i)
        genotype_test = genotypesbyn[countg].split("\t")
        countg=countg+1
        if len(genotype_test) > 1:
            if len(genotype_test[3]) + len(genotype_test[4]) == 2:
                genotype_forlist= genotype_test[9:]
                #genotypes= genotypes + genotype_forlist # i commented this out for now but have to come back and change later to include multiple snp changes, use gtex position : chr 12, 12:111754704-111754704
                genotypes= genotype_forlist 
    
    #have to do this command below for gtex vcf that has ./. for unannotated 
    #print("gentoypes", genotypes)
    count = 0
    indextodrop=[]
    for x in genotypes:
        if '.' in x:
            print(x)
            indextodrop.append(count)
        count = count+1
    print("indextodrop from within gtex", indextodrop)
    genotypes = [x for x in genotypes if "." not in x ] #(here we are removing them)
    #genotypes = [x if "." not in x else '0/0' for x in genotypes] #use this command instead 
    #print("Genotypes:", genotypes)
    genotypes = [x.split(":")[0] for x in genotypes]
    print("Genotypes again:", genotypes)
    dosageSplit = [[int(y) for y in x.split("/")]    for x in genotypes]
    dosage =  [sum(x)    for x in dosageSplit]
    
    print("\tGetting header")
    command = ["tabix", vcfFilePath ,  "-H", str(chrNum) + ":"  + str(pos) + "-" + str(pos)]
    print('tabix command for getting indivs', command)
    print("\tStillok")
    proc = subprocess.Popen(command, stdout = subprocess.PIPE)
    #ajay
    #print('proc', proc)
    output, err = proc.communicate()
    #print('output', output)
    print('err', err)
    indivs = output.decode().split("\n")[-2].split("\t")[9:]
    namestodelete=[indivs[i] for i in indextodrop]
    #added this instead of names to drop because i think names to drop is wrong
    for i in reversed(indextodrop):
        del indivs[i]
   # print("dosage:", dosage) 
    #i added indextodrop
    return dosage, indivs, indextodrop, namestodelete


def filter_and_sort_genotype_vector(genoV, indivGenoV, indivExpV, indextodrop, namestodelete):
    '''
    Filters and sorts the genotype Vector so that the individuals matching the indivExpVector

    '''
    #ajay if needed: use this command to match if exists : s_indx = [s.index(x) for x in s2 if x in s]
    #what i can do is send index of the genotype with dot here and remove it from GenoV and IndivExpV
    
    #find and drop the indexes that had no expression or ./. for expression
    #first get the names so you can delete from indivExpV
    print('len of genoV is ', len(genoV))
    print('len of indivGenoV is ', len(indivGenoV))
    print('len of indivExpV is ', len(indivExpV))
    #namestodelete=[indivGenoV[i] for i in indextodrop]
    indivGenoV2 = [x for x in indivGenoV if x not in namestodelete]
    indextodropforexp = [indivExpV.index(x) for x in indivExpV if x in namestodelete]
    #indivExpV2 = [x for x in indivExpV if x not in namestodelete]
    indextokeepforexp = [indivExpV.index(x) for x in indivExpV if x in indivGenoV2]
    #print("index to keep", indextokeepforexp)
    #print("index to keep len", len(indextokeepforexp))
    #print('len of indivExpV at first is ', len(indivExpV))
    #added this line
    indivExpV2 = [x for x in indivExpV if x in indivGenoV2]
    #print('len of indivGenoV2 is ', len(indivGenoV2))
    #print('len of indivExpV2 second is ', len(indivExpV2))
    #indivGenoVIdx = [indivGenoV.index(x) for x in indivExpV]
    #changed above to this 
    indivGenoVIdx = [indivGenoV2.index(x) for x in indivExpV2]
    print('indivGenoVIdx is ', indivGenoVIdx)
    print('indivGenoV2 ', indivGenoV2)
    print('genoV is ', genoV)    
    print('indextodropforexp is', indextodropforexp)
    #print('indivGenoVIdx is ', indivGenoVIdx)
    #print('indivGenoV', indivGenoV)
    print('indivExpV2 is', indivExpV2)
    filteredAndSortedGenoV = [genoV[x]  for x in indivGenoVIdx]
    #return filteredAndSortedGenoV, indivExpV
    return filteredAndSortedGenoV, indivExpV2, indextodropforexp, indextokeepforexp


#added this function for def calculate pp
def filter_and_sort_genotype_vector2(genoV, indivGenoV, indivExpV, indextodrop):
    '''
    Filters and sorts the genotype Vector so that the individuals matching the indivExpVector

    '''
    #ajay if needed: use this command to match if exists : s_indx = [s.index(x) for x in s2 if x in s]
    #what i can do is send index of the genotype with dot here and remove it from GenoV and IndivExpV
    
    #find and drop the indexes that had no expression or ./. for expression
    #first get the names so you can delete from indivExpV
    print('len of genoV is ', len(genoV))
    print('len of indivGenoV is ', len(indivGenoV))
    print('len of indivExpV is ', len(indivExpV))
    namestodelete=[indivGenoV[i] for i in indextodrop]
    indivGenoV2 = [x for x in indivGenoV if x not in namestodelete]
    indextodropforexp = [indivExpV.index(x) for x in indivExpV if x in namestodelete]
    #indivExpV2 = [x for x in indivExpV if x not in namestodelete]
    indextokeepforexp = [indivExpV.index(x) for x in indivExpV if x in indivGenoV2]
    print("index to keep", indextokeepforexp)
    print("index to keep len", len(indextokeepforexp))
    print('len of indivExpV at first is ', len(indivExpV))
    #added this line
    indivExpV2 = [x for x in indivExpV if x in indivGenoV2]
    print('len of indivGenoV2 is ', len(indivGenoV2))
    print('len of indivExpV2 second is ', len(indivExpV2))
    #indivGenoVIdx = [indivGenoV.index(x) for x in indivExpV]
    #changed above to this 
    indivGenoVIdx = [indivGenoV2.index(x) for x in indivExpV2]
    print('genoV is ', genoV)    
    print('indextodropforexp is', indextodropforexp)
    print('indivGenoVIdx is ', indivGenoVIdx)
    print('indivGenoV', indivGenoV)
    print('indivExpV is', indivExpV)
    filteredAndSortedGenoV = [genoV[x]  for x in indivGenoVIdx]
    #return filteredAndSortedGenoV, indivExpV
    return filteredAndSortedGenoV, indivExpV2, indextodropforexp, indextokeepforexp








def filter_and_sort_genotype_vector4(genoV, indivGenoV, indivExpV, indextodrop):
    '''
    Filters and sorts the genotype Vector so that the individuals matching the indivExpVector

    '''
    #ajay if needed: use this command to match if exists : s_indx = [s.index(x) for x in s2 if x in s]
    #what i can do is send index of the genotype with dot here and remove it from GenoV and IndivExpV
    
    #find and drop the indexes that had no expression or ./. for expression
    #first get the names so you can delete from indivExpV
    print('len of genoV is ', len(genoV))
    print('len of indivGenoV is ', len(indivGenoV))
    print('len of indivExpV is ', len(indivExpV))
    namestodelete=[indivGenoV[i] for i in indextodrop]
    indivGenoV2 = [x for x in indivGenoV if x not in namestodelete]
    indextodropforexp = [indivExpV.index(x) for x in indivExpV if x in namestodelete]
    #indivExpV2 = [x for x in indivExpV if x not in namestodelete]
    indextokeepforexp = [indivExpV.index(x) for x in indivExpV if x in indivGenoV2]
    print("index to keep", indextokeepforexp)
    print("index to keep len", len(indextokeepforexp))
    print('len of indivExpV at first is ', len(indivExpV))
    #added this line
    indivExpV2 = [x for x in indivExpV if x in indivGenoV2]
    print('len of indivGenoV2 is ', len(indivGenoV2))
    print('len of indivExpV2 is ', len(indivExpV2))
    #indivGenoVIdx = [indivGenoV.index(x) for x in indivExpV]
    #changed above to this 
    indivGenoVIdx = [indivGenoV2.index(x) for x in indivExpV2]
    print('genoV is ', genoV)    
    print('indextodropforexp is', indextodropforexp)
    print('indivGenoVIdx is ', indivGenoVIdx)
    print('indivGenoV', indivGenoV)
    print('indivExpV is', indivExpV)
    filteredAndSortedGenoV = [genoV[x]  for x in indivGenoVIdx]
    #return filteredAndSortedGenoV, indivExpV
    return filteredAndSortedGenoV, indivExpV2, indextodropforexp, indextokeepforexp


#reference
def filter_and_sort_genotype_vector3(genoV, indivGenoV, indivExpV, indextodrop):
    '''
    Filters and sorts the genotype Vector so that the individuals matching the indivExpVector

    '''
    #ajay if needed: use this command to match if exists : s_indx = [s.index(x) for x in s2 if x in s]
    #what i can do is send index of the genotype with dot here and remove it from GenoV and IndivExpV
    
    #find and drop the indexes that had no expression or ./. for expression
    #first get the names so you can delete from indivExpV
    print('len of genoV is ', len(genoV))
    print('len of indivGenoV is ', len(indivGenoV))
    print('len of indivExpV is ', len(indivExpV))
    namestodelete=[indivGenoV[i] for i in indextodrop]
    indivGenoV2 = [x for x in indivGenoV if x not in namestodelete]
    indextodropforexp = [indivExpV.index(x) for x in indivExpV if x in namestodelete]
    indivExpV2 = [x for x in indivExpV if x not in namestodelete]
    #added this line
    indivExpV2 = [x for x in indivExpV if x in indivGenoV2]
    print('len of indivGenoV2 is ', len(indivGenoV2))
    print('len of indivExpV2 is ', len(indivExpV2))
    #indivGenoVIdx = [indivGenoV.index(x) for x in indivExpV]
    #changed above to this 
    indivGenoVIdx = [indivGenoV2.index(x) for x in indivExpV2]
    print('genoV is ', genoV)    
    print('indextodropforexp is', indextodropforexp)
    print('indivGenoVIdx is ', indivGenoVIdx)
    print('indivGenoV', indivGenoV)
    print('indivExpV is', indivExpV)
    filteredAndSortedGenoV = [genoV[x]  for x in indivGenoVIdx]
    #return filteredAndSortedGenoV, indivExpV
    return filteredAndSortedGenoV, indivExpV2, indextodropforexp, indextokeepforexp



