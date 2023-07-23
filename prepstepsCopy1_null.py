import os
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import load_params as LP
#param_file = "param_files/ryo_local.csv"
#param_file = "param_files/binkley_local.csv" 
param_file = "param_files/trial_local.csv" 
paramDict = LP.load_param_file(param_file)


def ID_hotspots(outfilename): 
    fileIN = open(paramDict["gtexDir"] + 'coldspots.txt')
    hotspots = list()
    #fileOUT = open(outfilename, "w")
    for i in fileIN: 
        i = i.rstrip().split('\t')

        
        Chr = i[0]
        for num in i[0]: #Get rid of "chr in the file"
            if num in "chr":
                Chr = Chr.replace(num, '')  
        start = int(i[1])
        end = int(i[2])
        middle = str(i[3])
        bound = str(i[4])
        #need to do this for coldspots.txt in sherlock 
        #size = len(bound)
        # Slice string to remove last 3 characters from string
        #bound = bound[:size - 2] 

        hotspots.append([Chr, start, end, middle, bound])
        #hotspots.append(start)
        #hotspots.append(end)
        #hotspots.append(middle)
        #fileOUT.write(Chr + "\t" + str(start) + "\t" + str(end) + "\t" + str(middle)  +  "\t" + str(bound) + "\n")
    fileIN.close()
    #fileOUT.close()
    return hotspots

#####OUtputs cold spot for given eqtl
def find_variant_coldspot(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue):
    #fileIN = open( "/users/michaelbinkley/desktop/RTCstuffs/besteQTLs.txt" )
    
    fileOUT = open(outfilename, "w")
    hotspots = ID_hotspots('coldspots2.txt')
    eQTLhotspot = list()
    #print(hotspots[:1])
    #print(hotspots[:10])
    


    #for i in fileIN:
        
        #i = i.rstrip().split('\t')
        #gene = str(i[0])
        #chrnum = str(i[1])
        #snp = (i[3])
        #position = float(i[2])
        #pvalue = float(i[4])
    count1=0
    for  Chr, start, end, middle, bound in hotspots:
        if Chr == chrNum:
            count1 = count1+1
    count2=0
    somevalue=1
    for  Chr, start, end, middle, bound in hotspots:
        count2=count2+1
        if Chr  != chrNum :
            continue
        else: 
            #if float(pos) < float(start) and somevalue == 1: #add this for start 
            #    eQTLhotspot.append([Chr, 0, start])
            #    fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , geneName, snp, str(pos), str(pvalue)] ) + "\n")
            #    somevalue =2 #change count so that this only happens once 
            #print('pos is :', float(pos), 'start is :', float(start), 'bound is :', float(bound))
            if float(pos) < float(bound) and float(pos) > float(end): 
            #if float(pos) < float(bound) and float(pos) > float(start): #this ajay changed, ask Dr. Binkley because there are some that are not within start and end
                    #distance = abs(position - middle)
                #eQTLhotspot.append([Chr, end, bound])
                eQTLhotspot.append([Chr, end, bound]) #ajay changed to this because some were in between start and end 
                #print('pos is ', pos, 'needs to be less than bound : ', bound, 'and more than end :', end)
                #fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , gene, snp, str(position), str(pvalue)] ) + "\n")
                fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , geneName, snp, str(pos), str(pvalue)] ) + "\n") #changed gene name
            #elif count2 == count1 and float(pos) > float(end): #this is for end of chromosome but dont think i need it
            #    eQTLhotspot.append([Chr, end, bound])
            #    fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , geneName, snp, str(pos), str(pvalue)] ) + "\n")
            #    continue
            else: 
                continue
    #fileIN.close()
    fileOUT.close()
    return eQTLhotspot
    

#find_variant_coldspot('eQTLcoldspot.txt', paramDict, geneName, chrNum, pos, tissue)

#### This step can probably be truncated 
def output_selected_coldspots(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr):
    fileIN = open ('eQTLcoldspot' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue)
    
    fileOUT = open(outfilename, "w")
    selectedcoldspots = list()
    
    for i in fileIN:
        
        i = i.rstrip().split('\t')
        Chr = str(i[0])
        end = str(i[2])
        bound = str(i[4])
        print('we are in second method : output selected coldspots now')
    

        fileOUT.write("\t".join([ str(Chr),  str(end),  str(bound)] ) + "\n")
        selectedcoldspots.append([ str(Chr),  str(end),  str(bound)])                    
    
    fileIN.close()
    fileOUT.close()
    return selectedcoldspots
#output_selected_coldspots('selectedcoldspots.txt', paramDict, geneName, chrNum, pos, tissue)




# In[8]:

###THis outputs a file with GWAS variants in our cold spots of interest that colocalize. it will output the best GWAS
## first step takes 2 min to run
def find_variant_coldspot_gwas(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr, gwasfilename):
    #fileIN = open(paramDict["gtexDir"] +'/oncoarray_bcac_public_release_oct17.txt')
    #fileIN = open(paramDict["gtexDir"] +'/oncoarray_bcac_GTEXsnpsonly_public_release_oct17.txt')
    #fileIN = open(paramDict["gtexDir"] +'/RTCperm_null2_oncoarray_breast.txt')
    #fileIN = open(paramDict["gtexDir"] +'/RTCperm_int8_null_gtex_oncoarray_breast.txt')
    #fileIN = open(paramDict["gtexDir"] +'/RTCperm_null_gtex_oncoarray_colonsigmoid.txt')
    #instead of above we can put it in like this
    fileIN = open(gwasfilename)

    #for BreastGwas 
    #fileIN = open(paramDict["gtexDir"] +'BreastCancer_GWASlist_csv.txt')
    #for GwasCatalog, check r script in local to make this file
    #fileIN = open(paramDict["gtexDir"] +'CoronaryHeartGwasCatalogsubset.txt')
    
    #if using gtex null set 
    #fileIN = open(paramDict["gtexDir"] +'/RTCperm_int6_null_gtex_oncoarray_breast.txt')
    
    #fileIN.readline()
    fileOUT = open(outfilename, "w")

    cspots = find_variant_coldspot('eQTLcoldspot' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    #print("cspots is: ",cspots)
    #cspots = output_selected_coldspots('selectedcoldspots.txt')
    GWAScoldspots = list()
    for i in fileIN:
        #print(' line 1 in breastgwas read in file')
        #i = i.rstrip().split('\t')
        #gene = str(i[0])
        #chrnum = str(i[2])
        ##snp = (i[3]) #this was already commented out before
        #position = float(i[3])
        #pvalue = float(i[9])
        
        #below is for gtex null set of GWAS
        i = i.rstrip().split('\t')
        gene = str(i[7])[:-4]
        chrnum = str(i[5])
        position = float(int(i[6]))
        pvalue = .0000000000000001 #just giving fake pvalue 
        if geneName != str(i[2]):
            continue
        print("geneName is this : ", geneName, " and gene name of matched nullset is this : ", str(i[2]) )
        #below is for Coronary Heart Disease
        #i = i.rstrip().split(' ')
        #gene = str(i[2])
        #chrnum = str(i[0])
        #position = float(i[1])
        #pvalue = float(i[3])
        
        
        #below is for breast
        #i = i.rstrip().split(',')
        #gene = str(i[0])
        #chrnum = str(i[1])
        #position = float(i[2])
        #pvalue = float(i[7])
        if pvalue<10e-12:
            #print('pvalue is ', pvalue)
            for  (Chr), (end), (bound) in cspots:
                if Chr  != chrnum:
                    continue
                else: 
                    if position < float(bound) and position > float(end): 
                    #distance = abs(position - middle)
                        GWAScoldspots.append([Chr, end, bound, gene, position, pvalue])
                        #fileOUT.write("\t".join([ str(Chr), str(end),  str(bound) , gene, str(position), str(pvalue)] ) + "\n")
                    else: 
                        continue
    #print('this is Gwascoldspots :', GWAScoldspots)
    return GWAScoldspots
    fileIN.close()
    #fileOUT.close()
#find_variant_coldspot_gwas("gwasColdSpot.txt")

def sort_combined_gwas(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr, gwasfilename):
    #fileIN = open("gwasColdSpot.txt")
    GWAScoldspots = find_variant_coldspot_gwas("gwasColdSpot"+ prep_uuidstr + ".txt", paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr, gwasfilename)
    sortedData = list()
    fileOUT = open(outfilename, "w")
    #for i in fileIN:
    sortedlist = list()
    for i in GWAScoldspots:
        #i=i.rstrip().split("\t")
        i =[i[0], (i[1]), (i[2]), (i[3]), (i[4]), float(i[5])]
        sortedData.append( [(i[0]), (i[1]), (i[2]), i[3], (i[4]), (i[5])])
    sortedScoreData = sorted(sortedData, key = lambda x : (x[5]), reverse=False) #Sort by the element in the column index 
    #print('this is sortedScoreData', sortedScoreData)
    for i in sortedScoreData:
        #print(i)
        #sortedlist.append([(i[0]), (i[1]), (i[2]), i[3], (i[4]), (i[5])])
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")
    return sortedlist
    #return [x[0] for x in sortedScoreData]
    fileOUT.close()
    fileIN.close()
#sort_combined_gwas('sortedGWAS.txt')


### This outputs the best GWAS for a given coldspot
def output_header(outfilename, prep_uuidstr):
    print(outfilename)
    fileIN = open('sortedGWAS'+prep_uuidstr+'.txt')
    fileOUT = open(outfilename, "w")

    fileOUT.write(fileIN.readline() + "\n")
    #return list
    fileOUT.close()
    fileIN.close()
#output_header("bestGWAS.txt")


def prep_real_GWAS(outfilename, paramDict, prep_uuidstr):
    fileIN = open("bestGWAS" +prep_uuidstr + ".txt")
    import gzip, os, subprocess

    fileOUT = open(outfilename, "w")

    for i in fileIN:
        #i = i[0]
        
        i=i.rstrip().split("\t")
        if len(i)<6: continue
        #in chr22 you can just use this and fileOUT.write str(varG) 
        varG = str(i[3])
        
        #added this to get the input 
        chrnum= str(i[0])
        chrstart= str(i[1])
        chrend= str(i[2])
        vcfFilePath= paramDict["gtexDir"] + "GTEx_Analysis_v7_WGS.vcf.gz"
        print('\n this is vcf path', vcfFilePath)
        command = ["tabix", vcfFilePath ,  str(chrnum) + ":" + str(chrstart) + "-" + str(chrend)]
        print('command is ', command)
        proc = subprocess.Popen(command, stdout = subprocess.PIPE)
        output, err = proc.communicate()
        if len(output)==0:
            print("Error: There is a problem with the tabix output. Check to make sure the vcf file is tabix-indexed. May need to reindex. Exiting. ")
        #sys.exit()
        #variant = output.decode().split("\t")[2]
        genotypesbyn = output.decode().split("\n")
        #added this below, all it does is make sure it picks an SNP instead of indel
        countg=0
        genotypes=[]
        for i in genotypesbyn:
            print(i)
            genotype_test = genotypesbyn[countg].split("\t")
            countg=countg+1
            if len(genotype_test) > 1:
                if len(genotype_test[3]) + len(genotype_test[4]) == 2:
                    genotype_forlist= genotype_test
                    #genotypes= genotypes + genotype_forlist #check gtex func why it is commented out
                    genotypes= genotype_forlist 
        variant = genotypes[2]
        print(variant)
        fileOUT.write(str(variant) + "\t" + "Cancer")
        #for chr22 subset, use below
        #fileOUT.write(str(varG) + "_b37" + "\t" + "Cancer"  )

    fileOUT.close()
    fileIN.close()
#prep_real_GWAS("bestGWAS_var.txt")


# In[7]:

#####This step outputs the selected cold spot
def output_selected_coldspots(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr):
    #fileIN = open ('selectedcoldspots.txt')
    coldspots = find_variant_coldspot('eQTLcoldspot'+prep_uuidstr+'.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    fileOUT = open(outfilename, "w")
    selectedcoldspots2 = list()
    
    for i in coldspots:
        
        #i = i.rstrip().split('\t')
        Chr = str(i[0])
        end = float(i[1])
        bound = float(i[2])
        selectedcoldspots2.append([ (Chr),  (end),  (bound)])              
    

        fileOUT.write("\t".join([ str(Chr),  str(end),  str(bound)] ) + "\n")
              
    
    #fileIN.close()
    fileOUT.close()
    return selectedcoldspots2
#output_selected_coldspots('preVCF.txt')



# In[29]:

#####This step has many readlines, it outputs all variants in the cold spot

def filter_vcf(outfilename, paramDict, geneName, chrNum, pos, tissue, snp, pvalue, prep_uuidstr): 
    #fileIN = open(paramDict["gtexDir"] +'/chr22_subset_gtex_copy2.vcf')
    #fileIN = open(paramDict["gtexDir"] +'qtltools_test/chr22_subset_gtex_copy2.vcf')
    fileIN = open(paramDict["gtexDir"] + 'chr' + chrNum + '_GTEx_Analysis_v7_WGS.vcf')
   
    #selectedcoldspots2 = output_selected_coldspots('preVCF.txt') #already commented out
    selectedcoldspots2 = find_variant_coldspot('eQTLcoldspot' + prep_uuidstr + '.txt', paramDict, geneName, chrNum, pos, tissue, snp, pvalue)
    fileOUT = open(outfilename, "w")
    print("this is selectedcoldspots2")
    print(len(selectedcoldspots2))
    print(selectedcoldspots2)
    
    
    for i in fileIN:
        if i[0]=="#": continue
        i = i.rstrip().split('\t')
        if i[1]=="POS":
            continue      
        chrnum = str(i[0])
        #print(i[1])
        position = float(i[1])


        for  Chr, end, bound in selectedcoldspots2:
            if (Chr)  == chrnum:
                #print('chr is ', Chr, 'end is ', end, 'bound is ', bound)
                #print('position is ', position)
                if position > float(end)  and position < float(bound):  

                    fileOUT.write("\t".join([str(x) for x in i] ) + "\n")                    
                else: 
                    continue
      
    fileIN.close()
    fileOUT.close()
#filter_vcf('filteredgenotype.txt')
