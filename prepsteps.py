import os
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import load_params as LP
#param_file = "param_files/ryo_local.csv"
param_file = "param_files/binkley_local.csv" 
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

        hotspots.append([Chr, start, end, middle, bound])
        #hotspots.append(start)
        #hotspots.append(end)
        #hotspots.append(middle)
        #fileOUT.write(Chr + "\t" + str(start) + "\t" + str(end) + "\t" + str(middle)  +  "\t" + str(bound) + "\n")
    
    fileIN.close()
    #fileOUT.close()
    return hotspots

#####OUtputs cold spot for given eqtl
def find_variant_coldspot(outfilename):
    fileIN = open( "/users/michaelbinkley/desktop/RTCstuffs/besteQTLs.txt" )
    
    #fileOUT = open(outfilename, "w")
    hotspots = ID_hotspots('coldspots2.txt')
    eQTLhotspot = list()


    for i in fileIN:
        
        i = i.rstrip().split('\t')
        gene = str(i[0])
        chrnum = str(i[1])
        snp = (i[3])
        position = float(i[2])
        pvalue = float(i[4])
        
        for  Chr, start, end, middle, bound in hotspots:
            if Chr  != chrnum :
                continue
            else: 
                if position < float(bound) and position > float(end): 
                    #distance = abs(position - middle)
                    eQTLhotspot.append([Chr, end, bound])
                    #fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , gene, snp, str(position), str(pvalue)] ) + "\n")
                else: 
                    continue
    return eQTLhotspot
    
    fileIN.close()
    #fileOUT.close()

#find_variant_coldspot('eQTLcoldspot.txt')

#### This step can probably be truncated 
def output_selected_coldspots(outfilename):
    fileIN = open ('eQTLcoldspot.txt')
    
    fileOUT = open(outfilename, "w")
    selectedcoldspots = list()
    
    for i in fileIN:
        
        i = i.rstrip().split('\t')
        Chr = str(i[0])
        end = str(i[2])
        bound = str(i[4])
        
    

        fileOUT.write("\t".join([ str(Chr),  str(end),  str(bound)] ) + "\n")
        selectedcoldspots.append([ str(Chr),  str(end),  str(bound)])                    
    
    fileIN.close()
    fileOUT.close()
    return selectedcoldspots
#output_selected_coldspots('selectedcoldspots.txt')




# In[8]:

###THis outputs a file with GWAS variants in our cold spots of interest that colocalize. it will output the best GWAS
## first step takes 2 min to run
def find_variant_coldspot_gwas(outfilename):
    fileIN = open(paramDict["gtexDir"] +'/oncoarray_bcac_public_release_oct17.txt')
    fileIN.readline()
    fileOUT = open(outfilename, "w")

    cspots = find_variant_coldspot('eQTLcoldspot.txt')
    #cspots = output_selected_coldspots('selectedcoldspots.txt')
    GWAScoldspots = list()
    for i in fileIN:
        
        i = i.rstrip().split('\t')
        gene = str(i[0])
        chrnum = str(i[2])
        #snp = (i[3])
        position = float(i[3])
        pvalue = float(i[9])
        if pvalue<10e-9:
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
    return GWAScoldspots
    fileIN.close()
    #fileOUT.close()
#find_variant_coldspot_gwas("gwasColdSpot.txt")

def sort_combined_gwas(outfilename):
    #fileIN = open("gwasColdSpot.txt")
    GWAScoldspots = find_variant_coldspot_gwas("gwasColdSpot.txt")
    sortedData = list()
    fileOUT = open(outfilename, "w")
    #for i in fileIN:
    sortedlist = list()
    for i in GWAScoldspots:
        #i=i.rstrip().split("\t")
        i =[i[0], (i[1]), (i[2]), (i[3]), (i[4]), float(i[5])]
        sortedData.append( [(i[0]), (i[1]), (i[2]), i[3], (i[4]), (i[5])])
    sortedScoreData = sorted(sortedData, key = lambda x : (x[5]), reverse=False) #Sort by the element in the column index 
    
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
def output_header(outfilename):
    print(outfilename)
    fileIN = open('sortedGWAS.txt')
    fileOUT = open(outfilename, "w")

    fileOUT.write(fileIN.readline() + "\n")
    #return list
    fileOUT.close()
    fileIN.close()
#output_header("bestGWAS.txt")


def prep_real_GWAS(outfilename):
    fileIN = open("bestGWAS.txt")

    fileOUT = open(outfilename, "w")

    for i in fileIN:
        #i = i[0]
        
        i=i.rstrip().split("\t")
        if len(i)<6: continue
        print(i)
        varG = str(i[3]) 

        fileOUT.write(str(varG) + "_b37" + "\t" + "Cancer"  )

    fileOUT.close()
    fileIN.close()
#prep_real_GWAS("bestGWAS_var.txt")


# In[7]:

#####This step outputs the selected cold spot
def output_selected_coldspots(outfilename):
    #fileIN = open ('selectedcoldspots.txt')
    coldspots = find_variant_coldspot('eQTLcoldspot.txt')
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

def filter_vcf(outfilename): 
    fileIN = open(paramDict["gtexDir"] +'/chr22_subset_gtex_copy2.vcf')
    
    #selectedcoldspots2 = output_selected_coldspots('preVCF.txt')
    selectedcoldspots2 = find_variant_coldspot('eQTLcoldspot.txt')
    fileOUT = open(outfilename, "w")
    
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
                if position > float(end)  and position < float(bound):  

                    fileOUT.write("\t".join([str(x) for x in i] ) + "\n")                    
                else: 
                    continue
      
    fileIN.close()
    fileOUT.close()
#filter_vcf('filteredgenotype.txt')
