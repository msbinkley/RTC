import deprecated_funcs as DF
import eqtl_funcs as EF

def prep_real_GWAS(outfilename, paramDict):
    fileIN = open("bestGWAS.txt")
    import gzip, os, subprocess
    fileOUT = open(outfilename, "w")
    #fileOUTint = open("bestGWAS_var_int3.txt", "a+")
    for i in fileIN:
        #i = i[0]
        
        i=i.rstrip().split("\t")
        if len(i)<6: continue
        print(i)
        #in chr22 you can just use this and fileOUT.write str(varG) 
        varG = str(i[3])
        fileOUT.write(str(varG) + "_b37" + "\t" + "Cancer"  ) #added this from old
        #fileOUTint.write(str(varG) + "\t" + "Cancer" + "\n") #this is for null set for enrichment
        #added this to get the input 
        #chrnum= str(i[0])
        #chrstart= str(i[1])
        #chrend= str(i[2])
        #vcfFilePath= paramDict["gtexDir"] + "chr"+ str(chrnum)+ "_GTEx_Analysis_v7_WGS.vcf.gz"
        #print('\n this is vcf path', vcfFilePath)
        #command = ["tabix", vcfFilePath ,  str(chrnum) + ":" + str(chrstart) + "-" + str(chrend)]
        #print('command is ', command)
        #proc = subprocess.Popen(command, stdout = subprocess.PIPE)
        #output, err = proc.communicate()
        #if len(output)==0:
        #    print("Error: There is a problem with the tabix output. Check to make sure the vcf file is tabix-indexed. May need to reindex. Exiting. ")
        #genotypesbyn = output.decode().split("\n")
        ##added this below, all it does is make sure it picks an SNP instead of indel
        #countg=0
        #genotypes=[]
        #for i in genotypesbyn:
        #    print(i)
        #    genotype_test = genotypesbyn[countg].split("\t")
        #    countg=countg+1
        #    if len(genotype_test) > 1:
        #        if len(genotype_test[3]) + len(genotype_test[4]) == 2:
        #            genotype_forlist= genotype_test
        #            #genotypes= genotypes + genotype_forlist #check gtex funcs for why it is commented out
        #            genotypes= genotype_forlist 
        #variant = genotypes[2]
        #print(variant)
        #fileOUT.write(str(variant) + "\t" + "Cancer"  )
        
        #for chr 22 subset, use below 
        #fileOUT.write(str(varG) + "_b37" + "\t" + "Cancer"  )

    fileOUT.close()
    #fileOUTint.close()
    fileIN.close()
#prep_real_GWAS("bestGWAS_var.txt")

def prep_real_eQTL(paramDict, outfilename, snp): #this is the one original one that i also edited
    #DF.extract_best_eqtl_real(paramDict)
    #fileIN = open("besteQTLs.txt") #instead of reading the input file, we can just give it the snp from main script
    fileOUT = open(outfilename, "w")
    #for i in fileIN:
    #    i=i.rstrip().split("\t")
    #    varG = str(i[3])
    #    fileOUT.write(str(varG) + "\n")
    fileOUT.write(str(snp) + "\n")
    fileOUT.close()
    #fileIN.close()
    
def prep_real_eQTL_old(paramDict, outfilename): #this is the one I edited, dont use anymore
    DF.extract_best_eqtl_real(paramDict)
    fileIN = open("besteQTLs.txt")
    fileOUT = open(outfilename, "w")
    for i in fileIN:
        i=i.rstrip().split("\t")
        varG_tomatch = str(i[0])
        
        fileOUT.write(str(varG) + "\n")

    fileOUT.close()
    fileIN.close()
    
    
def prepare_perm_file_real(paramDict, outfilename):  
    #fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed')
    #fileIN = open('/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/GTExExpfile.bed')
    fileIN = open(paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed")
    causalvar = open("besteQTL_var.txt")
    fileIN.readline()
    fileOUT = open(outfilename, "w")
    for j in causalvar: 
        j = j.rstrip().split('\t')
        cvar = j[0]
    for i in fileIN: 
        i = i.rstrip().split('\t')
        gene = i[3]
        chrn = i[0]
        start = i[1]
        end = i[2]
        dist = i[4]
        sense = i[5]
        
        fileOUT.write("\t".join([ str(gene), str(chrn),  str(start) , str(end), sense, dist, dist, str(cvar), chrn , 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA' ] ) + "\n")

    fileIN.close()
    fileOUT.close()
    
### RUn this for the real eQTL, this does not use residuals or calculate residuals.
import load_params as LP
import eqtl_funcs as EF


def get_bed_real(chrNum, pos, geneName, tissue, paramDict):
    import deprecated_funcs as DF
    import eqtl_funcs as EF  
    expV2, indivVector2 = EF.get_eqtl_stats_real(chrNum, pos, geneName, tissue, paramDict)  

#print(slope)
#print(intercept)
#print(residuals)
    print("expv2 = ", expV2)
    print("indivVector2 = ",indivVector2)

    def output_gene_info_real(geneName, paramDict, chrNum, pos, tissue):
    
        fileIN = open(paramDict["tmpDir"] + '/chrname.txt')
        selected = list()
        for i in fileIN: 
            i = i.rstrip().split('\t')
        
            if i[3]== geneName:
                selected = i
                break
            
        return selected

  
    def output_header2_real(indivVector2):    
            
        header = "\t".join(['#chr',  'start', 'end' , 'gene' , 'length' , 'strand'] + [str(x) for x in indivVector2] )
        return header

        fileIN.close()

    def output_bed_real(indivVector2, expV2, filepath, geneName, paramDict, chrNum, pos, tissue):
        selected = output_gene_info_real(geneName, paramDict,chrNum, pos, tissue)
        header = output_header2_real(indivVector2)
        with open(filepath, 'w') as outfile:
            outfile.write(header + "\n")
            geneLocation = selected
            outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expV2]) + "\n")
    print(geneName)   
    #i uncommented this line below, need to create this file 
    output_bed_real(indivVector2, expV2, paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed", geneName, paramDict, chrNum, pos, tissue)

#get_bed_real(chrNum, pos, geneName, tissue, paramDict):
    


