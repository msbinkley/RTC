def output_coldspot(outfilename, paramDict):
    '''
    This step may not be necessary. But it filters out any eQTLs that are outside of coldspots
    '''

    #egene = open(paramDict["gtexDir"] + '/QTLinterest.txt')
    #egene = open(paramDict["gtexDir"] + '/Breast_Mammary_Tissue.v7.egenes.txt')
    #ajay changed to this
    egene = open(paramDict["gtexDir"] + 'Breast_Mammary_Tissue.v7.egenes.txt')
    egene.readline()

    fileOUT = open(outfilename, "w")
    #hotspots = open(paramDict["gtexDir"] + 'coldspots.txt')
    #ajay changed to this
    hotspots = open(paramDict["gtexDir"] + 'coldspots2.txt')
    for j in hotspots: 
        j = j.rstrip().split('\t')
        Chr = str(j[0])
        start = int(j[1])
        end = int(j[2])
        middle = str(j[3])
        bound = str (j[4])        
    
    for i in egene:
        
        i = i.rstrip().split('\t')
        chrnum = str(i[2])
        gene = i[0]
        position = float(i[14])
        if chrnum == Chr and position > start and position <end :
            continue
        else: 
            fileOUT.write("\t".join([str(x) for x in i])+ "\n")
            
    
    egene.close()
    fileOUT.close()
    hotspots.close()


def extract_best_eqtl(paramDict):
    '''
    Need a function to then identify the best eqtl for a given gene. Again, this step may not be necessary
    '''
    #fileIN = open(paramDict["gtexDir"] + '/QTLinterest.txt')    
    fileIN = open(paramDict["gtexDir"] + '/Breast_Mammary_Tissue.v7.egenes.txt')
    fileIN.readline()
    longestGeneCoordinates = list()
    currentLine = fileIN.readline().rstrip().split('\t') 
    currentGeneCoordinates = [str(currentLine[0]), str(currentLine[2]),float(currentLine[14]), str(currentLine[18]), float(currentLine[23])]#This is a variable that holds the information for the "current" gene.
    #"Current" gene is the gene that we are currently working to figuring out the longest coordinates.
    for counter, i in enumerate(fileIN): #Go through each line
        i = i.rstrip().split("\t")
        i =[str(i[0]), str(i[2]), float(i[14]), str(i[18]), float(i[23])]

          
        if i[0] == currentGeneCoordinates[0]: #If this line has the same gene as "current" SNP, then continue to compare the coordiantes.
            currentDistance = abs(currentGeneCoordinates[23])
            thisLinesDistance = abs(i[23])
            if thisLinesDistance < currentDistance: #If this line is shorter than our "current" gene, then update the "current" gene.
                currentGeneCoordinates = i
                
        if i[0] != currentGeneCoordinates[0]: 
            #If this line's gene does not match the currentGene's gene, then we must have moved on to the next gene.
            #Since we're now at a new gene, let's store the currentGene
            longestGeneCoordinates.append(currentGeneCoordinates )
            #Also, let's update the currentGene to match this gene. 
            currentGeneCoordinates = i
        if counter %100000==0:
            print(counter)

       
    fileOUT = open(paramDict["tmpSmall"] + "/besteQTLs.txt", "w")
    for i in longestGeneCoordinates:
        print(i)
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")  #Join will take a list, and create ONE string. where each element of the list is separated by the character in the quotation marks (here it is a "\t")
        #Also,we're converting the items into a string - because we can't output integers. 
    fileOUT.close() 





    
def extract_best_eqtl_real(paramDict):
    '''
    Need a function to then identify the best eqtl for a given gene. Again, this step may not be necessary
    '''
    #fileIN = open(paramDict["gtexDir"] + '/QTLinterest.txt')   
    #fileIN = open(paramDict["gtexDir"] + '/Breast_Mammary_Tissue.v7.egenes.txt')
    #fileIN = open(paramDict["gtexDir"] + '/Artery_Coronary.v7.egenes.txt')
    #fileIN = open(paramDict["gtexDir"] + '/egene_and_expbed_files/GTEx_Analysis_v7_eQTL/Colon_Sigmoid.v7.egenes.txt')
    fileIN = open(paramDict["gtexDir"] + '/egene_and_expbed_files/GTEx_Analysis_v7_eQTL/Prostate.v7.egenes.txt')
    fileIN.readline()
    longestGeneCoordinates = list()
    currentLine = fileIN.readline().rstrip().split('\t') 
    currentGeneCoordinates = [str(currentLine[0]), str(currentLine[2]), float(currentLine[14]), str(currentLine[11]),  str(currentLine[18]), float(currentLine[23])]#This is a variable that holds the information for the "current" gene.
    #"Current" gene is the gene that we are currently working to figuring out the longest coordinates.
    for counter, i in enumerate(fileIN): #Go through each line
        i = i.rstrip().split("\t")
        if str(i[11])== ' ' : continue
        i =[str(i[0]), str(i[2]), float(i[14]), str(i[11]),  str(i[18]), float(i[23])]

          
        if i[0] == currentGeneCoordinates[0]: #If this line has the same gene as "current" SNP, then continue to compare the coordiantes.
            currentDistance = abs(currentGeneCoordinates[23])
            thisLinesDistance = abs(i[23])
            if thisLinesDistance < currentDistance: #If this line is shorter than our "current" gene, then update the "current" gene.
                currentGeneCoordinates = i
                
        if i[0] != currentGeneCoordinates[0]: 
            #If this line's gene does not match the currentGene's gene, then we must have moved on to the next gene.
            #Since we're now at a new gene, let's store the currentGene
            longestGeneCoordinates.append(currentGeneCoordinates )
            #Also, let's update the currentGene to match this gene. 
            currentGeneCoordinates = i
        if counter %100000==0:
            print(counter)

       
    fileOUT = open(paramDict["tmpSmall"] + "/besteQTLs_real_Prostate.txt", "w")
    for i in longestGeneCoordinates:
        print(i)
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")  #Join will take a list, and create ONE string. where each element of the list is separated by the character in the quotation marks (here it is a "\t")
        #Also,we're converting the items into a string - because we can't output integers. 
    fileOUT.close() 



