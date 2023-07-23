
import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import allstepslist as ASL
import runRTC as rrtc
import numpy as np

#geneName = "ENSG00000172404.4"
#chrNum = "22"
#pos = "40407887"
#tissue = "Breast_Mammary_Tissue"

#def get_result(result):
#    global results
#    results.append(result)
    
def run_perm_H1(paramDict, geneName, chrNum, pos, tissue, numperm, prep_uuidstr):
    numperm = numperm
    #outlist = list()
    import multiprocessing as mp
    import numpy as np
    import time
    print("inside H1")
    #def the_H1_function(paramDict, geneName, chrNum, pos, tissue, numperm):
    outlist = list()
    import uuid
    uuidval = uuid.uuid4()
    print(uuidval)
    uuidvalstr = str('_uuid_' + str(uuidval))

    ASL.run_all_steps_H1(paramDict, chrNum, uuidvalstr, prep_uuidstr)
    rrtc.run_RTC(paramDict, geneName, chrNum, pos, tissue, uuidvalstr, prep_uuidstr)

        #fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/rtc_results2.txt')
        #for sherlock
        #fileIN = open('/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/rtc_results2.txt')
    fileIN = open(paramDict["gtexDir"] + 'rtc_results2' + uuidvalstr + '.txt')
    fileIN.readline()
    genes = list()
    for j in fileIN:
        g = j.rstrip().split(" ")
        name = g[19]

        #genes.append(name)
        outlist.append([geneName, chrNum, pos, tissue, name])
        #first = (j[0]) 
        
    fileIN.close()
    #ajay added
    print('outlist is this ', outlist)
    import os
    import glob
    myfiles1=paramDict["gtexDir"]+"tmpSmall/"+"*"+uuidvalstr+"*"
    myfiles2= paramDict["gtexDir"]+"*"+uuidvalstr+"*"
    myfiles = glob.glob(myfiles1) + glob.glob(myfiles2)
    filelist = myfiles
    print(myfiles)
    print(len(myfiles))
    if len(uuidvalstr) >=10 :
        print("good to go")
        if paramDict["gtexDir"][:-1] not in filelist:
            if paramDict["gtexDir"] not in filelist:
                if paramDict["gtexDir"]+"chr1_GTEx_Analysis_v7_WGS.vcf" not in filelist: #double check to make sure uuid is not null
                    print("we are good to remove uuid files")
                    #Iterate over the list of filepaths & remove each file.
                    for filePath in filelist:
                        print("outside if filePath:", filePath)
                        try:
                            if "_uuid_" in filePath :
                                if os.path.isfile(filePath):
                                    print("inside if filePath:", filePath)
                                    os.remove(filePath)
                        except:
                            print("Error while deleting file : ", filePath)
    return (outlist)
 
 #   return outlist

def run_perm_H0(paramDict, geneName, chrNum, pos, tissue, numperm, prep_uuidstr):
    numperm = numperm
    #outlist2 = list()
    #import multiprocessing as mp
    import numpy as np
    import time
    print("inside H0")
    #def the_H0_function(paramDict, geneName, chrNum, pos, tissue, numperm):
    outlist2 = list()
    import uuid
    uuidval = uuid.uuid4()
    print(uuidval)
    uuidvalstr = str('_uuid_' + str(uuidval))

    ASL.run_all_steps_H0(paramDict, chrNum, uuidvalstr, prep_uuidstr)
    rrtc.run_RTC(paramDict, geneName, chrNum, pos, tissue, uuidvalstr, prep_uuidstr)

        
        #fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/rtc_results2.txt')
        #for sherlock
        #fileIN = open('/oak/stanford/groups/emoding/scripts/RTCscripts/RTCstuffs/tmpSmall/rtc_results2.txt')
    fileIN = open(paramDict["gtexDir"] + 'rtc_results2' + uuidvalstr + '.txt')
    fileIN.readline()
    genes = list()
    for j in fileIN:
        g = j.rstrip().split(" ")
        name = g[19]
        #genes.append(name)
        outlist2.append([geneName, chrNum, pos, tissue, name])
        #first = (j[0]) 
        
    fileIN.close()
    #ajay added
    print('outlist2 is this ', outlist2)
    
    import os
    import glob
    myfiles1=paramDict["gtexDir"]+"tmpSmall/"+"*"+uuidvalstr+"*"
    myfiles2= paramDict["gtexDir"]+"*"+uuidvalstr+"*"
    myfiles = glob.glob(myfiles1) + glob.glob(myfiles2)
    filelist = myfiles
    print(myfiles)
    print(len(myfiles))
    if len(uuidvalstr) >=10 :
        print("good to go")
        if paramDict["gtexDir"][:-1] not in filelist:
            if paramDict["gtexDir"] not in filelist:
                if paramDict["gtexDir"]+"chr1_GTEx_Analysis_v7_WGS.vcf" not in filelist: #double check to make sure uuid is not null
                    print("we are good to remove uuid files")
                    #Iterate over the list of filepaths & remove each file.
                    for filePath in filelist:
                        print("outside if filePath:", filePath)
                        try:
                            if "_uuid_" in filePath :
                                if os.path.isfile(filePath):
                                    print("inside if filePath:", filePath)
                                    os.remove(filePath)
                        except:
                            print("Error while deleting file : ", filePath)
    return outlist2





