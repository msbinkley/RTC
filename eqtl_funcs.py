#!/usr/bin/env python3

import numpy as np
import gtex_funcs as GF
import ols_funcs as OF
import scipy.stats
import load_params as LP


#param_file = "param_files/binkley_local.csv" 
param_file = "param_files/trial_local.csv" 
paramDict = LP.load_param_file(param_file)


def get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict):
    '''

    '''
    #expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexExpDir"])

    print('this is gene name 1', geneName)
    #genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
    expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexDir"])
    #added indextodrop in both lines below
    genoVector, indivGenoVector, indextodrop,namestodelete = GF.get_dos_vector(str(chrNum), str(pos), paramDict["gtexDir"])

    #genoVector, __ = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector, indextodrop) 
    genoVector, indivExpVector, indextodropforexp, indextokeepforexp = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector, indextodrop, namestodelete) #changed above to this 
    print('indextodrop is', indextodrop) 
    print('indextodrop for exp is', indextodropforexp)
    print('indextokeep for exp is', indextokeepforexp)
    expVector=[i for n, i in enumerate(expVector) if n not in indextodropforexp]
    print('ef genovector: ', genoVector)
    print('ef expvector: ', expVector)
    results = scipy.stats.linregress(genoVector,expVector)

    residuals = OF.calculate_residuals(genoVector, expVector, results.slope, results.intercept)

    return results.slope, results.intercept, residuals, indivExpVector, indextodropforexp 

def get_eqtl_stats_real(chrNum, pos, geneName, tissue, paramDict):
    '''

    '''
    #expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexExpDir"])

    print('this is gene name 2', geneName)
    #genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
    expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexDir"])
    #added index to drop in both lines below 
    genoVector, indivGenoVector, indextodrop, namestodelete  = GF.get_dos_vector(str(chrNum), str(pos), paramDict["gtexDir"])
    #genoVector, __ = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector, indextodrop)
    genoVector, indivExpVector, indextodropforexp, indextokeepforexp = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector, indextodrop, namestodelete) #changed above to this
    expVector=[i for n, i in enumerate(expVector) if n not in indextodropforexp]
    #results = scipy.stats.linregress(genoVector,expVector)


    #residuals = OF.calculate_residuals(genoVector, expVector, results.slope, results.intercept)
    #results.slope, results.intercept,
    return  expVector, indivExpVector