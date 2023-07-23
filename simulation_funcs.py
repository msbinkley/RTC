#!/usr/bin/env python3

import pandas

def create_coldspot_file(hotspotFilePath, coldspotFilePath):
    '''
    Input:
	hotspotfilePath
    Output:
        coldspotFilePath

    IN PROGRESS
    '''
    fileOUT = open(outfilename, "w")
    hotspots = open('hotspots.txt')
    for j in hotspots:
        j = j.rstrip().split('\t')
        Chr = str(j[0])
        start = float(j[1])
        end = float(j[2])

    egene.close()
    fileOUT.close()
    hotspots.close()


def get_colocalized_coldspot(snpList1, snpList2, coldspotFilePath):
    '''
    Input:
        snpList1: List of GWAS SNPs.
        snpList2: List of eQTL SNPs.
    Output:
        coldspotList: Returns a list of coldspots with a colocalized GWAS and eQTL.

    IN PROGRESS
    '''
    pass


