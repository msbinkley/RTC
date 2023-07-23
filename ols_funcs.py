#!/usr/bin/env python3

import numpy
import scipy.stats

def calculate_residuals(xVec,yVec, beta, intercept):
    '''
    Calculates residuals using ordinary linear regression
    '''
    
    print('\t calculate residuals')
    xVec = numpy.array(xVec)
    residuals = yVec - (xVec*beta + intercept)
    return residuals


