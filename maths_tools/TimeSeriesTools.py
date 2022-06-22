#!/usr/bin/python
# ECP - various Time Series tools

import numpy as np
import scipy.stats as stats

# get overlap indexes
def overlap_indexes(Atime,Btime,verbose=False):
    start = max(Atime[0],Btime[0])
    end   = min(Atime[-1],Btime[-1])
    if verbose:
        print(start)
        print(end)

    if start <= end:
        overlapA = (Atime>=start) & (Atime<=end) 
        overlapB = (Btime>=start) & (Btime<=end)
    else:
        overlapA = None
        overlapB = None
        
    return overlapA, overlapB 
    
# Calcualte mena bias of two time-series with non common time vectors 
def meanbias(Adata,Bdata,Atime,Btime,verbose=False):
    #fetch overlap indexes
    overlapA,overlapB=overlap_indexes(Atime,Btime)
    A=Adata[overlapA]
    B=Bdata[overlapB]
    if verbose:
        print(A.shape, B.shape)
        print('JULES: ',Atime[overlapA][0:10])  #,Atime[overlapA][-1])
        print('Site: ',Btime[overlapB][0:10])  #,Btime[overlapB][-1])
    return np.mean(A-B)
    
# Calculate standard deviation of two time series with non common time vectors
def stddev(Adata,Bdata,Atime,Btime,verbose=False):
    #fetch overlap indexes
    overlapA,overlapB=overlap_indexes(Atime,Btime)
    A=Adata[overlapA]
    B=Bdata[overlapB]
    if verbose:
        print( A.shape, B.shape)
    return np.std(A-B)

# Calculate Pearson's R correlation of two time series with non common time vectors
def pearsonr(Adata,Bdata,Atime,Btime,verbose=False):
    #fetch overlap indexes
    overlapA,overlapB=overlap_indexes(Atime,Btime)
    A=Adata[overlapA]
    B=Bdata[overlapB]
    if verbose:
        print( A.shape, B.shape)
    return stats.pearsonr(A,B)









