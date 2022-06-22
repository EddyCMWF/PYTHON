import numpy as np


def Q10_func(data, Q10, rate, scale, offset):
    # Returns the Q10 function of data for input parameters
    return rate* ( Q10 ** ( scale*(data-offset) ))


# Round to given significant figures
def round2SignifFigs(vals,n):
    mags = 10.0**np.floor(np.log10(np.abs(vals)))  # order of mag's
    outvals = np.around(vals/mags,n-1)*mags             # round(val/omag)*omag
    try:
        outvals[np.where(np.isnan(vals))] = 0.0           # where order of mag = 0, set to zero
    except:
        if np.isnan(outvals):
            outvals=0.0
    #
    return outvals
#

