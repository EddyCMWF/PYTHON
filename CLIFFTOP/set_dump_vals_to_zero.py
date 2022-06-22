#!/bin/env python2.7


#import numpy as np
import netCDF4 as nc
import sys

infile=sys.argv[1]

zero_var_list=['frac_agr_prev','frac_past_prev','wood_prod_fast','wood_prod_slow','wood_prod_med']

inf=nc.Dataset(infile,'a')

for var in zero_var_list:
    inf.variables[var][:]=0.


inf.close()





