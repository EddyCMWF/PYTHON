#!/usr/bin/python2.7
#
############################################################
# 
# Program: JULES_with_SCequi.py
# Author: E. Comyn-Platt, edwcom@ceh.ac.uk
# Purpose: To perform a multi-stage spin up of JULES
#          incorporating SC simulator method
###########################################################

import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
import JULES_SC_functions as J_SC_F
import netCDF4 as nc


##############################################################################
# 0. Read in parsed arguments
##################################
#
infile = sys.argv[1]
outfile = sys.argv[2]

##############################################################################
# 1. Define directories and files etc.
##################################
#

##############################################################################
# 2. Define parameters:
##################################
#
print 'Defining parameters'
# 2.0 Define Variables:
# 2.0.1 Soil Specific Respiration Rate for RothC pools and single soil pool
kappa_s_RothC  = np.array([ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10 ])
kappa_s_1SC    = [0.5e-8]
# 2.0.2 Soil_Resp_Fractor beta_r in documetation, formulation from JULES src
clay_content=0.23   # as in JULES
beta_r  = 1.0 / (4.0895 + 2.672*  np.exp(-0.0786 * 100.0 * clay_content))
#
dz_soil = 0.1   # top soil layer is 10cm  
#
nPFTs=5
nPOOLs=4

################################
# Set initial guess parameters:
alpha=-1



##############################################################################
# 4. Calculate equi CS by fitting Cs=(C0-C_equi)*e(-a*kappa*t)+C_equi,
#     using the Eddy fit of exponential decay method     
#
##################################
#os.system('pwd')

# 4.1 Read in data from Spin 1 run
# squeeze and transpose soil/pft dim to position 0
inf      = nc.Dataset(infile,'r')
time_sec = inf.variables['time'].squeeze()
# reset time to start of run
time_sec -= time_sec[0] 
cs_in    = inf.variables['cs'][:].squeeze().transpose(1,0,2)
frac     = inf.variables['frac'][:].squeeze().transpose(1,0,2)
inf.close()

# Extract PFT frac and then calc Veg frac
if 'CHESS' in infile:
    ICEmask = np.where(frac[0,:]<-1)[0]
else:
    ICEmask = np.where(frac[8,0,:]>0.1)[0]


model_cs = np.zeros_like(cs_in)
for iPOOL in nPOOLs:
    for pt in ICEmask



# write new cs to spin2_startdump_file
print "Writing SC equilibrium output to: ", spin2_startdump_file
inf = nc.Dataset(spin1_dumpfile,'r')
outf= nc.Dataset(spin2_startdump_file,'w')
cs_1=inf.variables['cs'][:].squeeze()
for dim in inf.dimensions:
    outf.createDimension( str(dim),len(inf.dimensions[dim]) )
        
for var in inf.variables:
    outvar = outf.createVariable( str(var), \
                                  inf.variables[var].dtype, \
                                  inf.variables[var].dimensions )
    if str(var)=='cs':
        outvar[:]=cs_2
    else:
        outvar[:]=inf.variables[str(var)][:]

inf.close()
outf.close()


