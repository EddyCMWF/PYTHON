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
run = sys.argv[1]
start_year = sys.argv[2]
end_year = sys.argv[3]

##############################################################################
# 1. Define directories and files etc.
##################################
#
out_dir = '/prj/GREENHOUSE/SC_simulator/OPERATIONAL_output/work_'+run+'/' 

spin1_outfile = out_dir+'J4.3_'+run+'_spin1.monthly_mean.nc'
print('spin1_outfile = '+spin1_outfile)

spin1_dumpfile= out_dir+'J4.3_'+run+'_spin1.dump.'+end_year+'0101.0.nc'
print('spin1_dumpfile = '+spin1_dumpfile)

spin2_startdump_file = out_dir+'J4.3_'+run+'_spin2_startdump.nc'
print('spin2_startdump_file = '+spin2_startdump_file)

##############################################################################
# 2. Define parameters:
##################################
#
print('Defining parameters')
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


##############################################################################
# 4. Calculate equi CS
#     Using the Eleanor/Alberto/Eddy method calculate the equilibrium
#       soil carbon based on the soil resp and litter functions
#     
#     Output our first estimate of the equilibrium CS:
#       CS_2 = [CS_dpm_2,CS_rpm_2,CS_bio_2,CS_hum_2]
#
#     RPM and DPM should have reasonable equilibrium estimates now
#
##################################
#os.system('pwd')

# 4.1 Read in data from Spin 1 run
# Remove first year to avoid funny initial litterfall values
# squeeze and transpose soil/pft dim to position 0
inf      = nc.Dataset(spin1_outfile,'r')
loclit_c = inf.variables['lit_c'][12:,:].squeeze().transpose(1,0,2)
c_veg    = inf.variables['c_veg'][12:,:].squeeze().transpose(1,0,2)
resp_s   = inf.variables['resp_s'][12:,:].squeeze().transpose(1,0,2)
cs       = inf.variables['cs'][12:,:].squeeze().transpose(1,0,2)
frac     = inf.variables['frac'][12:,:].squeeze().transpose(1,0,2)
inf.close()

# Extract PFT frac and then calc Veg frac
PFTfrac = frac[:nPFTs,:]
VEGfrac = np.sum(PFTfrac,axis=0)
if 'CHESS' in run:
    ICEmask = np.where(VEGfrac[0,:]<-1)[0]
else:
    ICEmask = np.where(frac[8,0,:]>0.1)[0]

#print ICEmask

# 4.2 Calculate local litterfall
print 'Litterfall Calculations'
#Lit_C = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, loclit_c, c_veg, \
#                                           Per_PFT=True )
Lit_C=loclit_c*PFTfrac
#
# 4.3 Split the total litter fall into DPM and RPM components
Lit_SCpools = J_SC_F.JULES_LITTER_to_SCpool(Lit_C)

# 4.4 Soil Respiration Factor component
Soil_Resp_Fact = (resp_s)/(cs*(1-beta_r))

# 4.5 Calculate the soil carbon for each pool using the albmar/ECP equations
# DPM
C_dpm  = (np.sum(Lit_SCpools[0,:,:], axis=0)  /  \
              np.sum(Soil_Resp_Fact[0,:,:], axis=0) ) /  \
              (3600.*24.*360.)
# RPM
C_rpm  = (np.sum(Lit_SCpools[1,:,:], axis=0)       /  \
              np.sum(Soil_Resp_Fact[1,:,:], axis=0) ) /  \
              (3600.*24.*360.)

    # BIO and HUM factors:
dpm_rpm_term   = (kappa_s_RothC[0]*C_dpm)+(kappa_s_RothC[1]*C_rpm)
bio_factor =  0.46 / ( ( (1/beta_r)-1 )*kappa_s_RothC[2] )
hum_factor =  0.54 / (( (1/beta_r)-1 )*kappa_s_RothC[3] )
# BIO
C_bio   = bio_factor*dpm_rpm_term
# HUM
C_hum   = hum_factor*dpm_rpm_term
#
# Mask out ICE points
cs_2 = np.array( [C_dpm,C_rpm,C_bio,C_hum] )
cs_2[:,ICEmask]=0.0

# set values less than or equal to zero to 1e-6
cs_2[cs_2<=0]=1e-6

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


