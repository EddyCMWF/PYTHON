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
work_dir= '/home/users/ecomynplatt/SC_simulator/work_'+run+'/'
out_dir = '/work/scratch/ecomynplatt/SC_simulator/work_'+run+'/' 
os.chdir(work_dir)

spin2_dumpfile= out_dir+'J4.3_'+run+'_spin2.dump.'+end_year+'0101.0.nc'
print 'spin2_dumpfile = '+spin2_dumpfile
spin2_startdump_file = out_dir+'J4.3_'+run+'_spin2_startdump.nc'
print 'spin2_startdump_file = '+spin2_startdump_file

spin3_startdump_file = out_dir+'J4.3_'+run+'_spin3_startdump.nc'
print 'spin3_startdump_file = '+spin3_startdump_file

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


##############################################################################
# 5. 2nd CS modification
#      - We now shold have DPM/RPM and BIO in equilibrium/steady state
#      - the ratio of the equi CS_bio calculation to the the newly spun up
#           CS_bio is used to modify the CS_hum and CS_rpm
#      -  CS_dpm_4 = CS_dpm_3
#      -  CS_rpm_4 = CS_rpm_2 
#      -  CS_bio_4 = CD_bio_3
#      -  CS_hum_4 = CS_hum_2 * (CS_bio_3/CS_bio_2)
#      -  CS_4 = [CS_dpm_4,CS_rpm_4,CS_bio_4,CS_hum_4]
##################################

# open spin2 startdump file and read in soil carbon
cs_2 = nc.Dataset(spin2_startdump_file,'r').variables['cs'][:]

# open spin2 [end] dumpfile
#print spin2_dumpfile
inf=nc.Dataset(spin2_dumpfile,'r')
# read in cs_3 from spin2
cs_3 = inf.variables['cs'][:]
# calculate cs_4 as detailed above
cs_4=cs_3.copy()
cs_4[3,:]= cs_2[3,:] * (cs_3[2,:]/cs_2[2,:])
cs_4[1,:]= cs_2[1,:] ###* (cs_3[2,:]/cs_2[2,:])

print "Writing Hummus Adjust output to: "+spin3_startdump_file
outf=nc.Dataset(spin3_startdump_file,'w')

for dim in inf.dimensions:
    outf.createDimension( str(dim),len(inf.dimensions[dim]) )
    
for var in inf.variables:
    outvar = outf.createVariable( str(var), \
                                  inf.variables[var].dtype, \
                                  inf.variables[var].dimensions )
    if str(var)=='cs':
        outvar[:]=cs_4
    else:
        outvar[:]=inf.variables[str(var)][:]

inf.close()
outf.close()


