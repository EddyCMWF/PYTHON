#!/usr/bin/env python

################################################################################
# 
# Program: SC_simulator.py
# 
# Python Script to simulate soil carbon estimate for a long-term equilibrium
#   condition, i.e. dCs/dt=0 .'. litterfall=microbial_soil_respiraiton       
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
################################################################################
#
# 
#  
################################################################################
#
#
import pylab as plt
import numpy as np
import netCDF4 as nc
import JULES_SC_functions as J_SC_F
#

fill_value=-999.

#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/WFDEI/'
PreSpin_file = 'dumps/pre_run_output.dump.19800101.0.nc'
BigSpin_file = 'JULES_v4.JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'
#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
out_tag = 'WFDEI_simulated_SC_JULES-dpm-rpm'
#

#####################################################################################
# 2. Define or Read in the various data:
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
#
# 2.4 Read in QuickSpin from JULES output
print 'Reading in Soil Carbon from JULES Pre-Spin data'
inf = nc.Dataset(JULES_output_dir+PreSpin_file,'r')
lats_orig        = inf.variables['latitude'][:].squeeze()
lons_orig        = inf.variables['longitude'][:].squeeze()
cs_JULES         = inf.variables['cs'][:].squeeze()
inf.close()
#
# 2.4 Read in BigSpin CS from JULES output
print 'Reading in Soil Carbon from JULES Big-Spin dump'
inf         = nc.Dataset(JULES_output_dir+JULES_file_BS,'r')
cs_JULES_BS = inf.variables['cs'][:].squeeze()
inf.close()
#


# Calculate C_bio and C_hum from QuickSpin C_dpm and C_rpm
C_dpm_mean = np.mean(cs_JULES[-120:,0,:],axis=0)
C_rpm_mean = np.mean(cs_JULES[-120:,1,:],axis=0)


# See Eddy Equations:
# 6.3 Biomass:
#  C_bio =  ( 0.46 / ( (1/beta_r)-1 )*kappa_bio ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
# 6.4 Hummus:
#  C_hum =  ( 0.54 / ( (1/beta_r)-1 )*kappa_hum ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
#      Calculate Cpool based on Q10 temp function and RothC temp function
dpm_rpm_term   = (kappa_s_RothC[0]*C_dpm_mean)+(kappa_s_RothC[1]*C_rpm_mean)
bio_factor =  0.46 / ( ( (1/beta_r)-1 )*kappa_s_RothC[2] )
hum_factor =  0.54 / (( (1/beta_r)-1 )*kappa_s_RothC[3] )
# 6.3:
C_bio_mean   = bio_factor*dpm_rpm_term
# 6.4:
C_hum_mean   = hum_factor*dpm_rpm_term
#
# 6.5 append to single array:
C_4pools_mean   = np.array( [C_dpm_mean,C_rpm_mean,C_bio_mean,C_hum_mean] )
#


cs_JULES_mean=np.mean(cs_JULES[-120:,:,:],axis=0)
cs_JULES_BS_mean=np.mean(cs_JULES_BS[-120:,:,:],axis=0)

n96_landindexfile='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

# Open and read the land index file
inf=nc.Dataset(n96_landindexfile,'r')
index=inf.variables['index_2D'][:]
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
inf.close()

C_4pools_mean_2D= C_4pools_mean[:,index]
cs_JULES_2D=cs_JULES_mean[:,index]
cs_JULES_BS_2D=cs_JULES_BS_mean[:,index]
for i in range(4):
    C_4pools_mean_2D[i,:,:][np.where(index.mask==True)]=fill_value
    cs_JULES_2D[i,:,:][np.where(index.mask==True)]=fill_value
    cs_JULES_BS_2D[i,:,:][np.where(index.mask==True)]=fill_value

C_4pools_mean_2D=np.ma.masked_equal(C_4pools_mean_2D,-999.)
cs_JULES_2D=np.ma.masked_equal(cs_JULES_2D,-999.)
cs_JULES_BS_2D=np.ma.masked_equal(cs_JULES_BS_2D,-999.)

ABS_maxes = [ 0.1, 5, .6, 20. ]
diff_ranges= [ 0.05, 3, 0.3, 10. ]

SC_pools=['DPM','RPM','BIO','HUM']


ncols=5
nrows=4
for i in range(4):
    
    plt.subplot(nrows,ncols,(i*ncols)+1)
    plt.imshow(C_4pools_mean_2D[i,:,:],origin='bottom',vmin=0,vmax=ABS_maxes[i])
    if i==0:
        plt.title('Simulated SC')
    plt.ylabel(SC_pools[i])
    plt.colorbar(fraction=0.05)
    
    plt.subplot(nrows,ncols,(i*ncols)+2)
    plt.imshow(cs_JULES_2D[i,:,:],origin='bottom',vmin=0,vmax=ABS_maxes[i])
    if i==0:
        plt.title('JULES SC - Quick Spin')
    plt.colorbar(fraction=0.05)
        
    
    plt.subplot(nrows,ncols,(i*ncols)+3)
    plt.imshow(cs_JULES_BS_2D[i,:,:],origin='bottom',vmin=0,vmax=ABS_maxes[i])
    if i==0:
        plt.title('JULES SC - Big Spin')
    plt.colorbar(fraction=0.05)
        
    
    plt.subplot(nrows,ncols,(i*ncols)+4)
    plt.imshow(C_4pools_mean_2D[i,:,:]-cs_JULES_BS_2D[i,:,:],origin='bottom', \
               cmap='RdBu',vmin=(0-diff_ranges[i]),vmax=diff_ranges[i])
    if i==0:
        plt.title('Difference (Sim - JULES_BS)')
    plt.colorbar(fraction=0.05)
        
    
    plt.subplot(nrows,ncols,(i*ncols)+5)
    plt.imshow(cs_JULES_2D[i,:,:]-cs_JULES_BS_2D[i,:,:],origin='bottom', \
               cmap='RdBu',vmin=(0-diff_ranges[i]),vmax=diff_ranges[i])
    if i==0:
        plt.title('Difference (JULES_QS - JULES_BS)')
    plt.colorbar(fraction=0.05)
        

plt.show()


   
