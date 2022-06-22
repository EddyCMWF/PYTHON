#!/usr/bin/env python

################################################################################
# 
# Program: SC_simulator.py
# 
# Python Script to plot simulated soil carbon estimate 
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
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
#

JULES_output_DIR = '/users/eow/edwcom/CRUNCEP/n96/JULES_output/'
JULES_Quickspin_file = 'JULES_v4.3_TRIFFID_RsQ10_GLOBAL.monthly_mean.nc'
JULES_Bigspin_file   = 'JULES_v4.3_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'

SCsimulator_DIR = '/users/eow/edwcom/SC_simulator/output/'
SC_BS_FULLSERIES_file = 'CRUNCEPn96_simulation_JULESoutput_BigSpin_LAST10_monthly_funcmean_FastSoilCarbon.nc'

n96_landindexfile='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/monthly/'

# Open and read the land index file
inf=nc.Dataset(n96_landindexfile,'r')
index=inf.variables['index_2D'][:]
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
inf.close()

# Open and read in simulated CS estimates
# 
# Single timestep simulation
inf=nc.Dataset(SCsimulator_DIR+SC_BS_FULLSERIES_file,'r')
SC_C_4pools_Q10t_1D=inf.variables['C_4pools_Q10t'][:]
SC_C_4pools_RothCt_1D=inf.variables['C_4pools_RothCt'][:]
inf.close()

SC_C_4pools_Q10t_2D=SC_C_4pools_Q10t_1D[:,index]
SC_C_4pools_RothCt_2D=SC_C_4pools_RothCt_1D[:,index]
for i in range(4):
    SC_C_4pools_Q10t_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_C_4pools_RothCt_2D[i,:,:][np.where(index.mask==True)]=-999.0

SC_C_4pools_Q10t_2D   = np.ma.masked_equal(SC_C_4pools_Q10t_2D,-999.)
SC_C_4pools_RothCt_2D = np.ma.masked_equal(SC_C_4pools_RothCt_2D,-999.)


# Open and read in the JULES BigSpin data
inf=nc.Dataset(JULES_output_DIR+JULES_Bigspin_file)
JULES_lats=inf.variables['latitude'][:].squeeze()
JULES_lons=inf.variables['longitude'][:].squeeze()
JULES_time=nctime.num2date(inf.variables['time'][-120:],          \
                           units=inf.variables['time'].units, \
                           calendar='standard'                )

JULES_BS_SoilCarbon=inf.variables['cs'][-120:,:,:,:].squeeze()

inf.close()


JULES_BS_SoilCarbon_mean_1D = np.mean(JULES_BS_SoilCarbon,axis=0)

JULES_BS_SoilCarbon_mean_2D = JULES_BS_SoilCarbon_mean_1D[:,index]

for i in range(4):
    JULES_BS_SoilCarbon_mean_2D[i,:,:][np.where(index.mask==True)]=-999.0

JULES_BS_SoilCarbon_mean_2D=np.ma.masked_equal(JULES_BS_SoilCarbon_mean_2D,-999.)


# Plots to check if big spin up has spun up
# Plot time series for different regions

EUROPE_index  = np.where( (JULES_lats>=35.) & (JULES_lats<=70.) & \
                          (JULES_lons>=350.) | (JULES_lons<=45.))[0]
SIBERIA_index = np.where( (JULES_lats>=45.) & \
                          (JULES_lons>=45.) & (JULES_lons<=180.))[0]
AFRICA_index  = np.where( (JULES_lats>=-45.) & (JULES_lats<=35.) & \
                          (JULES_lons>=345.) | (JULES_lons<=50.))[0]
ASIA_index    = np.where( (JULES_lats>=0.) & (JULES_lats<=45.) & \
                          (JULES_lons>=50.) & (JULES_lons<=180.))[0]
NAMERI_index  = np.where( (JULES_lats>=40.) & (JULES_lats<=90.) & \
                          (JULES_lons>=240.) & (JULES_lons<=300.))[0]
CAMERI_index  = np.where( (JULES_lats>=-20.) & (JULES_lats<=40.) & \
                          (JULES_lons>=240.) & (JULES_lons<=300.))[0]
SAMERI_index  = np.where( (JULES_lats>=-50.) & (JULES_lats<=-20.) & \
                          (JULES_lons>=240.) & (JULES_lons<=330.))[0]

GEO_indexes = [ EUROPE_index, SIBERIA_index, AFRICA_index, ASIA_index, \
                NAMERI_index, CAMERI_index, SAMERI_index ]
GEO_names   = [ 'Europe', 'Siberia', 'Africa', 'Asia', \
                'North Ame', 'Central Ame', 'South Ame' ]

SC_pools=['DPM','RPM','BIO','HUM']


fig=plt.figure(figsize=(20,10))





plt.subplot(1,3,1)
plt.imshow(np.sum(SC_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0,vmax=30)
plt.title('Simulated SC - 4 Pools')
plt.colorbar(fraction=0.05)

plt.subplot(1,3,2)
plt.imshow(np.sum(JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',vmin=0,vmax=30)
plt.title('JULES-BigSpin - 4 Pools')
plt.colorbar(fraction=0.05)

plt.subplot(1,3,3)
plt.imshow(np.sum(SC_C_4pools_Q10t_2D-JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',cmap='RdBu',vmin=-15,vmax=15)
plt.title('Difference (Sim - JULES)')
plt.colorbar(fraction=0.05)

plt.show()


ABS_maxes = [ 0.3, 10, 1.2, 30. ]
diff_ranges= [ 0.15, 7, 0.8, 15. ]

for i in range(4):
    
    plt.subplot(4,3,(i*3)+1)
    plt.imshow(SC_C_4pools_Q10t_2D[i,:,:],origin='bottom',vmin=0,vmax=ABS_maxes[i])
    if i==0:
        plt.title('Simulated SC')
    plt.ylabel(SC_pools[i])
    plt.colorbar(fraction=0.05)
    
    plt.subplot(4,3,(i*3)+2)
    plt.imshow(JULES_BS_SoilCarbon_mean_2D[i,:,:],origin='bottom',vmin=0,vmax=ABS_maxes[i])
    if i==0:
        plt.title('JULES SC')
    plt.colorbar(fraction=0.05)
        
    
    plt.subplot(4,3,(i*3)+3)
    plt.imshow(SC_C_4pools_Q10t_2D[i,:,:]-JULES_BS_SoilCarbon_mean_2D[i,:,:],origin='bottom', \
               cmap='RdBu',vmin=(0-diff_ranges[i]),vmax=diff_ranges[i])
    if i==0:
        plt.title('Difference (Sim - JULES)')
    plt.colorbar(fraction=0.05)
        
plt.show()


