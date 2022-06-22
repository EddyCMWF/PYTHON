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
JULES_Bigspin_file   = 'JULES_v4.3_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'

SCsimulator_DIR = '/users/eow/edwcom/SC_simulator/output/'
SC_BS_1ts_file     = 'CRUNCEPn96_simulation_JULESoutput_BigSpin_LAST10_FastSoilCarbon.nc'

n96_landindexfile='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/single_timestep/'

# Open and read the land index file
inf=nc.Dataset(n96_landindexfile,'r')
index=inf.variables['index_2D'][:]
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
inf.close()

# Single timestep simulation
inf=nc.Dataset(SCsimulator_DIR+SC_BS_1ts_file,'r')
SC_BS_1ts_C_4pools_Q10t_1D=inf.variables['C_4pools_Q10t'][:]
SC_BS_1ts_C_4pools_RothCt_1D=inf.variables['C_4pools_RothCt'][:]
SC_BS_Local_Litterfall_C_1D=inf.variables['Local_Litterfall_C'][:]
SC_BS_Litterfall_C_1D=inf.variables['Litterfall_C'][:]
SC_BS_Total_Litterfall_C_1D=inf.variables['Total_Litterfall_C'][:]

SC_BS_cveg_1D=inf.variables['c_veg'][:]

SC_BS_PFTfrac_1D=inf.variables['PFT_frac'][:]

SC_BS_LAI_1D=inf.variables['LAI'][:]

SC_BS_Soil_Resp_Fact_Q10t_1D=inf.variables['Soil_Resp_Fact_Q10t'][:]

SC_BS_Lit_SCpools_1D=inf.variables['Lit_SCpools'][:]

inf.close()

#SC_BS_Litterfall_C_1D=np.sum(SC_BS_Litterfall_C_1D,axis=0)

SC_BS_1ts_C_4pools_Q10t_2D=SC_BS_1ts_C_4pools_Q10t_1D[:,index]
SC_BS_1ts_C_4pools_RothCt_2D=SC_BS_1ts_C_4pools_RothCt_1D[:,index]
SC_BS_Litterfall_C_2D=SC_BS_Litterfall_C_1D[:,index]
SC_BS_Total_Litterfall_C_2D=SC_BS_Total_Litterfall_C_1D[index]
SC_BS_Local_Litterfall_C_2D=SC_BS_Local_Litterfall_C_1D[:,index]
SC_BS_Lit_SCpools_2D=SC_BS_Lit_SCpools_1D[:,index]
SC_BS_cveg_2D=SC_BS_cveg_1D[:,index]
SC_BS_PFTfrac_2D=SC_BS_PFTfrac_1D[:,index]
SC_BS_LAI_2D=SC_BS_LAI_1D[:,index]
SC_BS_Soil_Resp_Fact_Q10t_2D=SC_BS_Soil_Resp_Fact_Q10t_1D[:,index]


for i in range(4):
    SC_BS_1ts_C_4pools_Q10t_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_1ts_C_4pools_RothCt_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_Soil_Resp_Fact_Q10t_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_Lit_SCpools_2D[i,:,:][np.where(index.mask==True)]=-999.0

for i in range(5):
    SC_BS_Litterfall_C_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_Local_Litterfall_C_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_cveg_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_PFTfrac_2D[i,:,:][np.where(index.mask==True)]=-999.0
    SC_BS_LAI_2D[i,:,:][np.where(index.mask==True)]=-999.0

SC_BS_Total_Litterfall_C_2D[np.where(index.mask==True)]=-999.0

SC_BS_1ts_C_4pools_Q10t_2D = np.ma.masked_equal(SC_BS_1ts_C_4pools_Q10t_2D,-999.)
SC_BS_1ts_C_4pools_RothCt_2D = np.ma.masked_equal(SC_BS_1ts_C_4pools_RothCt_2D,-999.)
SC_BS_Litterfall_C_2D = np.ma.masked_equal(SC_BS_Litterfall_C_2D,-999.)
SC_BS_Total_Litterfall_C_2D = np.ma.masked_equal(SC_BS_Total_Litterfall_C_2D,-999.)
SC_BS_Local_Litterfall_C_2D = np.ma.masked_equal(SC_BS_Local_Litterfall_C_2D,-999.)
SC_BS_Lit_SCpools_2D = np.ma.masked_equal(SC_BS_Lit_SCpools_2D,-999.)
SC_BS_cveg_2D = np.ma.masked_equal(SC_BS_cveg_2D,-999.)
SC_BS_PFTfrac_2D = np.ma.masked_equal(SC_BS_PFTfrac_2D,-999.)
SC_BS_LAI_2D = np.ma.masked_equal(SC_BS_LAI_2D,-999.)
SC_BS_Soil_Resp_Fact_Q10t_2D= np.ma.masked_equal(SC_BS_Soil_Resp_Fact_Q10t_2D,-999.0)


# Open and read in the JULES BigSpin data
inf=nc.Dataset(JULES_output_DIR+JULES_Bigspin_file)
JULES_lats=inf.variables['latitude'][:].squeeze()
JULES_lons=inf.variables['longitude'][:].squeeze()
JULES_time=nctime.num2date(inf.variables['time'][:],          \
                           units=inf.variables['time'].units, \
                           calendar='standard'                )

JULES_BS_SoilCarbon   = inf.variables['cs'][:,:,:,:].squeeze()
JULES_BS_Litterfall_C = inf.variables['lit_c'][:,:,:,:].squeeze()
JULES_BS_cveg = inf.variables['c_veg'][:,:,:,:].squeeze()
JULES_BS_resp_s = inf.variables['resp_s'][:,:,:,:].squeeze()

inf.close()

JULES_BS_SoilCarbon_total= np.sum(JULES_BS_SoilCarbon,axis=1)
JULES_BS_SoilCarbon_mean = np.mean(JULES_BS_SoilCarbon,axis=0)

JULES_BS_Litterfall_C_total= np.sum(JULES_BS_Litterfall_C,axis=1)
JULES_BS_Litterfall_C_mean = np.mean(JULES_BS_Litterfall_C,axis=0)

JULES_BS_cveg_total= np.sum(JULES_BS_cveg,axis=1)
JULES_BS_cveg_mean = np.mean(JULES_BS_cveg,axis=0)

JULES_BS_resp_s_total= np.sum(JULES_BS_resp_s,axis=1)
JULES_BS_resp_s_mean = np.mean(JULES_BS_resp_s,axis=0)

JULES_BS_SoilCarbon_mean_2D=JULES_BS_SoilCarbon_mean[:,index]
JULES_BS_Litterfall_C_mean_2D=JULES_BS_Litterfall_C_mean[:,index]
JULES_BS_cveg_mean_2D=JULES_BS_cveg_mean[:,index]
JULES_BS_resp_s_mean_2D=JULES_BS_resp_s_mean[:,index]

for i in range(4):
    JULES_BS_SoilCarbon_mean_2D[i,:,:][np.where(index.mask==True)]=-999.0
    JULES_BS_resp_s_mean_2D[i,:,:][np.where(index.mask==True)]=-999.0

for i in range(5):
    JULES_BS_Litterfall_C_mean_2D[i,:,:][np.where(index.mask==True)]=-999.0
    JULES_BS_cveg_mean_2D[i,:,:][np.where(index.mask==True)]=-999.0

JULES_BS_SoilCarbon_mean_2D=np.ma.masked_equal(JULES_BS_SoilCarbon_mean_2D,-999.)
JULES_BS_Litterfall_C_mean_2D=np.ma.masked_equal(JULES_BS_Litterfall_C_mean_2D,-999.)
JULES_BS_cveg_mean_2D=np.ma.masked_equal(JULES_BS_cveg_mean_2D,-999.)
JULES_BS_resp_s_mean_2D=np.ma.masked_equal(JULES_BS_resp_s_mean_2D,-999.)

JULES_BS_Soil_Resp_Fact_2D = JULES_BS_resp_s_mean_2D/JULES_BS_SoilCarbon_mean_2D



plt.subplot(1,2,1)
plt.imshow(SC_BS_Lit_SCpools_2D[0,:,:],origin='bottom')
plt.title('Litter in DPM pool')
plt.colorbar()

plt.subplot(1,2,2)
plt.imshow(SC_BS_Lit_SCpools_2D[0,:,:],origin='bottom')
plt.title('Litter in RPM pool')
plt.colorbar()

plt.show()


ncols=2
nrows=4

for i in range(nrows):
    plt.subplot(nrows,ncols,(i*ncols)+1)
    plt.imshow(SC_BS_1ts_C_4pools_Q10t_2D[i,:,:],origin='bottom')
    plt.colorbar()
    if i==0:
        plt.title('Simulated SC by pool')

    plt.subplot(nrows,ncols,(i*ncols)+2)
    plt.imshow(JULES_BS_SoilCarbon_mean_2D[i,:,:],origin='bottom')
    plt.colorbar()
    if i==0:
        plt.title('JULES SC by pool')


plt.show()



ncols=5
nrows=5

for i in range(nrows):
    plt.subplot(nrows,ncols,(i*ncols)+1)
    plt.imshow(SC_BS_Litterfall_C_2D[i,:,:],origin='bottom')
    if i==0:
        plt.title('Total Litterfall')
    plt.colorbar()
    
    plt.subplot(nrows,ncols,(i*ncols)+2)
    plt.imshow(SC_BS_Local_Litterfall_C_2D[i,:,:],origin='bottom')
    if i==0:
        plt.title('Local Litterfall')
    plt.colorbar()
    
    plt.subplot(nrows,ncols,(i*ncols)+3)
    plt.imshow(SC_BS_cveg_2D[i,:,:],origin='bottom')
    if i==0:
        plt.title('Vegetation Carbon')
    plt.colorbar()

    plt.subplot(nrows,ncols,(i*ncols)+4)
    plt.imshow(SC_BS_PFTfrac_2D[i,:,:],origin='bottom',vmin=0,vmax=1)
    if i==0:
        plt.title('PFT Fractional Cover')
    plt.colorbar()

    plt.subplot(nrows,ncols,(i*ncols)+5)
    plt.imshow(SC_BS_LAI_2D[i,:,:],origin='bottom')
    if i==0:
        plt.title('LAI')
    plt.colorbar()

plt.show()



plt.subplot(1,3,1)
plt.imshow(np.sum(SC_BS_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0,vmax=50.)
plt.title('Soil C - BigSpin_Simulator')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(np.sum(JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',vmin=0,vmax=50.)
#plt.imshow(np.sum(SC_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0)
plt.title('Soil - JULES-BigSpin')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(np.sum(SC_BS_1ts_C_4pools_Q10t_2D-JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',cmap='RdBu',vmin=-50,vmax=50)
plt.title('Difference (Simulator-JULES)')
plt.colorbar()

plt.show()



plt.subplot(1,3,1)
plt.imshow(np.sum(SC_BS_Litterfall_C_2D,axis=0),origin='bottom',vmin=0)
plt.title('Litterfall - BigSpin_Simulator')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(np.sum(JULES_BS_Litterfall_C_mean_2D,axis=0),origin='bottom',vmin=0)
#plt.imshow(np.sum(SC_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0)
plt.title('Litterfall - JULES-BigSpin')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(np.sum(SC_BS_Litterfall_C_2D-JULES_BS_Litterfall_C_mean_2D,axis=0),origin='bottom',cmap='RdBu',vmin=-5,vmax=5)
plt.title('Difference (Simulator-JULES)')
plt.colorbar()

plt.show()



plt.subplot(1,3,1)
plt.imshow(np.sum(SC_BS_Soil_Resp_Fact_Q10t_2D,axis=0)*(3600.*24.*365.),origin='bottom',vmin=0,vmax=30)
plt.title('SoilRespFactor - BigSpin_Simulator')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(np.sum(JULES_BS_Soil_Resp_Fact_2D,axis=0)*(3600.*24.*365.),origin='bottom',vmin=0,vmax=30)
plt.title('SoilRespFactor - JULES-BigSpin')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(np.sum(SC_BS_Soil_Resp_Fact_Q10t_2D*(3600.*24.*365.)-JULES_BS_Soil_Resp_Fact_2D*(3600.*24.*365.),axis=0),origin='bottom',cmap='RdBu',vmin=-10,vmax=10)
plt.title('Difference (Simulator-JULES)')
plt.colorbar()

plt.show()



ncols=5
nrows=1
i=0

plt.subplot(nrows,ncols,(i*ncols)+1)
plt.imshow(np.sum(SC_BS_Litterfall_C_2D,axis=0),origin='bottom')
plt.title('Total Litterfall')
plt.colorbar()

plt.subplot(nrows,ncols,(i*ncols)+2)
plt.imshow(np.sum(SC_BS_Local_Litterfall_C_2D,axis=0),origin='bottom')
plt.title('Local Litterfall')
plt.colorbar()

plt.subplot(nrows,ncols,(i*ncols)+3)
plt.imshow(np.sum(SC_BS_cveg_2D,axis=0),origin='bottom')
plt.title('Vegetation Carbon')
plt.colorbar()

plt.subplot(nrows,ncols,(i*ncols)+4)
plt.imshow(np.sum(SC_BS_PFTfrac_2D,axis=0),origin='bottom',vmin=0,vmax=1)
plt.title('PFT Fractional Cover')
plt.colorbar()

plt.subplot(nrows,ncols,(i*ncols)+5)
plt.imshow(np.sum(SC_BS_LAI_2D,axis=0),origin='bottom')
plt.title('LAI')
plt.colorbar()

plt.show()
