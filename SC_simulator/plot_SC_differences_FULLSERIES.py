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
#SC_monthly_file = 'CRUNCEPn96_simulation_JULESoutput_LAST10_monthly_FastSoilCarbon.nc'
SC_BS_FULLSERIES_file = 'CRUNCEPn96_simulation_JULESoutput_BigSpin_FULL_SERIES_FastSoilCarbon.nc'

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

SC_C_4pools_Q10t_2D=SC_C_4pools_Q10t_1D[:,:,index]
SC_C_4pools_RothCt_2D=SC_C_4pools_RothCt_1D[:,:,index]
for i in range(4):
    for j in range(12):
        SC_C_4pools_Q10t_2D[i,j,:,:][np.where(index.mask==True)]=-999.0
        SC_C_4pools_RothCt_2D[i,j,:,:][np.where(index.mask==True)]=-999.0

SC_C_4pools_Q10t_2D   = np.ma.masked_equal(SC_C_4pools_Q10t_2D,-999.)
SC_C_4pools_RothCt_2D = np.ma.masked_equal(SC_C_4pools_RothCt_2D,-999.)


# Open and read in the JULES BigSpin data
inf=nc.Dataset(JULES_output_DIR+JULES_Bigspin_file)
JULES_lats=inf.variables['latitude'][:].squeeze()
JULES_lons=inf.variables['longitude'][:].squeeze()
JULES_time=nctime.num2date(inf.variables['time'][:],          \
                           units=inf.variables['time'].units, \
                           calendar='standard'                )

JULES_BS_SoilCarbon=inf.variables['cs'][:,:,:,:].squeeze()

inf.close()


JULES_BS_SoilCarbon_mnth = JULES_BS_SoilCarbon[:,:,:]
JULES_BS_SoilCarbon_mnth = JULES_BS_SoilCarbon_mnth.transpose(1,0,2)

JULES_BS_SoilCarbon_mean_2D=JULES_BS_SoilCarbon_mnth[:,:,index]
for i in range(4):
    for j in range(12):
        JULES_BS_SoilCarbon_mean_2D[i,j,:,:][np.where(index.mask==True)]=-999.0

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

#fig = plt.figure()

for SCpool in range(4):
    print SC_pools[SCpool]
    plt.subplot(4,2,(SCpool*2)+1)
    for GEOindex,name in zip(GEO_indexes,GEO_names):
        #print name, GEOindex.shape
        #print SC_C_4pools_Q10t_1D[SCpool,:,GEOindex].shape
        SC_SoilCarbon_mean_TS=np.mean(SC_C_4pools_Q10t_1D[SCpool,:,GEOindex],axis=0)
        plt.plot(JULES_time,SC_SoilCarbon_mean_TS,label=name)
    #plt.ylim( (0,5) )
    plt.ylabel=SC_pools[SCpool]
    if SCpool==0:
        plt.title='Simulated Soil Carbon by Pool'
    #plt.set_ylabel(SC_pools[SCpool]+' (kgC)')
        
    plt.subplot(4,2,(SCpool*2)+2)
    for GEOindex,name in zip(GEO_indexes,GEO_names):
        #print name, GEOindex.shape
        #print JULES_BS_SoilCarbon_mnth[SCpool,:,GEOindex].shape
        JULES_BS_SoilCarbon_mean_TS=np.mean(JULES_BS_SoilCarbon_mnth[SCpool,:,GEOindex],axis=0)
        plt.plot(JULES_time,JULES_BS_SoilCarbon_mean_TS,label=name)
        
        print name, \
            np.mean(SC_C_4pools_Q10t_1D[SCpool,:,GEOindex]), \
            np.mean(JULES_BS_SoilCarbon_mnth[SCpool,:,GEOindex])
    
    
    #plt.ylim( (0,5) )

    if SCpool==0:
        plt.title='JULES Soil Carbon by Pool'
    #plt.set_ylabel(SC_pools[SCpool]+' (kgC)')
    
     

plt.show()




plt.subplot(1,3,1)
plt.imshow(np.sum(SC_BS_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0)
plt.title('BigSpin_Simulator')
plt.colorbar()

plt.subplot(1,3,2)
#plt.imshow(np.sum(JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',vmin=0)
plt.imshow(np.sum(SC_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',vmin=0)
plt.title('QuickSpin_Simulator')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(np.sum(SC_1ts_C_4pools_Q10t_2D-SC_BS_1ts_C_4pools_Q10t_2D,axis=0),origin='bottom',cmap='RdBu')
plt.title('Difference (QS-BS)')
plt.colorbar()

plt.show()







plt.subplot(1,3,1)
plt.imshow(np.sum(JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',vmin=0,vmax=25)
plt.title('JULES-BigSpin')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(np.sum(JULES_QS_SoilCarbon_mean_2D,axis=0),origin='bottom',vmin=0,vmax=25)
plt.title('JULES-QuickSpin')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(np.sum(JULES_QS_SoilCarbon_mean_2D-JULES_BS_SoilCarbon_mean_2D,axis=0),origin='bottom',cmap='RdBu',vmin=-10,vmax=10)
plt.title('Difference (QS-BS)')
plt.colorbar()

plt.show()


















colors=['y','g','orange','saddlebrown']
labels=['DPM','RPM','BIO','HUM']

fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothCt[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothCt[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

        ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()



fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10t[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10t[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()



lat_binsize=5.
minlat=-60.
maxlat=90.
lat_bins=np.arange(minlat,maxlat,lat_binsize)

Lit_SCpools_Zonal=np.zeros([2,len(lat_bins)])
Soil_Resp_Fact_Q10_Zonal=np.zeros([4,len(lat_bins)])
Soil_Resp_Fact_RothC_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_Q10_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_RothC_Zonal=np.zeros([4,len(lat_bins)])


for bin in range(len(lat_bins)):
    tempdex=np.where( (lats_JULES>=lat_bins[bin]        ) & \
                    (lats_JULES< lat_bins[bin]+lat_binsize) )[0]
    #
    Lit_SCpools_Zonal[:,bin]=np.mean(Lit_SCpools_ECP[:,tempdex],axis=1)
    Soil_Resp_Fact_Q10_Zonal[:,bin]=np.mean(Soil_Resp_Fact_Q10t[:,tempdex],axis=1)
    Soil_Resp_Fact_RothC_Zonal[:,bin]=np.mean(Soil_Resp_Fact_RothCt[:,tempdex],axis=1)
    C_4pools_Q10_Zonal[:,bin]=np.mean(C_4pools_Q10t[:,tempdex],axis=1)
    C_4pools_RothC_Zonal[:,bin]=np.mean(C_4pools_RothCt[:,tempdex],axis=1)




fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
#ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()


fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothC_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothC_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()

        
