#!/bin/env python3
import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

from PlotTools import plot_tools as PTs

L_sarah=False
Cs_EQ_dir='/prj/CLIFFTOP/Cs_Equilibrium/'
SARAH_file = Cs_EQ_dir+'carbon_withfrac.txt'   # 'test.txt' 
JULES_file = Cs_EQ_dir+'vn4.8_imogen.spinup_2000.dump.1850_zeroDtemp_zeroWP.nc'
ECP_file = Cs_EQ_dir+'/JULES_output/Equilibrium_Soil_Carbon.nc'
ANCILS_DIR = '/prj/CLIFFTOP/COMMON_DATA/ANCILS/'
plot_dir= Cs_EQ_dir+'plots/'

grid_file=ANCILS_DIR+'grid_info.nc'
grinf=nc.Dataset(grid_file,'r')
lats_2d = grinf.variables['latitude'][:]
lons_2d = grinf.variables['longitude'][:]
Area_2d = grinf.variables['Area'][:]       # I don't actually use this but it's here
land_index = grinf.variables['land_index'][:]
grinf.close()

# 1Dimension grid cell area data for calculating totals etc.
AREA_file=ANCILS_DIR+'Area_in_iris_format.nc'
Ainf=nc.Dataset(AREA_file,'r')
AREA_1D = Ainf.variables['area'][:].squeeze()
lats_1D = Ainf.variables['latitude'][:].squeeze()
lons_1D = Ainf.variables['longitude'][:].squeeze()
Ainf.close()

dz = np.array([0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293,0.21517585,
               0.23784142,0.25980762,0.28117066,0.30200527,0.32237098,0.34231625,0.36188121])
z = dz.cumsum()
nz = len(dz)

pools = ['DPM','RPM','BIO','HUM']
npools=len(pools)

if L_sarah:
    S_lines = open(SARAH_file).readlines()
    S_data = []
    for line in S_lines:
        split=line.replace('\n','').split(',')
        S_data.append([ float(cs) for cs in split ])
    
    S_data=np.array(S_data)
    S_data=S_data.reshape(1631,npools,nz).transpose(1,2,0)
    print(S_data[0,:])
    S_data[np.isnan(S_data)]=-999.
    S_data = np.ma.masked_less_equal(S_data,0)
    for iz in range(nz):
        S_data[:,iz,:] *= dz[iz]
    S_data[S_data<0.] = 0.

else:
    Einf=nc.Dataset(ECP_file,'r')
    S_data = Einf.variables['cs'][:].squeeze()
    Einf.close()
    S_data[np.isnan(S_data)]=1e-6
    plot_dir= Cs_EQ_dir+'ECP_plots/'

os.system('mkdir -p '+plot_dir)
Jinf = nc.Dataset(JULES_file,'r')
J_data = Jinf.variables['cs'][:].squeeze()
Jinf.close()

J_data = np.ma.masked_less_equal(J_data,0)


DIFF=S_data-J_data


Colours=('#f9f9ea','#fff68f','#ffa500','#5f1205')
pool_ranges = [ (0,0.2),(0,15.),(0,1.),(0,50) ]
tick_formats= [ '%.2f','%.1f','%.1f','%.1f' ]
fig,axes=plt.subplots(ncols=npools,nrows=1,figsize=[25,5])
for ipool in range(npools):
    ax=axes[ipool]
    plot_data = J_data[ipool,:].sum(axis=0)
    print(plot_data.shape)
    plot_data = np.ma.masked_array(plot_data[land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=6,
                 RESOLUTION='c', MAP_TYPE='Mesh',#SET_OVER='#5f1205',
                 AXIS=ax,TICK_FORMAT=tick_formats[ipool], CBAR_TICK_LENGTH=4,
                 CBAR_ORIENTATION='horizontal',FONTSIZES=[12,12,14,18]
                 )

    axes[ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
fig.savefig(plot_dir+'JULES_SoilC_Maps_byPool.png',bbox_inches='tight')

plt.show()

fig,axes=plt.subplots(ncols=npools,nrows=1,figsize=[25,5])
for ipool in range(npools):
    ax=axes[ipool]
    plot_data = S_data[ipool,:].sum(axis=0)
    plot_data = np.ma.masked_array(plot_data[land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=6,
                 RESOLUTION='c', MAP_TYPE='Mesh',#SET_OVER='#5f1205',
                 AXIS=ax,TICK_FORMAT=tick_formats[ipool], CBAR_TICK_LENGTH=4,
                 CBAR_ORIENTATION='horizontal',FONTSIZES=[12,12,14,18]
                 )

    axes[ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
fig.savefig(plot_dir+'SARAH_SoilC_Maps_byPool.png',bbox_inches='tight')

plt.show()


Colours=('blue','white','red')
pool_ranges = [ (-0.1,0.1),(-5.,5.),(-0.5,0.5),(-15,15) ]

fig,axes=plt.subplots(ncols=npools,nrows=1,figsize=[25,5])
for ipool in range(npools):
    ax=axes[ipool]
    plot_data = DIFF[ipool,:].sum(axis=0)
    plot_data = np.ma.masked_array(plot_data[land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=6,
                 RESOLUTION='c', MAP_TYPE='Mesh',SET_OVER='red',SET_UNDER='blue',
                 AXIS=ax,TICK_FORMAT=tick_formats[ipool], CBAR_TICK_LENGTH=4,
                 CBAR_ORIENTATION='horizontal',FONTSIZES=[12,12,14,18]
                 )

    axes[ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
fig.suptitle('Cs Difference - red=Sarah larger, blue=JULES larger',fontsize=30)
fig.savefig(plot_dir+'DIFFERNCE_SoilC_Maps_byPool.png',bbox_inches='tight')

plt.show()


#quit()


J_GlobalTot = np.sum(J_data*AREA_1D,axis=2)
S_GlobalTot = np.sum(S_data*AREA_1D,axis=2)
print(S_GlobalTot)
DIFF_GlobalTot = np.sum(DIFF*AREA_1D,axis=2)

fig,axes=plt.subplots(ncols=npools,nrows=1,figsize=[15,7.5])

for ipool in range(npools):
    ax = axes[ipool]
    ax.plot(J_GlobalTot[ipool,:]*1e-12,z,c='b',label='JULES',lw=2)
    ax.plot(S_GlobalTot[ipool,:]*1e-12,z,c='r',label='SARAH',lw=2)
    ax.plot(DIFF_GlobalTot[ipool,:]*1e-12,z,c='g',label='Difference',lw=2)
    ax.set_title(pools[ipool],fontsize=16)
    ax.set_xlabel('Global Total Soil Carbon (GtC)')
    if ipool==0:
        ax.set_ylabel('Depth (m)')
    ax.invert_yaxis()

handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc=8,fontsize=16, ncol=3)
fig.savefig(plot_dir+'SoilC_GlobalTotal_withDepth_byPool.png',bbox_inches='tight')

plt.show()

J_GlobalMean= np.mean(J_data,axis=2)
S_GlobalMean= np.mean(S_data,axis=2)
DIFF_GlobalMean= np.mean(DIFF,axis=2)

fig,axes=plt.subplots(ncols=npools,nrows=1,figsize=[15,7.5])

for ipool in range(npools):
    ax = axes[ipool]
    ax.plot(J_GlobalMean[ipool,:],z,c='b',label='JULES',lw=2)
    ax.plot(S_GlobalMean[ipool,:],z,c='r',label='SARAH',lw=2)
    ax.plot(DIFF_GlobalMean[ipool,:],z,c='g',label='Difference',lw=2)
    ax.set_title(pools[ipool],fontsize=16)
    ax.set_xlabel('Global Mean Soil Carbon (kgC m$^{-2}$)')
    if ipool==0:
        ax.set_ylabel('Depth (m)')
    ax.invert_yaxis()

handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc=8,fontsize=16, ncol=3)
fig.savefig(plot_dir+'SoilC_GlobalMean_withDepth_byPool.png',bbox_inches='tight')

plt.show()

quit()

Colours=('#f9f9ea','#fff68f','#ffa500','#5f1205')
pool_ranges = [ (0,0.03),(0,3.),(0,0.2),(0,10) ]

fig,axes=plt.subplots(ncols=npools,nrows=nz,figsize=[20,35])
for ipool in range(npools):
  for iz in range(nz):
    ax=axes[iz,ipool]
    plot_data = np.ma.masked_array(J_data[ipool,iz,land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=11,
                 RESOLUTION='c', MAP_TYPE='Mesh',SET_OVER='#5f1205',
                 AXIS=ax,
                 CBAR_ORIENTATION='off',FONTSIZES=[12,12,14,18]
                 )

  axes[0,ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
  TickLEVELS = np.linspace(data_range[0],data_range[1],num=6)
  CBAR=plt.colorbar(IMAGE,ax=axes[:,ipool].flatten().tolist(),orientation='horizontal',
                    ticks=TickLEVELS,pad=0.01,fraction=0.05) 
fig.savefig(plot_dir+'JULES_SoilC_Maps_byDepthandPool.png',bbox_inches='tight')

plt.show()

fig,axes=plt.subplots(ncols=npools,nrows=nz,figsize=[20,35])
for ipool in range(npools):
  for iz in range(nz):
    ax=axes[iz,ipool]
    plot_data = np.ma.masked_array(S_data[ipool,iz,land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=11,
                 RESOLUTION='c', MAP_TYPE='Mesh',SET_OVER='#5f1205',
                 AXIS=ax,
                 CBAR_ORIENTATION='off',FONTSIZES=[12,12,14,18]
                 )

  axes[0,ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
  TickLEVELS = np.linspace(data_range[0],data_range[1],num=6)
  CBAR=plt.colorbar(IMAGE,ax=axes[:,ipool].flatten().tolist(),orientation='horizontal',
                    ticks=TickLEVELS,pad=0.01,fraction=0.05) 
fig.savefig(plot_dir+'SARAH_SoilC_Maps_byDepthandPool.png',bbox_inches='tight')

plt.show()


Colours=('blue','white','red')
pool_ranges = [ (-0.02,0.02),(-2.,2.),(-0.1,0.1),(-5,5) ]

fig,axes=plt.subplots(ncols=npools,nrows=nz,figsize=[20,35])
for ipool in range(npools):
  for iz in range(nz):
    ax=axes[iz,ipool]
    plot_data = np.ma.masked_array(DIFF[ipool,iz,land_index],mask=land_index.mask)
    data_range=pool_ranges[ipool]
    IMAGE=PTs.plot_map(plot_data,lons_2d,lats_2d,
                 DATA_RANGE=data_range,
                 COLOURS=Colours,INTERPOLATE_COLOURS=True,NLEVELS=25,NTICKS=11,
                 RESOLUTION='c', MAP_TYPE='Mesh',SET_OVER='red',SET_UNDER='blue',
                 AXIS=ax,
                 CBAR_ORIENTATION='off',FONTSIZES=[12,12,14,18]
                 )

  axes[0,ipool].text(0.5,100,pools[ipool],ha='center',va='bottom',fontsize=16)
  TickLEVELS = np.linspace(data_range[0],data_range[1],num=6)
  CBAR=plt.colorbar(IMAGE,ax=axes[:,ipool].flatten().tolist(),orientation='horizontal',
                    ticks=TickLEVELS,pad=0.01,fraction=0.05) 
fig.suptitle('Cs Difference - red=Sarah larger, blue=JULES larger',fontsize=30)
fig.savefig(plot_dir+'DIFFERNCE_SoilC_Maps_byDepthandPool.png',bbox_inches='tight')

plt.show()

