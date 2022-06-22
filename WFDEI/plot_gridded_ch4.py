#!/bin/env python

import netCDF4 as nc
import numpy as np
from PlotTools import plot_tools as PT

data_dir='/users/eow/edwcom/WFD_EI/JULES_output/'
filename='JULES_WFDEI_nti_TRIFFID_gridded_daily_ch4.2009.nc'
filename='JULES_WFDEI_nti_NG-HWSD_gridded_daily_ch4.2009.nc'

inf=nc.Dataset(data_dir+filename,'r')
fwetl=inf.variables['fwetl'][:]
fch4=inf.variables['fch4_wetl'][:]
fch4_units='$'+inf.variables['fch4_wetl'].units+'$'
lat=inf.variables['latitude'][:]
lon=inf.variables['longitude'][:]
inf.close()
lon_2d,lat_2d=np.meshgrid(lon,lat)

fch4=np.max(fch4,axis=0)
#fch4_units='$10^{-9} kg C m^{-2} s^{-1}$'

#region='Global'
region='Scandinavia'

plot_dir=data_dir+'plots/'
plot_filename=filename[:-3]+'_'+region+'.png'
PLOT_TITLE=filename[:-3].replace('_',' ').replace('.',' - ')+' - '+region

if region=='Global':
    lat_range=[-90,90]
    lon_range=[-180,180]
    data_range=[0,20]
elif region=='Scandinavia':
    lat_range=[55,75]
    lon_range=[0,40]
    data_range=[0,15]

PT.plot_map(fch4,lon_2d,lat_2d, \
            DATA_RANGE=data_range, LON_RANGE=lon_range, LAT_RANGE=lat_range, \
            COLOURS=['aliceblue','powderblue','lightgreen','#f2dd0b','orange'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11, TICK_FORMAT='%.1f', \
            PLOT_TITLE='fch4-wetl '+PLOT_TITLE,CBAR_LABEL=fch4_units, \
            FILE_PLOT=plot_dir+plot_filename, \
            )

