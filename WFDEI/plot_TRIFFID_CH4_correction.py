#!/bin/env python

import netCDF4 as nc
import numpy as np
from PlotTools import plot_tools as PT

data_dir='/users/eow/edwcom/WFD_EI/JULES_output/'
plot_dir=data_dir+'plots/'

good_file='JULES_WFDEI_nti_TRIFFID_gridded_daily_soil.2002.nc'
bad_file='JULES_WFDEI_nti_TRIFFID_gridded_daily_soil.2002.nc.precorrected'

good_inf=nc.Dataset(good_file,'r')
good_fwetl=good_inf.variables['fwetl'][:]
good_fch4=good_inf.variables['fch4_wetl'][:]
fch4_units='$'+good_inf.variables['fch4_wetl'].units+'$'
lat=good_inf.variables['latitude'][:]
lon=good_inf.variables['longitude'][:]
good_inf.close()
lon_2d,lat_2d=np.meshgrid(lon,lat)

bad_inf=nc.Dataset(bad_file,'r')
bad_fwetl=bad_inf.variables['fwetl'][:]
bad_fch4=bad_inf.variables['fch4_wetl'][:]
bad_inf.close()

good_fch4=np.max(good_fch4,axis=0)
bad_fch4=np.max(bad_fch4,axis=0)



PT.plot_map(good_fch4,lon_2d,lat_2d, \
            DATA_RANGE=[0.,0.5],\
            COLOURS=['aliceblue','lightgreen','orange'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11, \
            PLOT_TITLE='Corrected TRIFFID fch4-wetl (maximum 2002)',CBAR_LABEL=fch4_units, \
            FILE_PLOT=plot_dir+'Corrected_TRIFFID_fch4_wetl_max2002.png', \
            )


PT.plot_map(bad_fch4,lon_2d,lat_2d, \
            DATA_RANGE=[0.,0.5],\
            COLOURS=['aliceblue','lightgreen','orange'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11, \
            PLOT_TITLE='Uncorrected TRIFFID fch4-wetl (maximum 2002)',CBAR_LABEL=fch4_units, \
            FILE_PLOT=plot_dir+'Uncorrected_TRIFFID_fch4_wetl_max2002.png', \
            )
