#!/usr/bin/python
#
# Python routine to plot methan time-series from JULES ncdf output
#
# Edward Comyn-Platt
# CEH Wallingford
# Jan 2015
# 
# Contains
#
import os
import sys
import numpy as np
#
import data_netCDF
from matplotlib import pyplot as plt
from matplotlib import dates  as mdates
#
dir      = '/users/eow/edwcom/JULES/test_runs/z_output/'
in_file  = 'TEST_SERIAL_job001_NCEP_CRU_v4_dom0001.monthly.nc'
plotname = 'Time_Series.png'


lat_dims,lats = data_netCDF.data_netCDF_array_var(dir+in_file,'latitude')
lon_dims,lons = data_netCDF.data_netCDF_array_var(dir+in_file,'longitude')
t_dims,time   = data_netCDF.data_netCDF_array_var(dir+in_file,'time')

#convert from seconds since start date to days since 0001-01-01 for mpl compatibility
time=(time/(3600.*24.)) + mdates.strpdate2num('%Y-%m-%d')('1980-01-01')

f_ch4_dims, f_ch4 = data_netCDF.data_netCDF_array_var(dir+in_file,'fch4_wetl')

f_ch4_TS = []

for i in range(0,t_dims[0]):
    f_ch4_TS.append(np.mean(f_ch4[i,0,:]))


plt.plot_date(x=time,y=f_ch4_TS)
plt.show()








