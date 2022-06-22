#!/usr/bin/python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import netcdftime as nctime


dir_cruncep_v4='/prj/ALANIS/jules_data/cru_ncep/0.5deg/meteo/v4/6hr/'
c_4_color = 'red'
c_4_name = 'cruncep_v4'

latlon_file='/users/eow/edwcom/CRUNCEP/cru_ncep_land.nc' 


params=['LWdown','Precip','PSurf','Qair','SWdown','Tair','Wind']

date='198010'

dict={}
for param in params:
    inf=nc.Dataset(dir_cruncep_v4+param+'/cruncep_'+param+'_'+date+'.nc','r')
    dict[param]=inf.variables[param][:]
    inf.close()


inf=nc.Dataset(latlon_file,'r')
lats=inf.variables['Latitude'][:]
lons=inf.variables['Longitude'][:]
inf.close()

poopoint=np.where( (lats==29.75)&(lons==47.75) )[0]


F=plt.figure()
plot_num=1

for param in params:
    AX = F.add_subplot(4,2,plot_num)
    AX.plot(dict[param][:,poopoint],ls='',marker='.')
    AX.set_title(param)
    plot_num+=1

