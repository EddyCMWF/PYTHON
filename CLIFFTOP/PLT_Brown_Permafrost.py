#!/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
#import pyproj

#Create Basemap projection for converting coordinates to lat/lon
#lambert_aea = {'proj': 'laea', 'lat_0':0, 'lon_0':0, 'x_0':0., 'y_0':0.,
#                'ellps': 'WGS84', 'datum': 'WGS84', 'R':6378137.0}
#lambProj = pyproj.Proj(lambert_aea)

M = Basemap(projection='nplaea',resolution='c',lon_0=0,boundinglat=1)
#M = Basemap(projection='laea',resolution='c',lat_0=90,lon_0=0)

infile='/prj/CLIFFTOP/PermaFrost/nhipa.nc'

inf=nc.Dataset(infile,'a')

data=inf.variables['nhipa'][:]
fill_value=inf.variables['nhipa']._FillValue
easting=inf.variables['easting'][:]
northing=inf.variables['northing'][:]

continuous_mask = (data==1)|(data==5)|(data==9)|(data==13)|(data==17)|(data==22)#|(data==21)
#discontinuous_mask = (data==2)|(data==3)|(data==6)|(data==7)|(data==10)|(data==11)|(data==14)|(data==15)|(data==18)|(data==19)
discontinuous_mask = (data==2)|(data==6)|(data==10)|(data==14)|(data==18)

data_minmax = np.zeros_like(data)
data_minmax[continuous_mask] = 2.
data_minmax[discontinuous_mask] = 1.
data_minmax = np.ma.masked_array(data_minmax,mask=(data==0)) #|(data==24))

data_minmax_grad = np.gradient(data_minmax)

data_cont = np.zeros_like(data)
data_cont[continuous_mask] = 1.
cont_lims=np.gradient(data_cont)
cont_lim = np.abs(cont_lims[0])+np.abs(cont_lims[1])
cont_lim[cont_lim>0.] = 1.


data_discont = np.zeros_like(data)
data_discont[discontinuous_mask] = 1.
discont_lims=np.gradient(data_discont)
discont_lim = np.abs(discont_lims[0])+np.abs(discont_lims[1])
discont_lim[discont_lim>0.] = 2.

data_minmax = cont_lim + discont_lim
data_minmax[data_minmax==3]=1.
data_minmax = np.ma.masked_array(data_minmax,mask=(data_minmax==0),fill_value=fill_value) #|(data==24))
data_minmax.data[data_minmax.mask==True]=fill_value

fig,axes=plt.subplots(ncols=2,nrows=1,figsize=(18,6))
raw_image=axes[0].imshow(data)
plt.colorbar(raw_image,ax=axes[0])
line_image=axes[1].imshow(data_minmax)
plt.colorbar(line_image,ax=axes[1])
plt.show()

# Save edge data field
var='perma_limits'
if var not in inf.variables:
    outvar=inf.createVariable(var,'float32',('northing','easting'),fill_value=fill_value)
    outvar.long_name='Continuous and discontinuous permafrost limit'
    outvar.units='-'
    outvar.values=[1,2]
    outvar.value_meanings='1:Continuous Permafrost Boundary; 2:Discontinuous Permafrost Boundary'
else:
    outvar=inf.variables[var]
outvar[:]=data_minmax



easting2D,northing2D=np.meshgrid(easting-easting[0],northing-northing[-1])
lon,lat = M(easting2D,northing2D,inverse=True)
#lat2=np.copy(lat)
#lat2[data==0]=fill_value
#lat2=np.ma.masked_equal(lat2,fill_value)

var='latitude'
if var not in inf.variables:
    outvar=inf.createVariable('latitude','float32',('northing','easting'),fill_value=fill_value)
    outvar.units='Degree North'
else:
    outvar=inf.variables[var]
outvar[:]=lat

var='longitude'
if var not in inf.variables:
    outvar=inf.createVariable('longitude','float32',('northing','easting'),fill_value=fill_value)
    outvar.units='Degree East'
else:
    outvar=inf.variables[var]
outvar[:]=lon


inf.close()

