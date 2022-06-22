#!/bin/env python2.7

import netCDF4 as nc
import numpy as np

grid_file='/prj/CLIFFTOP/COMMON_DATA/ANCILS/grid_info.nc'

inf=nc.Dataset(grid_file,'a')

lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]

lat_radias=lats*np.pi/180.
# calculate the area of the gid
rearth=6371000.
area=rearth*rearth*np.radians(3.75)*(np.sin(np.radians(lats+2.5))-np.sin(np.radians(lats)))

import math
#calculate the area of the gid
area2=np.zeros((56,96))
rearth=6371000.
#lats=inf.variables['latitude'][:]
#lons=inf.variables['longitude'][:]
for ilat in range(56):
    for ilon in range(96):
        area2[ilat,ilon]= rearth*rearth*math.radians(3.75)*\
                    (math.sin(math.radians(lats[ilat,ilon]+2.5))-math.sin(math.radians(lats[ilat,ilon])))


lats_above=lats+1.25
lats_below=lats-1.25
lons_east=lons+1.875
lons_west=lons-1.875

area3=2*np.pi*(rearth**2)*\
      np.abs(np.sin(np.radians(lats_above))-np.sin(np.radians(lats_below)))*\
      np.abs(lons_east-lons_west)


lat_degree_length = np.pi*rearth/180.

lat_gaps_m = (lats_above-lats_below)*lat_degree_length
lon_gaps_m_mid = (lons_east-lons_west)*lat_degree_length* \
                  np.cos(np.radians(lats))
lon_gaps_m = (lons_east-lons_west)*lat_degree_length* \
              (np.cos(np.radians(lats_above))+np.cos(np.radians(lats_below)))/2.0

area4=lat_gaps_m*lon_gaps_m


plt.subplot(1,3,1)
plt.imshow(area,origin='bottom',interpolation='none')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(area5,origin='bottom',interpolation='none')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow((area-area5)/area,origin='bottom',interpolation='none')
plt.colorbar()

plt.show()




outvar=inf.variable['Area']
outvar[:]=area4


inf.close()

