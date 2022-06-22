#!/bin/env python

import netCDF4 as nc
import numpy as np

Cdir='/users/eow/edwcom/CHESS/'

latlon_file = '/prj/chess/data/1km/v1.0/ancil/chess_soilparams_hwsd_bc.nc'
inf=nc.Dataset(latlon_file,'r')
lats_2D=inf.variables['lat'][:]
lons_2D=inf.variables['lon'][:]
inf.close()

site_loc_file=Cdir+'Active_EC_approx_locs_Xtra.csv'
site_loc_lines=open(site_loc_file).readlines()
header=site_loc_lines.pop(0)
names=[]
lats=[]
lons=[]
for line in site_loc_lines:
    split=line.split(',')
    names.append(split[0])
    lats.append(float(split[1]))
    lons.append(float(split[2]))

outf=open(site_loc_file+'.tmp','w')
outf.write(header)
for name,lat,lon in zip(names,lats,lons):
    distance= ( ( (lats_2D-lat)**2) + ( (lons_2D-lon)**2) ) ** 0.5
    y_index,x_index=np.where(distance==np.min(distance))
    print(name,lat,lats_2D[y_index,x_index],\
            lon,lons_2D[y_index,x_index], \
            x_index,y_index)
    outf.write('%15s, %9.5f, %9.5f, %9i, %9i \n' % \
               (name,lat,lon,x_index[0],y_index[0]) )

outf.close()

