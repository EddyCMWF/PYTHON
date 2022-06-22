#!/bin/env python

import netCDF4 as nc
import numpy as np

Cdir='/users/eow/edwcom/CHESS/'

latlon_file = Cdir+'chess_jules_land_index.nc'
inf=nc.Dataset(latlon_file,'r')
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
lats_1D=inf.variables['lats_1D'][:]
lons_1D=inf.variables['lons_1D'][:]
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
    distance=( ((lats_2D-lat)**2) + ((lons_2D-lon)**2) ) **0.5
    y_index,x_index=np.where(distance==np.min(distance))
    distance_1D=((lats_1D-lat)**2) + ((lons_1D-lon)**2)
    lp_index=np.where(distance_1D==np.min(distance_1D))[0]
    print('%15s, %9.5f, %9.5f, %9i, %9i, %9i \n' % \
          (name,lat,lon,x_index[0],y_index[0],lp_index[0]) )

    outf.write('%15s, %9.5f, %9.5f, %9i, %9i, %9i \n' % \
               (name,lat,lon,x_index[0],y_index[0],lp_index[0]) )

outf.close()

