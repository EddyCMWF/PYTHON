#!/usr/bin/python

import netCDF4 as nc
import numpy as np

CRUNCEP_gfile = '/users/eow/edwcom/CRUNCEP/cru_ncep_land.nc'
WFDEI_gfile   = '/users/eow/edwcom/WFD_EI/EI-Halfdeg-land-elevation.nc'

out_index_file = '/users/eow/edwcom/CRUNCEP/WFDEI_to_CRUNCEP_index.dat'

Cinf  = nc.Dataset(CRUNCEP_gfile,'r')
Clats = Cinf.variables['Latitude'][:]
Clons = Cinf.variables['Longitude'][:]
Cinf.close()
Clatlon = zip(Clats,Clons) ## [ (clat,clon) for clat,clon in zip(Clats,Clons) ]

Winf=nc.Dataset(WFDEI_gfile,'r')
Wlats = Winf.variables['latitude'][:]
Wlons = Winf.variables['longitude'][:]
Winf.close()
Wlatlon = zip(Wlats,Wlons) ## [ (wlat,wlon) for wlat,wlon in zip(Wlats,Wlons) ]

#outf = open(out_index_file,'w')
#index= 0#
#
#for wlatlon_pt in Wlatlon:
#    if (wlatlon_pt in Clatlon):
#        outf.write(str(index) + '\n')
#    
#    index+=1#
#
#outf.close()

outf = open(out_index_file,'w')
for clatlon_pt in Clatlon:
    
    try:
        index = Wlatlon.index(clatlon_pt)
        outf.write(str(index) + '\n')
        success=1
    except:
        outf.write('-9999\n')
        print clatlon_pt
        success=0
    
outf.close()




