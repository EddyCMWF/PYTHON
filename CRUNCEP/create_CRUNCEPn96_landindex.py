#!/usr/bin/env python
##################################################################################
#
# Program: fix_n96_land_frac.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: create index file for converting JULES output back to 2D
# 
##################################################################################
import numpy as np
import netCDF4 as nc
#import pylab as plt

infile_2D='/users/eow/edwcom/CRUNCEP/n96/ancil/n96-lat-lon.nc'
jules_file='/users/eow/edwcom/CRUNCEP/n96/JULES_output/JULES_v4.3_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'

output_file='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

# Read in original 2D file
inf=nc.Dataset(infile_2D,'r')
lats_2D=inf.variables['lat2d'][:]
lons_2D=inf.variables['lon2d'][:]
lats_1D=inf.variables['latitude'][:]
lons_1D=inf.variables['longitude'][:]
inf.close()

#lons_2D[np.where(lons_2D > 180.)]=lons_2D[np.where(lons_2D > 180.)]-360.

# Read in JULES output file
inf=nc.Dataset(jules_file,'r')
lats_jules=inf.variables['latitude'][:].squeeze()
lons_jules=inf.variables['longitude'][:].squeeze()
inf.close()

latlon_2D=zip(lats_2D.flat,lons_2D.flat)
latlon_jules=zip(lats_jules,lons_jules)

index_1D = np.array([ np.where((lats_2D.flat==lat_j)&(lons_2D.flat==lon_j))[0] for lat_j,lon_j in latlon_jules ]).squeeze()

index_2D = np.zeros_like(lats_2D)-999.0
index_2D.flat[index_1D] = np.arange(len(index_1D))

outf=nc.Dataset(output_file,'w')
outf.createDimension('land',len(index_1D))
outf.createDimension('lat',len(lats_1D))
outf.createDimension('lon',len(lons_1D))

outvar=outf.createVariable( 'lats_1D','float32',('land') )
outvar[:]=lats_jules
outvar=outf.createVariable( 'lons_1D','float32',('land') )
outvar[:]=lons_jules
outvar=outf.createVariable( 'index_1D','int32',('land'),fill_value=-999.)
outvar[:]=index_1D

outvar=outf.createVariable( 'lats_2D','float32',('lat','lon') )
outvar[:]=lats_2D
outvar=outf.createVariable( 'lons_2D','float32',('lat','lon') )
outvar[:]=lons_2D
outvar=outf.createVariable( 'index_2D','int32',('lat','lon'),fill_value=-999.)
outvar[:]=index_2D

outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.note='index for JULES runs on JASMIN using the CRUNCEP n96 data'
outf.note2='Python indexing, i.e. starts from 0'
outf.close()











