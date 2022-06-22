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

infile_2D='/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'
jules_file='/users/global/albmar/CHESS/outputs/chess_4.1.month.nc'

output_file='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

# Read in original 2D file
inf=nc.Dataset(infile_2D,'r')
lats_2D=inf.variables['lat'][:]
lons_2D=inf.variables['lon'][:]
inf.close()

#lons_2D[np.where(lons_2D > 180.)]=lons_2D[np.where(lons_2D > 180.)]-360.

# Read in JULES output file
inf=nc.Dataset(jules_file,'r')
lats_jules=inf.variables['latitude'][:].squeeze()
lons_jules=inf.variables['longitude'][:].squeeze()
inf.close()

latlon_2D=zip(lats_2D.flat,lons_2D.flat)
latlon_jules=zip(lats_jules, lons_jules)

#index = [ np.where( (lats_2D==lat_j)&

index_1D = np.array([ np.where((lats_2D.flat==lat_j)&(lons_2D.flat==lon_j))[0] \
                          for lat_j,lon_j in latlon_jules ]).squeeze()

index_2D = np.zeros_like(lats_2D)-999.0
index_2D.flat[index_1D] = np.arange(len(index_1D))

index_to1D = np.array( [ np.where( (lats_2D==lat_j)&(lons_2D==lon_j) ) \
                        lat_j, lon_j in latlon_jules ] ).squeeze()  
index_to_1D=index_to_1D.transpose(1,0)

outf=nc.Dataset(output_file,'w')
outf.createDimension('lat',lats_2D.shape[0])
outf.createDimension('lon',lats_2D.shape[1])
outf.createDimension('landpoints',lats_jules.shape[0])
outf.createDimension('dim',2)


outvar=outf.createVariable( 'lats_2D','float32',('lat','lon') )
outvar[:]=lats_2D
outvar=outf.createVariable( 'lons_2D','float32',('lat','lon') )
outvar[:]=lons_2D
outvar=outf.createVariable( 'index_2D','int32',('lat','lon'),fill_value=-999.)
outvar.longname='index to go from jules landpoints to 2D grid'
outvar[:]=index_2D


outvar=outf.createVariable( 'lats_1D','float32',('landpoints') )
outvar[:]=lats_jules
outvar=outf.createVariable( 'lons_1D','float32',('landpoints') )
outvar[:]=lons_jules
outvar=outf.createVariable( 'index_to1D','int32',('dim','landpoints') )
outvar.longname='index to go from 2D grid to jules landpoints'
outvar[:]=index_to1D

outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.note='CHESS land index for JULES output'
outf.note2='Python indexing, i.e. starts from 0'
outf.close()











