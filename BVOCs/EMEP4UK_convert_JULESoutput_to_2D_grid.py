#!/usr/bin/env python
##################################################################################
#
# Program: EMEP4UK_convert_JULESoutput_to_2D_grid.py
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: Module for converting 1D EME4UK data to original EMEP4UK 2D grid
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import argparse
#
#
##################################################################################
# conv_EMEP4UK_2D
# 
# function to convert to EMEP4UK 2D grid
#
# Required Inputs:
#   data - original data
#   lons - original lons
#   lats - original lats
#   
# Optional Inputs:
#   grid_file - netCDF file containing 2D grid coords 
#
def conv_EMEP4UK_2D(data,lons,lats,                                                 \
                    grid_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Landuse.nc',  \
                    grinf_latname='lat',grinf_lonname='lon',missing_value=-9999.0,  \
                    outfilename=None )
#
#
grinf=nc.Dataset(grid_file,'r')
grinf_lats=grinf.variables[grinf_latname][:]
grinf_lons=grinf.variables[grinf_lonname][:]
grinf.close()
grinf_dims=grinf_lons.shape

# create index by finding the min distance between input lat/lons and grid_file lat/lons
index = [ np.argmin( ((grinf_lats-lat)**2) + ((grinf_lons-lon)**2) ) \
          for lat,lon in zip(list(lats.flatten()),list(lons.flatten())) ]

new_data=np.zeros_like(grinf_lats)+missing_value
new_data.flat[index]=data.flatten()

if (outfilename!=None):
    outf=nc.Dataset(outfilename,'w')
    outf.createDimension('x',len(index))
    
    outvar=outf.createVariable('lat','float32',('x'))
    outvar[:]=lats
    outvar=outf.createVariable('lon','float32',('x'))
    outvar[:]=lons
    outvar=outf.createVariable('index','int',('x'))
    outvar[:]=lats
    
    outf.author='Edward Comyn-Platt; edwcom@ceh.ac.uk'
    
    outf.close()
