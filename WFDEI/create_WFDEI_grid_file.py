#!/usr/bin/python
#
# Python 
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# April 2015
#
# Contains
#
import os, sys
import numpy as np
#
import netCDF4 as nc
#

infile='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
outfile='/users/eow/edwcom/WFD_EI/wfdei-gridfile.nc'

# Read data
inf       = nc.Dataset(infile,'r')
lats      = inf.variables['latitude'][:]
lons      = inf.variables['longitude'][:]
land_frac = inf.variables['land_fraction'][:]

inf.close()


longitudes,latitudes = np.meshgrid(lons,lats)

longitudes=longitudes.flatten()
latitudes=latitudes.flatten()

lsmask=land_frac.flatten()
lsmask[lsmask>0.5]=1

rx=0.5    # x resolution
rx2=rx/2. # half of above
                                                                                             
ry=0.5    # y resolution
ry2=ry/2. # half of above 

loncorn=np.array([longitudes-rx2,longitudes+rx2,longitudes+rx2,longitudes-rx2],dtype='float32')
loncorn=loncorn.T
latcorn=np.array([latitudes-ry2,latitudes-ry2,latitudes+ry2,latitudes+rx2],dtype='float32')
latcorn=latcorn.T

grid_dims=land_frac.shape

# Open file ti write to
outf=nc.Dataset(outfile,'w')

# create dimensions
outf.createDimension('grid_size',len(latitudes))
outf.createDimension('grid_corners',4)
outf.createDimension('grid_rank',2)

# create variables and store data
outvar=outf.createVariable('grid_dims','int',('grid_rank'))
outvar[:]=grid_dims

outvar=outf.createVariable('grid_center_lat','float32',('grid_size'))
outvar.units='degrees'
outvar[:]=latitudes

outvar=outf.createVariable('grid_center_lon','float32',('grid_size'))
outvar.units='degrees'
outvar[:]=longitudes

outvar=outf.createVariable('grid_imask','int',('grid_size'))
outvar.units='unitless'
outvar[:]=lsmask

outvar=outf.createVariable('grid_corner_lat','float32',('grid_size','grid_corners'))
outvar.units='degrees'
outvar[:]=latcorn

outvar=outf.createVariable('grid_corner_lon','float32',('grid_size','grid_corners'))
outvar.units='degrees'
outvar[:]=loncorn

outf.title="WFD 0.5 degree grid"
outf.conventions="SCRIP"


outf.close()









