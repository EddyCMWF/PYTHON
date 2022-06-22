#!/usr/bin/python

import gdal
import numpy as np
import netCDF4 as nc

topo_infile='/users/eow/edwcom/WFD_EI/topoidx_WFDEI_0p5_lp_global.nc'

topo_outfile='/users/eow/edwcom/WFD_EI/topoidx_WFDEI_0p5_lp_global_filled.nc'


# Read in lat and lons from CRUNCEP land file
inf = nc.Dataset(topo_infile,'r')
land = inf.variables['land'][:]
lats = inf.variables['latitude'][:]
lons = inf.variables['longitude'][:]
ti_mean=inf.variables['ti_mean'][:]
ti_std=inf.variables['ti_std'][:]
fexp=inf.variables['fexp'][:]
inf.close()

badex = np.where(ti_mean.mask==True)[0]
goodex= np.where(ti_mean.mask==False)[0]

for pt in badex:
    temppoint = np.argmin(  ((lats[pt]-lats[goodex])**2) \
                          + ((lons[pt]-lons[goodex])**2) )
    new_point=goodex[temppoint]
    ti_mean[pt]=ti_mean[new_point]
    ti_std[pt]=ti_std[new_point]


badex = np.where(ti_std.mask==True)[0]
goodex= np.where(ti_std.mask==False)[0]
for pt in badex:
    temppoint = np.argmin(  ((lats[pt]-lats[goodex])**2) \
                          + ((lons[pt]-lons[goodex])**2) )
    new_point=goodex[temppoint]
    ti_mean[pt]=ti_mean[new_point]
    ti_std[pt]=ti_std[new_point]

outf=nc.Dataset(topo_outfile,'w')

outf.createDimension('land',len(lats))

#land index
outvar       = outf.createVariable('land','int',('land'))
outvar.units = "-"
outvar.note  = "Index for grid, row-wise from NW corner" 
outvar[:]    = land

#longitude
outvar                = outf.createVariable('longitude','float32',('land'))
outvar.units          = "degrees east"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre longitude"
outvar[:]             = lons

#latitude
outvar                = outf.createVariable('latitude','float32',('land'))
outvar.units          = "degrees north"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre latitude"
outvar[:]             = lats

#ti_mean
outvar                = outf.createVariable('ti_mean','float32',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index mean"
outvar[:]             = ti_mean

#ti_std
outvar                = outf.createVariable('ti_std','float32',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index standard deviation"
outvar[:]             = ti_std

#fexp
outvar                = outf.createVariable('fexp','float32',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic f exponent"
outvar[:]             = fexp

outf.setncattr('title','Aggregate topographic index data for WFD-EI')
outf.setncattr('institution','CEH - Wallingford')
outf.setncattr('source','Simon Dadson and Toby Marthews, U. Oxford')
outf.setncattr('contact', 'E. Comyn-Platt (edwcom@ceh.ac.uk)')
outf.setncattr('note','0.5 degree data provided by T. Marthews with Ice points filled with closest neighbour')

outf.close()



