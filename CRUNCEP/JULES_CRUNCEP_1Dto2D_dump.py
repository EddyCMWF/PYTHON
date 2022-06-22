#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
#import argparse
#import getpass
#import datetime as dt

out_mv=-9999.

CRU_1d_infile  = '/users/eow/edwcom/CRUNCEP/start_dumps/CRU_NCEP_v4_global.dump.20010101.0.nc'
CRU_indexfile   = '/users/eow/edwcom/CRUNCEP/cru_ncep_land.nc'
CRU_LandFracfile   = '/users/eow/edwcom/CRUNCEP/cru_ncep_LandFraction.nc'

CRU_2d_outfile = '/users/eow/edwcom/CRUNCEP/start_dumps/CRU_NCEP_v4_global_2D.dump.20010101.0.nc'
CRU_2d_landmask= '/users/eow/edwcom/CRUNCEP/cruncep-land-mask.nc'

#out_dim_lons=len(inf_grid.dimensions['longitude'])
#out_dim_lats=len(inf_grid.dimensions['latitude'])

lons     = np.arange(-180,180,0.5)+0.25
lats     = np.arange(-90,90,0.5)+0.25


lon_grid,lat_grid=np.meshgrid(lons,lats)

# read gridfile
inf_index = nc.Dataset(CRU_indexfile,'r')
index = inf_index.variables['land'][:]
inf_index.close()

# read landfrac file
inf_landfrac = nc.Dataset(CRU_LandFracfile,'r')
lsmask=inf_landfrac.variables['lsmask'][:]
inf_landfrac.close()

LAND_IND=np.zeros_like(lon_grid,dtype='int')-9999
LAND_IND.flat[index]=np.arange(len(index))
LAND_IND=LAND_IND[::-1,:]
LAND_IND=np.ma.masked_equal(LAND_IND,-9999.0)

LAND_FRAC=np.zeros_like(lon_grid)-9999.0
LAND_FRAC.flat[index]=lsmask
LAND_FRAC=LAND_FRAC[::-1,:]
LAND_FRAC=np.ma.masked_equal(LAND_FRAC,-9999.0)



# Output index and land frac to land mask file
# Same conventions as Phils wfd-ei equivalent
outf_LM=nc.Dataset(CRU_2d_landmask,'w')
outf_LM.createDimension('longitude',len(lons))
outf_LM.createDimension('latitude',len(lats))
outvar=outf_LM.createVariable('longitude','float32',('longitude'))
outvar.units='degrees_east'
outvar[:]=lons
outvar=outf_LM.createVariable('latitude','float32',('latitude'))
outvar.units='degrees_north'
outvar[:]=lats
outvar=outf_LM.createVariable('land_fraction','float32',('latitude','longitude'),fill_value=-9999.0)
outvar.units='unitless'
outvar[:]=LAND_FRAC
outvar=outf_LM.createVariable('land_index','int',('latitude','longitude'),fill_value=-9999)
outvar.units='unitless'
outvar[:]=LAND_IND+1
outf_LM.note='Land-sea mask used with CRU-NCEP forcing data.'
outf_LM.indexing='Fortran'
outf_LM.author='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
outf_LM.close()

#plt.subplot(2,2,1)
#plt.imshow(LAND_IND,vmin=0,origin='bottom')
#plt.colorbar()
#plt.subplot(2,2,2)
#plt.imshow(LAND_FRAC,vmin=0,origin='bottom')
#plt.colorbar()
#plt.subplot(2,2,3)
#plt.imshow(lat_grid,origin='bottom')
#plt.colorbar()
#plt.subplot(2,2,4)
#plt.imshow(lon_grid,origin='bottom')
#plt.colorbar()
#plt.show()


# open 1D file to read data
inf=nc.Dataset(CRU_1d_infile,'r')

#read dimensions
in_dim_land   = len(inf.dimensions['land'])
in_dim_tile   = len(inf.dimensions['tile'])
in_dim_scpool = len(inf.dimensions['scpool'])
in_dim_soil   = len(inf.dimensions['soil'])
in_dim_snow   = len(inf.dimensions['snow'])

# open 2D file to write data
outf=nc.Dataset(CRU_2d_outfile,'w')
# Write dimensions to file
outf.createDimension('lon',len(lons))
outf.createDimension('lat',len(lats))
outf.createDimension('tile',in_dim_tile)
outf.createDimension('scpool',in_dim_scpool)
outf.createDimension('soil',in_dim_soil)
outf.createDimension('snow',in_dim_snow)


# read in canopy data
var='canopy'
in_dat=inf.variables[var][:]
#convert to 2D
out_dat=in_dat[:,LAND_IND]
#write to file
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat

var='cs'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND]
outfoutvar=outf.createVariable(var,'float32',('scpool','lat','lon'))
outfoutvar[:]=out_dat

var='gs'
in_dat=inf.variables[var][:]
out_dat=in_dat[LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('lat','lon'))
outfoutvar[:]=out_dat


var='snow_tile'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='sthuf'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('soil','lat','lon'))
outfoutvar[:]=out_dat


var='t_soil'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('soil','lat','lon'))
outfoutvar[:]=out_dat


var='tstar_tile'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='sthzw'
in_dat=inf.variables[var][:]
out_dat=in_dat[LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('lat','lon'))
outfoutvar[:]=out_dat


var='zw'
in_dat=inf.variables[var][:]
out_dat=in_dat[LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('lat','lon'))
outfoutvar[:]=out_dat


var='rho_snow'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='snow_depth'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='snow_grnd'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='nsnow'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat


var='snow_ds'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('snow','tile','lat','lon'))
outfoutvar[:]=out_dat


var='snow_ice'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('snow','tile','lat','lon'))
outfoutvar[:]=out_dat


var='snow_liq'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('snow','tile','lat','lon'))
outfoutvar[:]=out_dat


var='tsnow'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,:,LAND_IND-1]
outfoutvar=outf.createVariable(var,'float32',('snow','tile','lat','lon'))
outfoutvar[:]=out_dat


inf.close()
outf.close()
