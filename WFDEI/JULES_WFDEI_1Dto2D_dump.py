#!/usr/bin/python2.7

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
#import argparse
#import getpass
#import datetime as dt

out_mv=-9999.

WFDEI_DIR='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/'

WFD_1d_infile = WFDEI_DIR+'dumps/WFDEI_newtopo_startdump_withSC_ECP20150731.nc'
WFD_gridfile  = WFDEI_DIR+'ancils/wfdei-land-mask.nc'

WFD_2d_outfile= WFDEI_DIR+'/dumps/WFDEI_newtopo_startdump_withSC_ECP20150731_2D.nc'

# read gridfile
inf_grid = nc.Dataset(WFD_gridfile,'r')

out_dim_lons=len(inf_grid.dimensions['longitude'])
out_dim_lats=len(inf_grid.dimensions['latitude'])

lons     = inf_grid.variables['longitude'][:]
lats     = inf_grid.variables['latitude'][:]
LAND_FRAC= inf_grid.variables['land_fraction'][:]
LAND_IND = inf_grid.variables['land_index'][:]

inf_grid.close()

index=LAND_IND.data[LAND_IND.mask==False]-1

out_dim_land=out_dim_lons*out_dim_lats

# open 1D file to read data
inf=nc.Dataset(WFD_1d_infile,'r')

#read dimensions
in_dim_land   = len(inf.dimensions['land'])
in_dim_tile   = len(inf.dimensions['tile'])
in_dim_scpool = len(inf.dimensions['scpool'])
in_dim_soil   = len(inf.dimensions['soil'])
in_dim_snow   = len(inf.dimensions['snow'])

# open 2D file to write data
outf=nc.Dataset(WFD_2d_outfile,'w')
# Write dimensions to file
outf.createDimension('lon',out_dim_lons)
outf.createDimension('lat',out_dim_lats)
outf.createDimension('tile',in_dim_tile)
outf.createDimension('scpool',in_dim_scpool)
outf.createDimension('soil',in_dim_soil)
outf.createDimension('snow',in_dim_snow)


# read in canopy data
var='canopy'
in_dat=inf.variables[var][:]
#convert to 2D
out_dat=in_dat[:,LAND_IND-1]
#write to file
outfoutvar=outf.createVariable(var,'float32',('tile','lat','lon'))
outfoutvar[:]=out_dat

var='cs'
in_dat=inf.variables[var][:]
out_dat=in_dat[:,LAND_IND-1]
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
