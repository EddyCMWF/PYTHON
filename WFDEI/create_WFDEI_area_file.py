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

resolution='0.25'

if resolution=='0.5':
    grid_file='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
    frac_name='land_fraction'
    index_name='land_index'
    area_file_2D='/prj/ALANIS/UM_Modelling/EMISSIONS/GridCell_Area_0.5.nc'
    wfdei_area_file='/users/eow/edwcom/WFD_EI/wfdei-land-area.nc'
elif resolution=='0.25':
    grid_file='/users/eow/edwcom/WFD_EI/e2obs_wrr2_landfrac_025_land_only.nc'
    frac_name='land_frac'
    index_name='land_index'
    area_file_2D='/prj/ALANIS/UM_Modelling/EMISSIONS/GridCell_Area_0.25.nc'
    wfdei_area_file='/users/eow/edwcom/WFD_EI/wfdei-land-area_0.25.nc'


# Read data
grinf = nc.Dataset(grid_file,'r')
lats  = grinf.variables['latitude'][:]
lons  = grinf.variables['longitude'][:]
index = grinf.variables[index_name][:]
frac  = grinf.variables[frac_name][:]
grinf.close()

nLandPts=index.max()+1

Ainf  = nc.Dataset(area_file_2D,'r')
area_data=Ainf.variables['area'][:].squeeze()
Ainf.close()

area_data=np.ma.masked_array(area_data*frac,mask=index.mask)

area_data_1D=np.zeros(nLandPts)
area_data_1D[index[area_data.mask==False]-1]=area_data[area_data.mask==False]

if resolution in ['0.5']:
    longitudes,latitudes = np.meshgrid(lons,lats)
else:
    longitudes,latitudes=lons,lats

lats_1D=np.zeros(nLandPts)
lats_1D[index[area_data.mask==False]-1]=latitudes[area_data.mask==False]
lons_1D=np.zeros(nLandPts)
lons_1D[index[area_data.mask==False]-1]=longitudes[area_data.mask==False]


# Open file to write to
outf=nc.Dataset(wfdei_area_file,'w')

# create dimensions
outf.createDimension('land',len(lats_1D))

outvar=outf.createVariable('lat','float32',('land'))
outvar.units='degrees'
outvar[:]=lats_1D

outvar=outf.createVariable('lon','float32',('land'))
outvar.units='degrees'
outvar[:]=lons_1D

outvar=outf.createVariable('area','float32',('land'))
outvar.units='m^2'
outvar[:]=area_data_1D

outf.close()









