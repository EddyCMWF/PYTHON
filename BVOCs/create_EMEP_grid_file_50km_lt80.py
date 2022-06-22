#!/usr/bin/python
#
# Python module to convert 50 km EMEP grid from ascii to ncdf
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# September 2013
#
# Contains
#
import os, sys
import numpy as np
#
import netCDF4 as nc
#


infile  = '/users/eow/edwcom/EMEP/EMEPgrid.csv'
LSM_file= '/users/eow/edwcom/EMEP/Unified_EMEP_model.OpenSource2011/input/landsea_mask.dat'
outfile = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/EMEPgrid_lt80.nc'

data = np.loadtxt(infile, delimiter=',')

data=data[data[:,3]<80.,:]
data=data[::-1,:]

i=data[:,0]
j=data[:,1]

lon_cen=data[:,2]
lat_cen=data[:,3]

lon_cors=data[:,[5,7,9,11]]
lat_cors=data[:,[6,8,10,12]]

LSM_data = np.loadtxt(LSM_file)
LSM = LSM_data[data[:,3]<80.,2]
LSM=LSM[::-1]
LSM[LSM>1]=1

NDIMS = 3

dim_name = ['grid_size','grid_corners','grid_rank']
dim_size = [len(i),4,2]
dim_type = ['int','int','int']

NC_ID = nc.Dataset(outfile,'w')
setattr(NC_ID, 'title', "EMEP 50km grid")
setattr(NC_ID, 'conventions', "SCRIP" )
for cnt in range(len(dim_name)):
    NC_ID.createDimension(dim_name[cnt],dim_size[cnt])

VAR_ID = NC_ID.createVariable('grid_dims',np.dtype('int32').char,('grid_rank'))
VAR_ID[:] = np.array([len(i),4],dtype=int)
del VAR_ID

VAR_ID = NC_ID.createVariable('grid_center_lat','float64',('grid_size'))
VAR_ID[:] = lat_cen
setattr(VAR_ID,'units','degrees')
del VAR_ID

VAR_ID = NC_ID.createVariable('grid_center_lon','float64',('grid_size'))
VAR_ID[:] = lon_cen
setattr(VAR_ID,'units','degrees')
del VAR_ID

VAR_ID = NC_ID.createVariable('grid_imask','int32',('grid_size'))
VAR_ID[:] = LSM
setattr(VAR_ID,'units','unitless')
del VAR_ID

VAR_ID = NC_ID.createVariable('grid_corner_lat','float64',('grid_size','grid_corners'))
VAR_ID[:] = lat_cors
setattr(VAR_ID,'units','degrees')
del VAR_ID

VAR_ID = NC_ID.createVariable('grid_corner_lon','float64',('grid_size','grid_corners'))
VAR_ID[:] = lon_cors
setattr(VAR_ID,'units','degrees')
del VAR_ID

NC_ID.close()


