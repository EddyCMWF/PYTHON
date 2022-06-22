#!/usr/bin/python
#
# Python module to extract a segment of 50 km EMEP grid (in ascii)
# and store as netcdf grid file for SCRIP. 
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# September 2013
#
# Contains
#
import os, sys
import numpy as np
import
#
import netCDF4 as nc
#


################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class extract:
        
	def parse_input(self):
                
		parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
                
		# optional
		parser.add_argument('--inmissing',type=float,help='Missing value in input file',required=False,default=None)
		parser.add_argument('--outmissing',type=float,help='Missing value in output file',required=False,default=None)
                parser.add_argument('--newlatres',type=float,help='New Resolution for Latitudes',required=False,default=None)
		parser.add_argument('--newlonres',type=float,help='New Resolution for Longitudes',required=False,default=None)
                
		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('latmin',help='Minimum latitude',type=float)
		parser.add_argument('latmax',help='Maximum latitude',type=float)
		parser.add_argument('lonmin',help='Minimum longitude',type=float)
		parser.add_argument('lonmax',help='Maximum longitude',type=float)
                
		# Parse the arguments
		args=parser.parse_args()
                
		return args.infile, \
				args.outfile, \
				args.latmin, \
				args.latmax, \
				args.lonmin, \
				args.lonmax, \
				args.inmissing, \
				args.outmissing, \
                                args.newlatres, \
                                args.newlonres



infile  = '/users/eow/edwcom/EMEP/EMEPgrid.csv'
LSM_file= '/users/eow/edwcom/EMEP/Unified_EMEP_model.OpenSource2011/input/landsea_mask.dat'
outfile = '/users/eow/edwcom/EMEP/EMEPgrid.nc'

data = np.loadtxt(infile, delimiter=',')

i=data[:,0]
j=data[:,1]

lon_cen=data[:,2]
lat_cen=data[:,3]

lon_cors=data[:,[4,6,8,10]]
lat_cors=data[:,[5,7,9,11]]

LSM_data = np.loadtxt(LSM_file)
LSM = LSM_data[:,2]
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


