#!/usr/bin/env python

################################################################################
#
# Script to create a LAI phenology file for the new CHESS runs
#
# ELR 04/03/2014
#
################################################################################
#
# Copied from make_new_phenol_from_old.py and edited to
#   a) create a monthly climatology instead of daily
#   b) use a 2d mask instead of doing the landvector
#
# ELR 14/05/2014
#
################################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import calendar as cal
import mpl_toolkits.basemap.pyproj as pyproj 
import datetime as dt


################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class lai:

################################################################################
# Useful stuff
################################################################################
	def __init__(self):
		# projection from:
		# http://spatialreference.org/ref/epsg/27700/proj4/
		self.osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

		# http://spatialreference.org/ref/epsg/4258/
		self.etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')


################################################################################
# Conversion stuff
################################################################################

	def bng_to_ll(self,x,y):
		lonlat = [ pyproj.transform(self.osgb36, self.etrs89, x[i], y[i] ) for i in range(len(x)) ]

		lon = np.array([l[0] for l in lonlat])
		lat = np.array([l[1] for l in lonlat])

		return lon, lat

################################################################################
# Parse input
################################################################################

	def parse_input(self):

		parser=argparse.ArgumentParser(description='Create LAI file')

		# optional 
		parser.add_argument('-m','--outmissing',help='Missing value for output',required=False,default=-1e20,type=float)

		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('maskfile',help='Mask file')
		parser.add_argument('maskvar',help='Mask var')
		parser.add_argument('outfile',help='Output file')

		# Parse the arguments
		args=parser.parse_args()

		return args.infile, \
				args.maskfile, \
				args.maskvar, \
				args.outfile, \
				args.outmissing

################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################


if __name__=='__main__':

	# Call the class
	l=lai()

	infile, maskfile, maskvar, outfile, outmv = l.parse_input()

	nt=366
	nz=5

	inf=np.fromfile(infile,dtype='float32')
	inf=inf.byteswap().reshape((nt,nz))

	ntout=12
	mnend=0
	clim=[]
	for t in range(ntout):
		mnstart=mnend
		mnend=mnstart+cal.monthrange(1964,t+1)[1]
		clim.append(inf[mnstart:mnend].mean(axis=0))



	mf=nc.Dataset(maskfile,'r')
	maskdata=mf.variables[maskvar][:]
	maskindx=np.where(maskdata>0)
	eastings=mf.variables['x'][:]
	northings=mf.variables['y'][:]

	nx=len(eastings)
	ny=len(northings)


	outf=nc.Dataset(outfile,'w')

	outf.createDimension('time',ntout)
	outf.createVariable('time','f',('time',),fill_value=outmv)
	outf.variables['time'].units = "days since 1961-01-01"
	outf.variables['time'].calendar = "gregorian" ;
	outf.variables['time'].long_name = "Time in days days since 1961-01-01"
	outf.variables['time'][:] = np.array([0,]+[ cal.monthrange(1961,mn)[1] for mn in range(1,12)]).cumsum()


	outf.createDimension('pft',nz)
	outf.createVariable('pft','f',('pft',),fill_value=outmv)
	outf.variables['pft'].units=''
	outf.variables['pft'].long_name = "plant functional type"
	outf.variables['pft'].standard_name = "plant_functional_type" 
	outf.variables['pft'].axis = 'z'
	outf.variables['pft'][:]=np.arange(1,nz+1)

	outf.createDimension('x',nx)
	outf.createVariable('x','f',('x',),fill_value=outmv)
	outf.variables['x'].units='m'
	outf.variables['x'].long_name = "easting - OSGB36 grid reference" 
	outf.variables['x'].standard_name = "projection_x_coordinate" 
	outf.variables['x'].units='m'
	outf.variables['x'].point_spacing = "even" 
	outf.variables['x'][:]=np.arange(0,nx*1000,1000)+500.0

	outf.createDimension('y',ny)
	outf.createVariable('y','f',('y',),fill_value=outmv)
	outf.variables['y'].units='m'
	outf.variables['y'].long_name = "northing - OSGB36 grid reference" 
	outf.variables['y'].standard_name = "projection_y_coordinate" 
	outf.variables['y'].units='m'
	outf.variables['y'].point_spacing = "even" 
	outf.variables['y'][:]=np.arange(0,ny*1000,1000)+500.0

	outf.createVariable('lat','f',('y','x'),fill_value=outmv)
	outf.variables['lat'].long_name = "latitude of grid box centre" 
	outf.variables['lat'].standard_name = "latitude" 
	outf.variables['lat'].units = "degrees_north" 

	outf.createVariable('lon','f',('y','x'),fill_value=outmv)
	outf.variables['lon'].long_name = "longitude of grid box centre" 
	outf.variables['lon'].standard_name = "longitude" 
	outf.variables['lon'].units = "degrees_east" 

	xg,yg=np.meshgrid(np.arange(0,nx*1000,1000),np.arange(0,ny*1000,1000))
	xg+=500
	yg+=500

	glon,glat = l.bng_to_ll(xg.flatten(),yg.flatten())
	glon=glon.reshape([ny,nx])
	glat=glat.reshape([ny,nx])

	outf.variables['lat'][:]=glat[:]
	outf.variables['lon'][:]=glon[:]

	# Coord variable
	outf.createVariable('crs','i',[])
	outf.variables['crs'].long_name = "coordinate_reference_system" 
	outf.variables['crs'].grid_mapping_name = "transverse_mercator" 
	outf.variables['crs'].EPSG_code = "EPSG:27700" 
	outf.variables['crs'].semi_major_axis = 6377563.396 
	outf.variables['crs'].semi_minor_axis = 6356256.91 
	outf.variables['crs'].inverse_flattening = 299.3249646 
	outf.variables['crs'].latitude_of_projection_origin = 49. 
	outf.variables['crs'].longitude_of_projection_origin = -2. 
	outf.variables['crs'].false_easting = 400000. 
	outf.variables['crs'].false_northing = -100000. 
	outf.variables['crs'].scale_factor_at_projection_origin = 0.9996012717 

	# data variable
	outf.createVariable('lai','f',('time','pft','y','x'),fill_value=outmv,zlib=True,complevel=6,shuffle=True)
	outf.variables['lai'].standard_name='LAI'
	outf.variables['lai'].long_name='leaf_area_index'
	outf.variables['lai'].units='none'
	outf.variables['lai'].grid_mapping='crs'
	outf.variables['lai'].coordinates='pft lat lon'


	lai=np.ones([ntout,nz,ny,nx])*outmv
	for t in range(ntout):
		for z in range(nz):
			lai[t,z,:,:][maskindx]=clim[t][z]

	lai=np.ma.masked_equal(lai,outmv)

	outf.variables['lai'][:]=lai[:]

	outf.title = "Leaf area index (LAI)" 
	outf.description = "Leaf area index (LAI) for five PFTs at 1km resolution over Great Britain" 
	outf.institution = "CEH Wallingford - NERC" 
	now=dt.datetime.now()
	outf.history = "Created "+now.strftime('%Y-%M-%d %h:%m:%s')
	outf.date_created = now.strftime('%Y-%M-%d')


	outf.geospatial_lat_min = glat.min()
	outf.geospatial_lat_max = glat.max()
	outf.geospatial_lon_min = glon.min()
	outf.geospatial_lon_max = glon.max()
	outf.version = "beta_version" 
	outf.Conventions = "CF-1.6" 
	outf.creator_name = "Emma L. Robinson" 
	outf.creator_email = "emrobi@ceh.ac.uk" 

	outf.close()



