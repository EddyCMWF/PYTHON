#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import getpass
import datetime as dt
import mpl_toolkits.basemap.pyproj as pyproj 


################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class gridc:

################################################################################
# Useful stuff
################################################################################
	def __init__(self):
		# projection from:
		# http://spatialreference.org/ref/epsg/27700/proj4/
		self.osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 
		self.osgb36_mine=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=455000 +y_0=-88000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

		# http://spatialreference.org/ref/epsg/4258/
		self.etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')




################################################################################
# Conversion stuff
################################################################################

	def bng_to_ll(self,x,y):
		newx, newy = pyproj.transform(gridc().osgb36, gridc().etrs89, x, y )

		return newx, newy


################################################################################
# Get the corners
################################################################################
	def corners(self,lats,lons,resx,resy,redcy=False):

		rx2=resx/2.0
		ry2=resy/2.0
		
		loncorn = np.array([ [lon-rx2,lon+rx2,lon+rx2,lon-rx2] for lon in lons ])
		latcorn = np.array([ [lat-ry2,lat-ry2,lat+ry2,lat+ry2] for lat in lats ])

		# Make the cells that intersect the poles into triangles, with
		# a repeated point for redundancy
		if redcy:
			for i in np.where(lats==90.0)[0]:
				# For the north pole, we have a single point on 
				# the pole at (lon, lat+r2) and this is 
				# duplicated as it's the last in the array
				loncorn[i,2:]=lons[i]

			for i in np.where(lats==-90.0)[0]:
				# For the south pole, we have a single point on
				# the pole at (lon, lat-r2), but we duplicate
				# the last corner point at (lon-r2,lat+r2
				# These changes to the array make that happen,
				loncorn[i,0]=lons[i]
				loncorn[i,2]=lons[i]-r2
				latcorn[i,1]=lats[i]+r2


		return latcorn, loncorn

################################################################################
# Parse input
################################################################################

	def parse_input(self):

		parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')

		# optional

		# positional
		parser.add_argument('landfile',help='Landmask file')
		parser.add_argument('outfile',help='Output file')

		# Parse the arguments
		args=parser.parse_args()

		return args.landfile, \
				args.outfile

################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################


if __name__=='__main__':

	# Call the class
	g=gridc()

	landfile, outfile = g.parse_input()

	# read the land mask
	lf=nc.Dataset(landfile,'r')
	landfrac=lf.variables['landfrac'][:]

	# get the coordinates
	# The input file now contains grid centres, not lower lefts
	ct_east = lf.variables['x'][:]
	ct_north = lf.variables['y'][:]
	nx=len(ct_east)
	ny=len(ct_north)
	ct_east_g,ct_north_g=np.meshgrid(ct_east,ct_north)
	ct_east_v=ct_east_g.flatten()
	ct_north_v=ct_north_g.flatten()
	nland=len(ct_east_v)

	# Get the 'mask'
	gmask=np.zeros(nland)
	gmask[np.where(landfrac.flatten()>0)]=1.0

	# Get the full grid centres and corners arrays
	grid_centres={}
	grid_centres['lat']=np.zeros(nland)
	grid_centres['lon']=np.zeros(nland)
	grid_centres['east']=np.zeros(nland)
	grid_centres['north']=np.zeros(nland)

	grid_corners={}
	grid_corners['lat']=np.zeros(nland*4)
	grid_corners['lon']=np.zeros(nland*4)
	grid_corners['east']=np.zeros(nland*4)
	grid_corners['north']=np.zeros(nland*4)

	# Centres were read
	grid_centres['east']=ct_east_v
	grid_centres['north']=ct_north_v

	# Squars are 1km, so centres are 500m away from lower left corner
	ll_east_v=ct_east_v-500
	ll_north_v=ct_north_v-500

	# Corners are anti-clockwise from lower left
	grid_corners['east'][0::4]=ll_east_v
	grid_corners['east'][1::4]=ll_east_v+1000
	grid_corners['east'][2::4]=ll_east_v+1000
	grid_corners['east'][3::4]=ll_east_v

	grid_corners['north'][0::4]=ll_north_v
	grid_corners['north'][1::4]=ll_north_v
	grid_corners['north'][2::4]=ll_north_v+1000
	grid_corners['north'][3::4]=ll_north_v+1000

	# Convert coords to lat/lon
	for i in range(nland):
		grid_centres['lon'][i],grid_centres['lat'][i] = g.bng_to_ll(grid_centres['east'][i],grid_centres['north'][i])

	for i in range(nland*4):
		grid_corners['lon'][i],grid_corners['lat'][i] = g.bng_to_ll(grid_corners['east'][i],grid_corners['north'][i])


	print "Output file: %s"%outfile
	outf=nc.Dataset(outfile,'w')

	# create dimensions
	outf.createDimension('grid_size',len(gmask))
	outf.createDimension('grid_corners',4)
	outf.createDimension('grid_rank',2)

	# create variables
	gdims=outf.createVariable('grid_dims','i',('grid_rank',))
	#gdims[:]=[nland,1]
	gdims[:]=[nx,ny]


	gctr_lat=outf.createVariable('grid_center_lat','f',('grid_size',))
	gctr_lat.units='degrees'
	gctr_lat[:]=grid_centres['lat'][:]

	gctr_lon=outf.createVariable('grid_center_lon','f',('grid_size',))
	gctr_lon.units='degrees'
	gctr_lon[:]=grid_centres['lon'][:]

	gim=outf.createVariable('grid_imask','i',('grid_size',))
	gim[:]=gmask.astype('int')
	gim.units='unitless'

	gcrn_lat=outf.createVariable('grid_corner_lat','f',('grid_size','grid_corners'))
	gcrn_lat.units='degrees'
	gcrn_lat[:]=grid_corners['lat'][:]

	gcrn_lon=outf.createVariable('grid_corner_lon','f',('grid_size','grid_corners'))
	gcrn_lon.units='degrees'
	gcrn_lon[:]=grid_corners['lon'][:]

	outf.title="CHESS 1km grid"

	outf.conventions='SCRIP'

	outf.close()






