#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import getpass
import datetime as dt

################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class gridc:
################################################################################
# Flatten grid
################################################################################
	def grid_to_vect(self, grid):

		vect=grid.flatten()

		return vect

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
		parser.add_argument('infile',help='Input file')
		parser.add_argument('outfile',help='Output file')

		# Parse the arguments
		args=parser.parse_args()

		return args.infile, \
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

	infile, outfile = g.parse_input()

	# Open input file
	print "Read input file: %s"%infile
	inf=nc.Dataset(infile,'r')

	# Set up grid
	# These will be the CENTRES of the grid squares
	lon=inf.variables['Longitude'][:]
	lat=inf.variables['Latitude'][:]

	nx=len(lon)
	ny=len(lat)

	glon,glat=np.meshgrid(lon,lat)

	glat=g.grid_to_vect(glat)
	glon=g.grid_to_vect(glon)

	indata=g.grid_to_vect(inf.variables['Z'][:])
	gmask=np.ones_like(glon)
	gmask[np.where(indata.mask)]=0.0

	# get corners
	resx=np.abs(lon[1]-lon[0])
	resy=np.abs(lat[1]-lat[0])
	glatcorn,gloncorn=g.corners(glat,glon,resx,resy)
	
	# Open output file
	print "opening file: "+outfile
	outf=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')

	# create dimensions
	outf.createDimension('grid_size',len(gmask))
	outf.createDimension('grid_corners',4)
	outf.createDimension('grid_rank',2)

	# create variables
	gdims=outf.createVariable('grid_dims','i',('grid_rank',))
	gdims[:]=[nx,ny]


	gctr_lat=outf.createVariable('grid_center_lat','double',('grid_size',))
	gctr_lat.units='degrees'
	gctr_lat[:]=glat

	gctr_lon=outf.createVariable('grid_center_lon','double',('grid_size',))
	gctr_lon.units='degrees'
	gctr_lon[:]=glon

	gim=outf.createVariable('grid_imask','i',('grid_size',))
	gim[:]=gmask.astype('int')
	gim.units='unitless'

	gcrn_lat=outf.createVariable('grid_corner_lat','double',('grid_size','grid_corners'))
	gcrn_lat.units='degrees'
	gcrn_lat[:]=glatcorn

	gcrn_lon=outf.createVariable('grid_corner_lon','double',('grid_size','grid_corners'))
	gcrn_lon.units='degrees'
	gcrn_lon[:]=gloncorn

	outf.title="%.9f x %.9f degree"%(resx,resy)

	outf.conventions='SCRIP'

	outf.close()











	

