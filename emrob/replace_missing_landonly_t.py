#!/usr/bin/env python

################################################################################
#
# Script to replace missing values in a soil file with the nearest soil
#
# ELR 17/02/2014
#
################################################################################
#
# Copied from replace_missing.py and edited to allow a mask that tells us if
# the point is supposed to be missing (ie. that it's sea)
#
# ELR 14/05/2014
#
################################################################################
#
# Also mask out based on landmask (as it didn't do it in the remapping, grr)
#
# ELR 15/05/2014
#
################################################################################
# 
# Copied from ../soils/src/replace_missing_landonly.py and added 
# time-templating of input data. Note that we still assume that t and z
# dimensions are missing or degenerate
# Reverted to using pre-computed points for filling, rather than calculating
# on the fly
#
# ELR 22/05/2014
#
################################################################################

import netCDF4 as nc
import sys
import numpy as np
import argparse
import matplotlib
#matplotlib.use('Agg')
import mpl_toolkits.basemap.pyproj as pyproj 
import copy
import string


################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################

class replace:

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
# Convert lat/lon back to BNG
################################################################################
	def ll_to_bng(self,x,y):
		newx, newy = pyproj.transform(replace().etrs89, replace().osgb36, x, y )

		return newx, newy

	
################################################################################
# Get distance between points (assuming 1km rectiliniar grid holds)
# Input values are eastings/northings
################################################################################
	def dist_bng(self,coords1,coords2):

		diff = coords1-coords2
		dist = np.sqrt((diff*diff).sum())

		return dist

################################################################################
# Get distance between points squared (avoids spending computing power on sqrt
# but as dist is 1:1 with dist squared, we can still test on dist squared 
# instead of dist
# Input values are eastings/northings
################################################################################
	def dist_sq_bng(self,coords1,coords2):

		diff = coords1-coords2
		dist = (diff*diff).sum()

		return dist

################################################################################
# Spherical law of cosines to get distances in lat/lon space, assuming 
# spherical earth
# http://www.movable-type.co.uk/scripts/latlong.html
# Input are lat/lon pairs
################################################################################
	def dist_ll(self,coords1,coords2):

		rads1=np.pi*coords1.astype('double')/180.0
		rads2=np.pi*coords2.astype('double')/180.0
		R=6371.0 # radius of earth in km
		dlon=np.min(np.abs(rads1[1]-rads2[1]),np.abs(rads2[1]-rads1[1]))
		dist = np.arccos( (np.sin(rads1[0])*np.sin(rads2[0])) + \
				(np.cos(rads1[0])*np.cos(rads2[0])* \
				np.cos(dlon)) )*R

		return dist


################################################################################
# Parse input
################################################################################

	def parse_input(self):

		parser=argparse.ArgumentParser(description='Replace missing values with nearest value')

		# optional
		parser.add_argument('--missingval',help='Missing value',default=-1e20,type=float,required=False)
		parser.add_argument('--latthresh',help='Maximum distance in lat to search for values',default=90.0,type=float,required=False)
		parser.add_argument('--lonthresh',help='Maximum distance in lon to search for values',default=180.0,type=float,required=False)
		parser.add_argument('--monthlyfiles',help='Files are monthly (otherwise assume annual)',required=False,action='store_true')


		# positional
		parser.add_argument('infiletmplt',help='Input file template')
		parser.add_argument('ptfile',help='Points file')		
		parser.add_argument('outfiletmplt',help='Output file template')
		parser.add_argument('startyr',help='First year',type=int)
		parser.add_argument('endyr',help='Last year',type=int)
		parser.add_argument('varnames',nargs='+',help='Name of variables to replace values in')

		# Parse the arguments
		args=parser.parse_args()

		return string.Template(args.infiletmplt), \
				args.ptfile, \
				string.Template(args.outfiletmplt), \
				args.varnames, \
				args.missingval, \
				args.startyr, \
				args.endyr, \
				args.latthresh, \
				args.lonthresh, \
				args.monthlyfiles


################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################


if __name__=='__main__':

	# Call the class
	r=replace()

	infiletmplt, ptfile, outfiletmplt, varnames, mv, startyr, endyr, lat_thresh, lon_thresh, monthly = r.parse_input()

	# First open the points file
	ptf=open(ptfile,'r')
	missing=[]
	replacements={}
	for lin in ptf.readlines():
		els=lin.strip().split()
		m=int(els[0])
		missing.append(m)
		replacements[m]=[ int(r) for r in els[3::3] ]

	if monthly:
		nmn=12
	else:
		nmn=1

	# Loop over years
	for yr in range(startyr, endyr+1):
		print "Year = ", yr
		for mn in range(1,nmn+1):
			# First open the map file and get the soil indices required
			infile=infiletmplt.substitute(yr=yr,mn="%02d"%mn)
			print "Opening file: ", infile
			inf=nc.Dataset(infile,'r')

			indata={}
			for var in varnames:
				indata[var]=inf.variables[var][:]

			nt,ny,nx=indata[varnames[0]].shape

			if yr==startyr and mn==1:
				inlat=inf.variables['lat'][:]
				inlon=inf.variables['lon'][:]
				inlon[np.where(inlon>180)]-=360.0

			outdata={}
			for var in varnames:
				outdata[var]=[indata[var][t,:].flatten() for t in range(nt)]

				for t in range(nt):
					for m in missing:
					# Just replace with the first nearest point
						outdata[var][t][m]=outdata[var][t][replacements[m][0]]
			



			outfile=outfiletmplt.substitute(yr=yr,mn="%02d"%mn)
			print "Opening output file: ", outfile

			outf=nc.Dataset(outfile,'w')

			for dim in inf.dimensions:
				outf.createDimension(dim,len(inf.dimensions[dim]))

			for var in inf.variables:
				if 'missing_value' in inf.variables[var].ncattrs():
					outf.createVariable(var,inf.variables[var].dtype, \
							inf.variables[var].dimensions, \
							fill_value=inf.variables[var].missing_value)
				elif 'fill_value' in inf.variables[var].ncattrs():
					outf.createVariable(var,inf.variables[var].dtype, \
							inf.variables[var].dimensions, \
							fill_value=inf.variables[var].fill_value)
				elif '_FillValue' in inf.variables[var].ncattrs():
					outf.createVariable(var,inf.variables[var].dtype, \
							inf.variables[var].dimensions, \
							fill_value=inf.variables[var]._FillValue)
				else:
					outf.createVariable(var,inf.variables[var].dtype, \
							inf.variables[var].dimensions )

				for attr in inf.variables[var].ncattrs():
					if attr not in ('missing_value','fill_value','_FillValue'):
						outf.variables[var].setncattr(attr,inf.variables[var].getncattr(attr))

				if var in varnames:
					outf.variables[var][:]=np.array(outdata[var][:])
				else:
					outf.variables[var][:]=inf.variables[var][:]


			for attr in inf.ncattrs():
				outf.setncattr(attr,inf.getncattr(attr))

			outf.history="Missing values in %s replaced with nearest neighbours"%infile
			outf.close()
			inf.close()

	print "Finished"
