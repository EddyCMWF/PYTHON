#!/usr/bin/env python

################################################################################
#
# Script to find the nearest points with which to replace missing data
# This only does the search, not the replacement
#
# ELR 17/04/2014
#
################################################################################
#
# Added a land mask file, and accounted for time/level dimensions
# 
# ELR 21/5/2014
#
################################################################################

import netCDF4 as nc
import sys
import numpy as np
import argparse
import matplotlib
import copy


################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################

class replace:

################################################################################
# Spherical law of cosines to get distances in lat/lon space, assuming 
# spherical earth
# http://www.movable-type.co.uk/scripts/latlong.html
# Input are lat/lon pairs
################################################################################
	def dist_ll(self,coords1,coords2):

		if np.all(coords1==coords2):
			dist=0.0
		else:
			rads1=np.pi*coords1.astype('double')/180.0
			rads2=np.pi*coords2.astype('double')/180.0
			R=6371.0 # radius of earth in km
			dlon=rads2[1]-rads1[1]
			dist = np.arccos( (np.sin(rads1[0])*np.sin(rads2[0])) +\
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
		parser.add_argument('--maskname',help='Name of mask variable to read from the mask file',default='landfrac',required=False)


		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('maskfile',help='Land mask file')		
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('var',help='Variable to fill')

		# Parse the arguments
		args=parser.parse_args()

		return args.infile, \
				args.maskfile, \
				args.outfile, \
				args.missingval, \
				args.var, \
				args.latthresh, \
				args.lonthresh, \
				args.maskname


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

	infile, maskfile, outfile, mv, var, latthresh, lonthresh, maskname = r.parse_input()

# First open the map file and get the variables required
	print "Opening file: ", infile
	inf=nc.Dataset(infile,'r')
	if len(inf.variables[var].shape)==2:
		indata=inf.variables[var][:]
	elif len(inf.variables[var].shape)==3:
		indata=inf.variables[var][0,:]
	elif len(inf.variables[var].shape)==4:
		indata=inf.variables[var][0,0,:]
	else:
		print "Error: too many or too few dimensions for variable "+var
		sys.exit(2)

	indata=indata.flatten()

	inlat=inf.variables['lat'][:].flatten()
	inlon=inf.variables['lon'][:].flatten()
	inlon[np.where(inlon>180)]-=360.0
	npts=len(inlon)

	# Open the mask file
	print "Opening mask: ", maskfile
	mf=nc.Dataset(maskfile,'r')
	maskdata=mf.variables[maskname][:].flatten()
	inmask=np.ma.masked_equal(maskdata,0).mask

	#ineasting=np.zeros_like(inlon)
	#innorthing=np.zeros_like(inlon)

	# Convert back to east/northings
	print "Converting lat/lon to grid ref"
	coords = [ np.array([inlat[i],inlon[i]]) for i in range(npts) ]

	# Find the missing values and loop over them
	print "Finding missing data"
	# Get the mask from the first timestep of the first variable
	# subject to the input land mask
	# (assuming that the missing points aren't time dependent
	datamask=indata.flatten().mask
	missing=np.where(np.logical_and(datamask,~inmask))[0]
	replaced=[]

	# Find points within certain distance of each missing point
	print "Finding points for which to calculate distances"
	include = [ np.where(np.logical_and( 
		np.abs(inlon[m]-inlon)<=lonthresh,
		np.abs(inlat[m]-inlat)<=latthresh ))[0] for m in missing ]

	# Get the distances 
	print "Calculating distances"
	dists=[np.array([ r.dist_ll(coords[missing[mi]],coords[i]) for i in include[mi] ]) for mi in range(len(missing)) ]

	print "Find replacement points"

	notfound=0
	replacements={}
	for mi in range(len(missing)):
		m=missing[mi]
		minc=np.where(include[mi]==m)

		# Make sure we don't pick the point itself
		dists[mi][minc]=1e20

		# Check we've not got any other points on top of the original
		if np.any(dists[mi]==0):
			print "Error: degenerate points, how did that happen?"
			sys.exit(2)

		# Find some non-missing data
		found=False
		ended=False
		while not found and not ended:
			# Find the minimum distance, and where it is
			mindist=min(dists[mi])
			if mindist==1e20:
				ended=True
			else:
				#print "==="
				#print "mindist = %f"%mindist
				mindistpoints=np.array(np.where(dists[mi]==mindist)[0])

				mindistmask=datamask[include[mi][mindistpoints]]


				if all(mindistmask):
					# All are masked, we spiral out to the
					# next nearest
					dists[mi][mindistpoints]=1e20
				else:
					# There's at least one actual value

					# We just take the first index
					replacements[m]=include[mi][mindistpoints[np.where(~mindistmask)]]

					found=True
					replaced.append(m)


		if ended and not found:
			print "Warning: no unmasked data found for this point"
			notfound+=1

		mi+=1

	print "Not found: %d"%notfound
	print "Opening output file: ", outfile

	outf=open(outfile,'w')
	
	for m in replaced:
		outf.write('%d %f %f    '%(m,inlat[m],inlon[m]))
		for r in replacements[m]:
			outf.write('%d %f %f '%(r,inlat[r],inlon[r]))
		outf.write('\n')

	outf.close()
	print "Finished"
