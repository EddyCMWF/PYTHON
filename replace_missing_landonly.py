#!/usr/bin/env python

################################################################################
#
# Script to replace missing values in a soil file with the nearest soil
#
# ELR 17/02/2014
#
# Modifications:
#
# ECP - 03/02/1985 - to work with EMEP 50km conversion
#  
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

import netCDF4 as nc
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plot
#matplotlib.use('Agg')
import mpl_toolkits.basemap.pyproj as pyproj 
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
# Useful stuff
################################################################################
	def __init__(self):
		# projection from:
		# http://spatialreference.org/ref/epsg/27700/proj4/
		self.osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 
#
		# http://spatialreference.org/ref/epsg/4258/
		self.etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')
#
################################################################################
# Convert lat/lon back to BNG
################################################################################
	def ll_to_bng(self,x,y):
		newx, newy = pyproj.transform(replace().etrs89, replace().osgb36, x, y )
#
		return newx, newy
#	
################################################################################
# Get distance between points (assuming 1km rectiliniar grid holds)
# Input values are eastings/northings
################################################################################
	def dist_bng(self,coords1,coords2):
#
		diff = coords1-coords2
		dist = np.sqrt((diff*diff).sum())
#
		return dist
#
################################################################################
# Get distance between points squared (avoids spending computing power on sqrt
# but as dist is 1:1 with dist squared, we can still test on dist squared 
# instead of dist
# Input values are eastings/northings
################################################################################
	def dist_sq_bng(self,coords1,coords2):
#
		diff = coords1-coords2
		dist = (diff*diff).sum()
#
		return dist
#
################################################################################
# Spherical law of cosines to get distances in lat/lon space, assuming 
# spherical earth
# http://www.movable-type.co.uk/scripts/latlong.html
# Input are lat/lon pairs
################################################################################
	def dist_ll(self,coords1,coords2):
#
		rads1=np.pi*coords1.astype('double')/180.0
		rads2=np.pi*coords2.astype('double')/180.0
		R=6371.0 # radius of earth in km
		dlon=np.min(np.abs(rads1[1]-rads2[1]),np.abs(rads2[1]-rads1[1]))
		dist = np.arccos( (np.sin(rads1[0])*np.sin(rads2[0])) + \
				(np.cos(rads1[0])*np.cos(rads2[0])* \
				np.cos(dlon)) )*R
#
		return dist
#
################################################################################
# Parse input
################################################################################
	def parse_input(self):
#
		parser=argparse.ArgumentParser(description='Replace missing values with nearest value')
#
		# optional
		parser.add_argument('--mapindxname',help='Name of index variable to read from the map file',default='mu',required=False)
		parser.add_argument('--maskname',help='Name of mask variable to read from the mask file',default='landfrac',required=False)
		parser.add_argument('--missingval',help='Missing value',default=-9999.0,type=float,required=False)
		parser.add_argument('--EMEP',help='EMEP grid flag',default=False,type=bool,required=False)
#
		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('maskfile',help='Land mask file')
		parser.add_argument('outfile',help='Output file')
#
		# Parse the arguments
		args=parser.parse_args()
#
		return args.infile, \
			args.maskfile, \
			args.outfile, \
			args.mapindxname, \
			args.maskname, \
			args.missingval, \
			args.EMEP
#

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

	infile, maskfile, outfile, miname, maskname, mv, EMEP = r.parse_input()

        # First open the map file and get the soil indices required
	print "Opening file: ", infile
        inf=nc.Dataset(infile,'r')
        indata=inf.variables[miname][:].flatten()
        inlat=inf.variables['lat'][:].flatten()
        inlon=inf.variables['lon'][:].flatten()
        inlon[np.where(inlon>180)]-=360.0
        npts=len(inlon)
        
	#ineasting=np.zeros_like(inlon)
	#innorthing=np.zeros_like(inlon)
        
        # Then open the mask file
	print "Opening mask: ", maskfile
        mf=nc.Dataset(maskfile,'r')
        maskdata=mf.variables[maskname][:].flatten()
        
        indata.mask[maskdata==0] = True
        #indata.data[maskdata==0] = -9999.
        
	# Get coordinates from the mask file
        if EMEP == True:
                eastings,northings=mf.variables['i'][:],mf.variables['j'][:]
                i_EMEP,j_EMEP = mf.variables['grid_EMEP_i_coord'][:],mf.variables['grid_EMEP_j_coord'][:]
                print i_EMEP.shape, j_EMEP.shape
        else:
                eastings,northings=np.meshgrid(mf.variables['x'][:],mf.variables['y'][:])
	
	coords = [ np.array([eastings.flatten()[i],northings.flatten()[i]]) for i in range(npts) ]

	#for i in range(npts):
		##ineasting[i],innorthing[i]=coords[i]
#
		## These are grid centres, so round to nearest 100 to avoid
		## some of the precision issues
		##ineasting[i]=np.round(ineasting[i],-2)
		##innorthing[i]=np.round(innorthing[i],-2)
		#coords[i][0]=np.round(coords[i][0],-2)
		#coords[i][1]=np.round(coords[i][1],-2)
#
	# Find the missing values and loop over them
	print "Finding missing data"
	# Data need to be filled if they're missing in the input data and have
	# a non-zero value in the mask
	missing=np.where(np.logical_and(indata.mask,maskdata>0,indata.data>=0.))[0]

	# Find points within certain distance of each missing point
	print "Finding points for which to calculate distances"
	latlon_thresh=5.
	include = [ np.where(np.logical_and( 
                np.abs(inlon[m]-inlon)<=latlon_thresh,
                np.abs(inlat[m]-inlat)<=latlon_thresh ))[0] for m in missing ]

	# Get the distances (squared)
	print "Calculating distances"
	dists=[np.array([ r.dist_sq_bng(coords[missing[mi]],coords[i]) for i in include[mi] ]) for mi in range(len(missing)) ]

	print "Replace missing data (%d points)"%len(missing)
	outdata=copy.deepcopy(indata)
	for mi in range(len(missing)):
		m=missing[mi]
		minc=np.where(include[mi]==m)
		print "  %d [%d]"%(mi,m)
		print "--------------------------------------------------------------------------------"
		#print "indata[%d]=%f"%(m,indata[m])



		# Make sure we don't pick the point itself
		dists[mi][minc]=1e16

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
                        if (len(dists[mi])==1):        #if only 1 point found within lat/lon threshold jump to next point
                                ended=True
			elif mindist>=1e16:
				ended=True
			else:
				print "==="
				print "mindist = %f"%mindist
				mindistpoints=np.where(dists[mi]==mindist)

				mindistdata=indata[include[mi][mindistpoints]]
				print "mindistdata = ", mindistdata
				print "mask = ", all(mindistdata.mask)
				if all(mindistdata.mask):
					dists[mi][mindistpoints]=1e16
				else:
					count=np.array([len(np.where(mindistdata==d)[0]) for d in np.unique(mindistdata[np.where(np.logical_not(mindistdata.mask))])])
					if len(count)!=len(np.unique(count)):
						print "Warning: multiple possible values for point %d:"%m
						print np.unique(mindistdata[np.where(np.logical_not(mindistdata.mask))])

					outdata[m]=np.unique(mindistdata[np.where(np.logical_not(mindistdata.mask))])[np.argmax(count)]


					print "outdata[%d]=%f"%(m,outdata[m])
					found=True


		if ended and not found:
			print "Warning: no unmasked data found for this point"

		mi+=1

	# Extra step is to mask out the sea points, just to be sure

	outdata=np.ma.masked_where(maskdata==0,outdata)
        print 'Writing to: '+outfile
	outf=nc.Dataset(outfile,'w')
        
	for dim in inf.dimensions:
		outf.createDimension(dim,len(inf.dimensions[dim]))
	for var in inf.variables:
                print var
		if 'missing_value' in inf.variables[var].ncattrs():
                        print inf.variables[var].missing_value
			outf.createVariable(var,inf.variables[var].dtype, \
					inf.variables[var].dimensions, \
					fill_value=inf.variables[var].missing_value)
		elif '_FillValue' in inf.variables[var].ncattrs():
			outf.createVariable(var,inf.variables[var].dtype, \
					inf.variables[var].dimensions, \
					    fill_value=inf.variables[var]._FillValue)
                else:
			outf.createVariable(var,inf.variables[var].dtype, \
					inf.variables[var].dimensions )

		if var==miname:
			outf.variables[var][:]=outdata[:]
		else:
			outf.variables[var][:]=inf.variables[var][:]
        if (EMEP==True):
                if not ('i_EMEP' in outf.variables):
                        outf.createVariable('i_EMEP','float32',inf.variables['lon'].dimensions )
                        outf.variables['i_EMEP'][:]=i_EMEP[:]
                if not ('j_EMEP' in outf.variables):
                        outf.createVariable('j_EMEP','float32',inf.variables['lat'].dimensions )
                        outf.variables['j_EMEP'][:]=j_EMEP[:]
        
	outf.history="Missing values in %s replaced with LAF of nearest neighbours"%infile
	outf.close()
	print "Finished"
