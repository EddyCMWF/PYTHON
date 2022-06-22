#/usr/bin/python

################################################################################
################################################################################
#
# Script to create a regular lat/lon grid file for use in SCRIP, using a land 
# mask
#
# ELR 9/9/2013
#
################################################################################
#
# Edited to go back to (-180,180) rather than (0,360) for longitude
#
# ELR 25/9/2013
#
################################################################################
# 
# Edited to allow selection of one hemisphere only
#
# ELR 22/11/2013
#
################################################################################

# Import useful things
import netCDF4 as nc
import matplotlib.pyplot as plt
import argparse
import numpy as np

################################################################################
# Class definition
################################################################################
class gridc:

################################################################################
# Flatten grid
################################################################################
	def grid_to_vect(self, grid):

		vect=grid.flatten()

		return vect

################################################################################
# Vect to grid
################################################################################
	def vect_to_grid(self, vect, nx, ny):

		grid = np.array([ vect[j*nx:(j+1)*nx] for j in range(ny) ])

		return grid

################################################################################
# Create full global grid
################################################################################
	def global_grid(self,res,yrev=False,poslon=False,hem='g'):

		if hem=='g':
			lat=np.arange(-90.0+(res/2.0),90.0,res)
		elif hem=='n':
			lat=np.arange(0.0+(res/2.0),90.0,res)
		elif hem=='s':
			lat=np.arange(-90.0+(res/2.0),0.0,res)
		else:
			print "Invalid hemisphere selection: %s"%hem
			sys.exit(2)

		if yrev:
			lat=lat[::-1]

		if poslon:
			lon=np.arange(0.0+(res/2.0),360.0,res)
		else:
			lon=np.arange(-180.0+(res/2.0),180.0,res)

		ny=len(lat)
		nx=len(lon)

		lon,lat=np.meshgrid(lon,lat)

		lat=gridc().grid_to_vect(lat)
		lon=gridc().grid_to_vect(lon)

		data=np.zeros_like(lon)
	
		return data, lat, lon, nx, ny

################################################################################
# Get the corners
################################################################################
	def corners(self,lats,lons,res,redcy=False):

		r2=res/2.0
		
		loncorn = np.array([ [lon-r2,lon+r2,lon+r2,lon-r2] for lon in lons ])
		latcorn = np.array([ [lat-r2,lat-r2,lat+r2,lat+r2] for lat in lats ])

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
# Read the pdef file
################################################################################
	def read_pdef(self,pdefname,nx,ny,bswap=False):

		# read the index (4byte ints)
		indx=np.fromfile(pdefname,dtype='int32',count=nx*ny)
		indx.byteswap(bswap)
		indx=indx.astype('float32')
		indx=indx.reshape([ny,nx])

		# read the weights (floats) 
		wt=np.fromfile(pdefname,dtype='float32',count=-1)
		wt=wt[nx*ny:(2*nx*ny)]
		wt.byteswap(bswap)
		wt=wt.reshape([ny,nx])

		return indx,wt

################################################################################
# Parse the input
################################################################################
        def parse_input(self):

		parser=argparse.ArgumentParser(description='Create a grid file for SCRIP')

		# Optional
		parser.add_argument('--infile','-i',help='Input file (containing lat/lon values on a land-vector, to make a mask)',required=False,default=None)
		parser.add_argument('--yrev','-y',help='Reverse y-axis in input file',action='store_true')
		parser.add_argument('--poslon','-p',help='Put longitude on range (0,360). If not selected, default is longitude on range (-180,180)',action='store_true')
		parser.add_argument('--nhem','-n',help='Select only northern hemisphere',action='store_true')
		parser.add_argument('--shem','-s',help='Select only southern hemisphere',action='store_true')
		parser.add_argument('--redundancy','-r',help='Make cells that intersect the poles into triangles ratehr than squares',action='store_true')
		# positional
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('resolution',help='Resolution of the output grid (because it''s not necesarrily obvious from the input file)',type=float)

		# Parse the arguments
		args=parser.parse_args()

		if args.nhem and args.shem:
			print "Error: can only restrict to either northern hemisphere or southern hemisphere, not both!"
			sys.exit(2)

		if args.nhem:
			hem='n'
		elif args.shem:
			hem='s'
		else:
			hem='g'

		return args.infile, \
				args.outfile, \
				args.resolution, \
				args.yrev, \
				args.poslon, \
				args.redundancy, \
				hem


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

	infile, outfile, res, yrev, poslon, redcy, hem = g.parse_input()

	# Full global lat/lon
	gmask, glat, glon, nx, ny = g.global_grid(res,yrev,poslon,hem)

	if infile: 
		inf=nc.Dataset(infile,'r')
		ilat=inf.variables['latitude'][:]
		ilon=inf.variables['longitude'][:]
		#clat,clon = g.corners(ilat,ilon,res)
		npts=len(ilat)

		print "finding indxs"
		indxs = [ np.where(np.logical_and(glat==ilat[i],glon==ilon[i])) for i in range(npts) ]

		for i in indxs:
			if len(i[0])>1:
				print "Error: degenerate grid"
				sys.exit(1)
			elif len(i[0])<1:
				print "Error: Point not found"
				sys.exit(1)
		
		print "getting gmask"

		for indx in indxs:
			gmask[indx] = 1
	else:
		gmask[:]=1

	print "getting corners"
	glatcorn,gloncorn=g.corners(glat,glon,res,redcy)

	print "opening file: "+outfile
	outf=nc.Dataset(outfile,'w')

	outf.createDimension('grid_size',len(gmask))
	outf.createDimension('grid_corners',4)
	outf.createDimension('grid_rank',2)

	gdims=outf.createVariable('grid_dims','i',('grid_rank',))
	gdims[:]=[nx,ny]


	gctr_lat=outf.createVariable('grid_center_lat','d',('grid_size',))
	gctr_lat.units='degrees'
	gctr_lat[:]=glat

	gctr_lon=outf.createVariable('grid_center_lon','d',('grid_size',))
	gctr_lon.units='degrees'
	gctr_lon[:]=glon

	gim=outf.createVariable('grid_imask','i',('grid_size',))
	gim[:]=gmask.astype('int')
	gim.units='unitless'

	gcrn_lat=outf.createVariable('grid_corner_lat','d',('grid_size','grid_corners'))
	gcrn_lat.units='degrees'
	gcrn_lat[:]=glatcorn

	gcrn_lon=outf.createVariable('grid_corner_lon','d',('grid_size','grid_corners'))
	gcrn_lon.units='degrees'
	gcrn_lon[:]=gloncorn

	outf.title="%.2f degree"%res

	outf.conventions='SCRIP'

	outf.close()


