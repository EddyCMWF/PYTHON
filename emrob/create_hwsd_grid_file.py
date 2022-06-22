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
	def corners(self,lats,lons,resx,resy,redcy=False,method='Emma'):

		rx2=resx/2.0
		ry2=resy/2.0
		
                if (method == 'Emma'):
                        # Emma looping method
                        loncorn = np.array([ [lon-rx2,lon+rx2,lon+rx2,lon-rx2] for lon in lons ])
                        latcorn = np.array([ [lat-ry2,lat-ry2,lat+ry2,lat+ry2] for lat in lats ])
                elif (method == 'Eddy'):
                        # Eddy array method
                        loncorn=np.array([lons-rx2,lons+rx2,lons+rx2,lons-rx2],dtype='float32')
                        loncorn=loncorn.T
                        latcorn=np.array([lats-ry2,lats-ry2,lats+ry2,lats+rx2],dtype='float32')
                        latcorn=latcorn.T
                else:
                        print 'Unrecognised Corner Method, using Emma method'
                        # Emma looping method
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
        lon=inf.variables['lon'][:]
        lon=lon.astype('float32')
        lat=inf.variables['lat'][:]
        lat=lat.astype('float32')

        nx=len(lon)
        ny=len(lat)

        glon,glat=np.meshgrid(lon,lat)
        
        glon=glon.astype('float32')
        glat=glat.astype('float32')
        
        glon=g.grid_to_vect(glon)
        glat=g.grid_to_vect(glat)
        
        #index area of interest, i.e. EMEP grid
        index =  ( (glon<=-50.) & (glat>=79) ) | \
                 ( (glon>-50.) & (glon<=-40.) & (glat>=65.) )  | \
                 ( (glon>-40.) & (glon<=-12.) & (glat>=60.) )  | \
                 ( (glon>-30.) & (glon<=0.)   & (glat>=30.) )  | \
                 ( (glon>0.)   & (glon<=40.)  & (glat>=20.) )  | \
                 ( (glon>40.)  & (glon<=90.)  & (glat>=30.) )  | \
                 ( (glon>90.)  & (glon<=110.) & (glat>=45.) )  | \
                 ( (glon>110.) & (glon<=160.) & (glat>=60.) )  | \
                 ( (glon>160.) & (glat>=70.) )                 & \
                   (glat < 85.)
        
        glon=glon[index]
        glat=glat[index]
        print len(glon)
        
#        indata=inf.variables['mu'][:]
#        indata=(g.grid_to_vect(indata.mask))[index]
#        gmask=np.ones_like(glon,dtype='int')
#        gmask[np.where(indata)]=0.0
        inf.close()
#        del indata
        
        
	# get corners
        resx=np.abs(lon[1]-lon[0])
        resy=np.abs(lat[1]-lat[0])
        del lon
        del lat
        glatcorn,gloncorn=g.corners(glat,glon,resx,resy,method='Eddy')
        
	# Open output file
#        print "opening file: "+outfile
#        outf=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
#
#	# create dimensions
        grid_size = len(glon)
#        outf.createDimension('grid_size',grid_size)
#        outf.createDimension('grid_corners',4)
#        outf.createDimension('grid_rank',2)
#
#        # create variables
#        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
#        gdims[:]=[nx,ny]
#
#        gctr_lat=outf.createVariable('grid_center_lat','float32',('grid_size',))
#        gctr_lat.units='degrees'
#        gctr_lat[:]=glat
#        del glat
#        
#        gctr_lon=outf.createVariable('grid_center_lon','float32',('grid_size',))
#        gctr_lon.units='degrees'
#        gctr_lon[:]=glon
#        del glon
#        
#        gim=outf.createVariable('grid_imask','i',('grid_size',))
#        gim.units='unitless'
#        gim[:]=gmask.astype('int')
#        del gmask
#
#	outf.title="%.9f x %.9f degree"%(resx,resy)
#	outf.conventions='SCRIP'
#	outf.close()

     
	# Open corners output file
        print "opening file: "+outfile[:-3]+'_latcorners.nc'
        outf=nc.Dataset(outfile[:-3]+'_latcorners.nc','w',format='NETCDF3_CLASSIC')

	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)

        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=[nx,ny]

        gcrn_lat=outf.createVariable('grid_corner_lat','float32',('grid_size','grid_corners'))
        gcrn_lat.units='degrees'
        gcrn_lat[:]=glatcorn
        del glatcorn
        
        outf.title="%.9f x %.9f degree"%(resx,resy)
	outf.conventions='SCRIP'
	outf.close()


	# Open corners output file
        print "opening file: "+outfile[:-3]+'_loncorners.nc'
        outf=nc.Dataset(outfile[:-3]+'_loncorners.nc','w',format='NETCDF3_CLASSIC')

	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)

        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=[nx,ny]

        gcrn_lon=outf.createVariable('grid_corner_lon','float32',('grid_size','grid_corners'))
        gcrn_lon.units='degrees'
        gcrn_lon[:]=gloncorn
        del gloncorn
        
	outf.title="%.9f x %.9f degree"%(resx,resy)
	outf.conventions='SCRIP'
	outf.close()











	

