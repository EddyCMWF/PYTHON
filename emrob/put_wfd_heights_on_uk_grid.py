#!/usr/bin/env python
################################################################################
#
# This reads the WFD height data then extracts to a 2d grid covering the UK
#
# ELR 15/05/2014
#
################################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

class heights:

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
# Reshape the vector to the full grid
################################################################################
	def land_to_grid(self,data,indx,wt,missingval):

		outdata=np.ones(indx.shape)*missingval

		iy,ix=np.where(indx!=missingval)

		for i in range(len(ix)):
			outdata[iy[i],ix[i]]=data[int(indx[iy[i],ix[i]])-1]

		outdata=np.ma.masked_equal(outdata,missingval)
		outdata=outdata[::-1,:]
		return outdata

################################################################################
################################################################################
#
# Start the main routine
#
################################################################################
################################################################################

if __name__=='__main__':

	h=heights()

	htfile='/prj/WATCH/WB1/data/0p5deg/ancil/data/WFD-land-lat-long-z.nc'

	print "Reading heights from "+htfile
	htf=nc.Dataset(htfile,'r')

	ht=htf.variables['Z'][:]
	lat=htf.variables['Latitude'][:]
	lon=htf.variables['Longitude'][:]

	pdef='/prj/WATCH/WB1/data/0p5deg/ancil/data/WFD_0p5deg_pdefData.gra'

	nx=720
	ny=280
	indx,wt=h.read_pdef(pdef,nx,ny,bswap=True)

	ht_grid=h.land_to_grid(ht,indx,wt,-999)
	lat_grid=h.land_to_grid(lat,indx,wt,-999)
	lon_grid=h.land_to_grid(lon,indx,wt,-999)

	# Get the vectors by assuming that it's been put on the grid correctly
	lat_v=lat_grid.mean(axis=1)
	lon_v=lon_grid.mean(axis=0)

	lonmin=-8.25
	lonmax=2.25
	xi=np.where(np.logical_and(lon_v>=lonmin,lon_v<=lonmax))

	latmin=49.25
	latmax=62.25
	yi=np.where(np.logical_and(lat_v>=latmin,lat_v<=latmax))

	xi_grid,yi_grid=np.meshgrid(xi,yi)

	ht_grid_out=ht_grid[yi_grid,xi_grid]

	nx_out=len(xi[0])
	ny_out=len(yi[0])

	lon_out=lon_v[xi]
	lat_out=lat_v[yi]

	outfile='/local/localscratch/emrobi/CHESS/data/1km/v3/input_data/wfd/WFD-lat-lon-z-uk_grid.nc'
	outf=nc.Dataset(outfile,'w')

	outf.createDimension('x',nx_out)
	outf.createDimension('y',ny_out)
	outf.createDimension('t',1)

	mv=-1e20
	outf.createVariable('Z','f',('t','y','x'),fill_value=mv)
	# Flip the output to match the others
	outf.variables['Z'][:]=ht_grid_out[::-1,:]
	outf.variables['Z'].units="m above mean sea level"
	outf.variables['Z'].valid_min=ht_grid_out.min()
	outf.variables['Z'].valid_max=ht_grid_out.max()


	outf.createVariable('time','i',('t',))
	outf.variables['time'][:]=0

	outf.createVariable('Latitude','f',('y',))
	outf.variables['Latitude'][:]=lat_out[::-1]
	outf.variables['Latitude'].units='degrees north'
	outf.variables['Latitude'].note='pixel centre latitudes'

	outf.createVariable('Longitude','f',('x',))
	outf.variables['Longitude'][:]=lon_out[:]
	outf.variables['Longitude'].units='degrees east'
	outf.variables['Longitude'].note='pixel centre longitudes'

	outf.note="Data extracted from "+htfile


