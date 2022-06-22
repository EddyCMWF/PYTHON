#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse
import string
import calendar as cal
import mpl_toolkits.basemap.pyproj as pyproj 
import datetime as dt


class rainfall:

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

		parser=argparse.ArgumentParser(description='Trim rainfall files to new grid')

		# optional
		parser.add_argument('-s','--scale',help='Scale to mm/s from mm/day', required=False, action='store_true')
		parser.add_argument('-o','--outvarname',help='Variable name for output data', required=False, default=False)
		parser.add_argument('-O','--outmissing',help='Missing value for output data', required=False, default=False,type=float)
		parser.add_argument('-u','--units',help='Units for output', required=False)
		parser.add_argument('-e','--extramask',help='Extra values for masking',default='', required=False)

		# positional
		parser.add_argument('inftmplt',help='Template for input files')
		parser.add_argument('landmaskf',help='Land mask file')
		parser.add_argument('outftmplt',help='Template for output files')
		parser.add_argument('startyr',help='First year to convert',type=int)
		parser.add_argument('endyr',help='Last year to convert',type=int)

		# Parse the arguments
		args=parser.parse_args()

		if args.extramask=='':
			extramask=[]
		else:
			extramask=[float(e) for e in args.extramask.split(',')]

		return string.Template(args.inftmplt), \
				args.landmaskf, \
				string.Template(args.outftmplt), \
				args.startyr, \
				args.endyr, \
				args.scale, \
				args.outvarname, \
				args.outmissing, \
				args.units, \
				extramask

################################################################################
################################################################################
#
# Start the main routine
#
################################################################################
################################################################################

if __name__=='__main__':

	r=rainfall()


	inftmplt, landmaskf, outftmplt, startyr, endyr, scale, outvar, outmv, outunits, extramask = r.parse_input()

	# open land mask
	lf=nc.Dataset(landmaskf,'r')
	landfrac=np.ma.masked_less_equal(lf.variables['landfrac'][:],0.0)
	ny,nx=landfrac.shape

	eastings=lf.variables['x'][:]
	northings=lf.variables['y'][:]

	eastings_g, northings_g = np.meshgrid(eastings,northings)

	#print "Get lat/lons"
	lon, lat = r.bng_to_ll(eastings_g.flatten(),northings_g.flatten())

	lon = lon.reshape([ny,nx])
	lat = lat.reshape([ny,nx])


	datavar='rainfall_amount'
	if not outvar:
		outvar=datavar

	gridvars=['x','y','time','lat','lon']
	griddims=['x','y','time']
	dimlens={ 'x':nx, \
			'y':ny, \
			'time':0 }

	torig=dt.datetime(startyr,1,1)
	origorig=dt.datetime(1800,1,1)
	tshift=(torig-origorig).days

	for yr in range(startyr, endyr+1):

		print "Year = ", yr

		infname=inftmplt.substitute(yr=yr)

		inf=nc.Dataset(infname,'r')

		if not outmv:
			outmv=inf.variables[datavar]._FillValue

		iend=0
		for mn in range(1,13):

			print "  %02d"%mn
			ndays=cal.monthrange(yr,mn)[1]
			istart=iend
			iend=istart+ndays

			nt=iend-istart
			dimlens['time']=nt

			dmask=np.array([landfrac.mask for t in range(nt)])
			# Don't need the scale factor as python uses it
			# implicitly when reading
			#indata=np.ma.masked_where(dmask,inf.variables[datavar][istart:iend,:-(ny+1):-1,:nx]*inf.variables[datavar].scale_factor*scale)

			indata=np.ma.masked_where(dmask,inf.variables[datavar][istart:iend,:-(ny+1):-1,:nx])
			for em in extramask:
				indata=np.ma.masked_equal(indata,em)

			if scale:
				indata=indata/86400.0 

			timearr=inf.variables['time'][istart:iend]


			outfname=outftmplt.substitute(yr=yr,mn="%02d"%mn)
			outf=nc.Dataset(outfname,'w')

			# Set up dimensions
			for dim in griddims:
				outf.createDimension(dim,dimlens[dim])

			# Set up coordinate variables
			outf.createVariable('time','f',('time',),fill_value=outmv)
			outf.variables['time'].units = "days since "+torig.strftime('%Y-%m-%d')
			outf.variables['time'].calendar = "gregorian" ;
			outf.variables['time'].long_name = "Time in days days since "+torig.strftime('%Y-%m-%d')
			outf.variables['time'][:]=timearr-tshift

			outf.createVariable('y','f',('y',),fill_value=outmv)
			outf.variables['y'].units='m'
			outf.variables['y'].long_name = "northing - OSGB36 grid reference" 
			outf.variables['y'].standard_name = "projection_y_coordinate" 
			outf.variables['y'].units='m'
			outf.variables['y'].point_spacing = "even" 
			outf.variables['y'][:]=np.arange(0,dimlens['y']*1000,1000)+500.0

			outf.createVariable('x','f',('x',),fill_value=outmv)
			outf.variables['x'].units='m'
			outf.variables['x'].long_name = "easting - OSGB36 grid reference" 
			outf.variables['x'].standard_name = "projection_x_coordinate" 
			outf.variables['x'].units='m'
			outf.variables['x'].point_spacing = "even" 
			outf.variables['x'][:]=np.arange(0,dimlens['x']*1000,1000)+500.0

			outf.createVariable('lat','f',('y','x'),fill_value=outmv)
			outf.variables['lat'].long_name = "latitude of grid box centre" 
			outf.variables['lat'].standard_name = "latitude" 
			outf.variables['lat'].units = "degrees_north" 

			outf.createVariable('lon','f',('y','x'),fill_value=outmv)
			outf.variables['lon'].long_name = "longitude of grid box centre" 
			outf.variables['lon'].standard_name = "longitude" 
			outf.variables['lon'].units = "degrees_east" 
					
			outf.variables['lat'][:]=lat[:]
			outf.variables['lon'][:]=lon[:]

			# Set up remapping variable
			outf.createVariable('crs','i')
			for attr in inf.variables['crs'].ncattrs():
				outf.variables['crs'].setncattr(attr,inf.variables['crs'].getncattr(attr))
				

			# Set up data variable(s)
			if inf.variables[datavar].dtype=='float64':
				ndtype='float32'
			elif inf.variables[datavar].dtype=='in64':
				ndtype='int32'
			else:
				ndtype=inf.variables[datavar].dtype
			outf.createVariable(outvar,ndtype,inf.variables[datavar].dimensions,fill_value=outmv,zlib=True,complevel=6,shuffle=True)

			for attr in inf.variables[datavar].ncattrs():
				if attr == 'standard_name':
					outf.variables[outvar].standard_name = 'precipitation'
				elif attr == 'units':
					if outunits:
						outf.variables[outvar].units = outunits
					else:
						outf.variables[outvar].units = inf.variables[datavar].units
				elif attr not in ['_FillValue','valid_min','valid_max','scale_factor']:
					outf.variables[outvar].setncattr(attr,inf.variables[datavar].getncattr(attr))

			outf.variables[outvar][:]=indata[:]


			# copy some global attributes
			for attr in inf.ncattrs():
				if attr not in [ 'geospatial_lat_min', 'geospatial_lat_max', 'geospatial_lon_min', 'geospatial_lon_max', 'version', 'creator_name', 'creator_email' , 'history', 'date_created', 'date_modified', 'date_issued']:
					outf.setncattr(attr,inf.getncattr(attr))
					
			now=dt.datetime.now()
			outf.history = "Created "+now.strftime('%Y-%m-%d %H:%M:%S')
			outf.date_created = now.strftime('%Y-%m-%d')


			outf.geospatial_lat_min = lat.min().astype('float32')
			outf.geospatial_lat_max = lat.max().astype('float32')
			outf.geospatial_lon_min = lon.min().astype('float32')
			outf.geospatial_lon_max = lon.max().astype('float32')
			outf.version = "beta_version" 
			outf.Conventions = "CF-1.6" 
			outf.creator_name = "Emma L. Robinson" 
			outf.creator_email = "emrobi@ceh.ac.uk" 


			outf.close()





