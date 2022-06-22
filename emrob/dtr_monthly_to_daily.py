#!/usr/bin/env python
################################################################################
#
# Expand the monthly DTR files to be daily
#
# ELR 20/05/2014
#
################################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse
import string
import calendar as cal
import datetime as dt
import mpl_toolkits.basemap.pyproj as pyproj 


class dtr:

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

		parser=argparse.ArgumentParser(description='Make monthly DTR files into daily')

		# optional
		parser.add_argument('--datavar',help='Variable name for extending', required=False, default='dtr')
		parser.add_argument('--timevar',help='Time variable name', required=False, default='time')
		parser.add_argument('--latvar',help='Latitude variable name', required=False, default='lat')
		parser.add_argument('--lonvar',help='Longitude variable name', required=False, default='lon')
		parser.add_argument('--timedim',help='Time dimension name', required=False, default='time')
		parser.add_argument('--compress',help='Make compressed netcdf files', required=False, action='store_true')
		parser.add_argument('--yearorig',help='Year for time origin', required=False, type=int, default=-1)
		parser.add_argument('--mv',help='Missing value for output', required=False, type=float, default=-99999.0)

		# These default to the DTR data, but can be overwritten if
		# necessary
		parser.add_argument('--standardname',help='standard_name attribute for output data', required=False, default='daily_temperature_range')
		parser.add_argument('--units',help='units attribute for output data', required=False, default='K')
		parser.add_argument('--longname',help='long_name attribute for output data', required=False, default='Daily temperature range')
		parser.add_argument('--comment',help='comment attribute for output data', required=False)
		parser.add_argument('--description',help='discription attribute for file', required=False, default="Daily temperature range (K) at 1km resolution over Great Britain, regridded from the monthly CRU data") 
		# positional
		parser.add_argument('inftmplt',help='Template for input files')
		parser.add_argument('outftmplt',help='Template for output files')
		parser.add_argument('startyr',help='First year to convert',type=int)
		parser.add_argument('endyr',help='Last year to convert',type=int)


		# Parse the arguments
		args=parser.parse_args()

		if args.yearorig == -1:
			yearorig=args.startyr
		elif args.yearorig < 0:
			print "Error: year of origin must be positive integer"
			sys.exit(2)
		else:
			yearorig=args.yearorig

		return string.Template(args.inftmplt), \
				string.Template(args.outftmplt), \
				args.startyr, \
				args.endyr, \
				args.datavar, \
				args.timevar, \
				args.latvar, \
				args.lonvar, \
				args.timedim, \
				args.compress, \
				yearorig, \
				args.mv, \
				args.standardname, \
				args.units, \
				args.longname, \
				args.comment, \
				args.description


################################################################################
################################################################################
#
# Start the main routine
#
################################################################################
################################################################################

if __name__=='__main__':

	d=dtr()

	inftmplt, outftmplt, startyr, endyr, datavar, timevar, latvar, lonvar, timedim, compress, yearorig, mv, standardname, units, longname, comment, description = d.parse_input()

	torig=dt.datetime(yearorig,1,1)
	for yr in range(startyr, endyr+1):

		print "Year = ", yr

		for mn in range(1,13):
			print "   month = ", mn

			infname=inftmplt.substitute(yr=yr,mn="%02d"%mn)
			inf=nc.Dataset(infname,'r')

			# Get lat, lon, x, y for all files
			if yr==startyr and mn==1:
				ny,nx=inf.variables['lat'].shape

				eastings=np.arange(0,nx*1000,1000)+500
				northings=np.arange(0,ny*1000,1000)+500
				
				eastings_g, northings_g = np.meshgrid(eastings,northings)

				lon, lat = d.bng_to_ll(eastings_g.flatten(),northings_g.flatten())

				lon = lon.reshape([ny,nx])
				lat = lat.reshape([ny,nx])



			outfname=outftmplt.substitute(yr=yr,mn="%02d"%mn)
			outf=nc.Dataset(outfname,'w')

			nt=cal.monthrange(yr,mn)[1]

			for dim in inf.dimensions:
				if dim==timedim:
					outf.createDimension(dim,nt)
				else:
					outf.createDimension(dim,len(inf.dimensions[dim]))

			# First create x and y and crs (which aren't in the 
			# original file)
			outf.createVariable('x','f',('x',),fill_value=mv)
			outf.variables['x'].units='m'
			outf.variables['x'].long_name = "easting - OSGB36 grid reference" 
			outf.variables['x'].standard_name = "projection_x_coordinate" 
			outf.variables['x'].units='m'
			outf.variables['x'].point_spacing = "even" 
			outf.variables['x'][:]=eastings

			outf.createVariable('y','f',('y',),fill_value=mv)
			outf.variables['y'].units='m'
			outf.variables['y'].long_name = "northing - OSGB36 grid reference" 
			outf.variables['y'].standard_name = "projection_y_coordinate" 
			outf.variables['y'].units='m'
			outf.variables['y'].point_spacing = "even" 
			outf.variables['y'][:]=northings

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


			for var in inf.variables:
				if compress and var==datavar:
					zlib=True
					complevel=6
					shuffle=True
				else:
					zlib=False
					complevel=0
					shuffle=False


				outf.createVariable(var,inf.variables[var].dtype,inf.variables[var].dimensions,fill_value=mv,zlib=zlib,complevel=complevel,shuffle=shuffle)

				if var == datavar:
					indata=inf.variables[var][:]
					outdata=np.array([ indata for t in range(nt) ])
					outf.variables[var][:]=outdata[:]

					outf.variables[var].standard_name=standardname
					outf.variables[var].units=units
					outf.variables[var].long_name=longname
					outf.variables[var].grid_mapping='crs'
					if comment:
						outf.variables[var].comment=comment
				elif var == timevar:
					
					for t in range(nt):
						outf.variables[var][:]=[ (dt.datetime(yr,mn,t+1)-torig).days for t in range(nt) ]

						outf.variables['time'].units = "days since %d-01-01"%yearorig 
						outf.variables['time'].calendar = "gregorian" ;
						outf.variables['time'].long_name = "Time in days days since %d-01-01"%yearorig 

				elif var == latvar:
					outf.variables[var][:]=lat[:]
					outf.variables[var].long_name = "latitude of grid box centre" 
					outf.variables[var].standard_name = "latitude" 
					outf.variables[var].units = "degrees_north" 

				elif var == lonvar:
					outf.variables[var][:]=lon[:]
					outf.variables[var].long_name = "longitude of grid box centre" 
					outf.variables[var].standard_name = "longitude" 
					outf.variables[var].units = "degrees_east" 

			outf.title = longname
			outf.description = description
			outf.institution = "CEH Wallingford - NERC" 
			now=dt.datetime.now()
			outf.history = "Created "+now.strftime('%Y-%M-%d %h:%m:%s')
			outf.date_created = now.strftime('%Y-%M-%d')


			outf.geospatial_lat_min = lat.min()
			outf.geospatial_lat_max = lat.max()
			outf.geospatial_lon_min = lon.min()
			outf.geospatial_lon_max = lon.max()
			outf.version = "beta_version" 
			outf.Conventions = "CF-1.6" 
			outf.creator_name = "Emma L. Robinson" 
			outf.creator_email = "emrobi@ceh.ac.uk" 

			outf.close()
