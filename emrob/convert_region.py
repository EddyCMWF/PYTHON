#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import string
import mpl_toolkits.basemap.pyproj as pyproj 
import datetime as dt



# projection from:
# http://spatialreference.org/ref/epsg/27700/proj4/
osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

# http://spatialreference.org/ref/epsg/4258/
etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')


################################################################################
# Conversion stuff
################################################################################

def bng_to_ll(x,y):
	lonlat = [ pyproj.transform(osgb36, etrs89, x[i], y[i] ) for i in range(len(x)) ]

	lon = np.array([l[0] for l in lonlat])
	lat = np.array([l[1] for l in lonlat])

	return lon, lat



nxin=700
nyin=1300

nxout=656
nyout=1057

regions=['eng_low1','england','scotland','wales']

raster_tmplt=string.Template('/users/global/emrobi/Data-Area/projects/HELM/runoff_files/raster/${reg}_shft.dat')

masks={}
masks['gb']=np.zeros([nyout,nxout])

for reg in regions:
	rfname=raster_tmplt.substitute(reg=reg)

	rf=np.fromfile(rfname,dtype='float32')
	rf=rf.reshape([nyin,nxin])
	masks['gb']+=rf[:nyout,:nxout]
	masks[reg]=np.ma.masked_equal(rf[:nyout,:nxout],0)

masks['gb']=np.ma.masked_equal(masks['gb'],0)
masks['gb']/=masks['gb']

northings=np.arange(500,(nyout*1000)+500,1000)
eastings=np.arange(500,(nxout*1000)+500,1000)
eastings_g, northings_g = np.meshgrid(eastings,northings)

#print "Get lat/lons"
lon, lat = bng_to_ll(eastings_g.flatten(),northings_g.flatten())

lon = lon.reshape([nyout,nxout])
lat = lat.reshape([nyout,nxout])

outfname='/prj/chess/data/1km/v3/ancil/chess_regional_masks.nc'

outf=nc.Dataset(outfname,'w')

outf.createDimension('x',nxout)
outf.createDimension('y',nyout)

mv=-99999.0 

outf.createVariable('y','f',('y',),fill_value=mv)
outf.variables['y'].units = 'm'
outf.variables['y'].long_name = \
		'northing - OSGB36 grid reference'
outf.variables['y'].standard_name = \
		'projection_y_coordinate'
outf.variables['y'].point_spacing = 'even'
outf.variables['y']=northings

outf.createVariable('x','f',('x',),fill_value=mv)
outf.variables['x'].units = 'm'
outf.variables['x'].long_name = \
		'easting - OSGB36 grid reference'
outf.variables['x'].standard_name = \
		'projection_x_coordinate'
outf.variables['x'].point_spacing = 'even'
outf.variables['x']=eastings

outf.createVariable('lat','f',('y','x'),fill_value=mv)
outf.variables['lat'].long_name = "latitude of grid box centre" 
outf.variables['lat'].standard_name = "latitude" 
outf.variables['lat'].units = "degrees_north" 

outf.createVariable('lon','f',('y','x'),fill_value=mv)
outf.variables['lon'].long_name = "longitude of grid box centre" 
outf.variables['lon'].standard_name = "longitude" 
outf.variables['lon'].units = "degrees_east" 
		
outf.variables['lat'][:]=lat[:]
outf.variables['lon'][:]=lon[:]



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


regions=regions+['gb',]

reg_long_names={'eng_low1':'English Lowlands', \
		'england':'England', \
		'wales':'Wales', \
		'scotland':'Scotland', \
		'gb':'Great Britain' }
for reg in regions:
	outf.createVariable(reg,'f',('y','x'),fill_value=mv)
	outf.variables[reg].long_name='Mask of '+reg_long_names[reg]
	outf.variables[reg].standard_name='mask_of_reg'
	outf.variables[reg].units=''
	outf.variables[reg].grid_mapping='crs'
	outf.variables[reg][:]=masks[reg][:]

outf.title='Regional land masks at 1km resolution over Great Britain'
outf.description = 'Regional land masks at 1km resolution over Great Britain'
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
