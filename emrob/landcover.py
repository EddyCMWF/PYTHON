#!/usr/bin/env python

# Assemble the original extracted land cover from Jon to a netcdf file

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import string

import mpl_toolkits.basemap.pyproj as pyproj 


osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

# http://spatialreference.org/ref/epsg/4258/
etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')

def bng_to_ll(x,y):
	lonlat = [ pyproj.transform(osgb36, etrs89, x[i], y[i] ) for i in range(len(x)) ]

	lon = np.array([l[0] for l in lonlat])
	lat = np.array([l[1] for l in lonlat])

	return lon, lat

################################################################################

# read the landcover data
nxin=700
nyin=1300
nzin=9

nxout=656
nyout=1057
nzout=8

mv=-99999.0


ftmplt=string.Template('/users/rec/jon/JULESdrv/ancil/lcov2000_${p}.bin')

indata=[]
for i in range(nzin):
	if i!=3:
		ida=np.fromfile(ftmplt.substitute(p='%02d'%(i+1)),dtype='float32').byteswap().reshape([nyin,nxin])
		indata.append(ida[:nyout,:nxout])

indata=np.ma.masked_equal(np.array(indata),mv)


################################################################################
# read the land mask
lf=nc.Dataset('/local/localscratch/emrobi/CHESS/data/1km/v3/ancil/chess_landfrac.nc','r')
mask=np.ma.masked_equal(lf.variables['landfrac'][:],0)

for i in range(nzout):
	indata[i,:]=np.ma.masked_where(mask.mask,indata[i,:])

################################################################################
# write output
outf=nc.Dataset('/local/localscratch/emrobi/CHESS/data/1km/v3/ancil/chess_landcover_2000.nc','w')

outf.createDimension('z',nzout)
outf.createDimension('y',nyout)
outf.createDimension('x',nxout)

outf.createVariable('z','f',('z',),fill_value=mv)
outf.variables['z'].units=''
outf.variables['z'].long_name = "land surface type"
outf.variables['z'].standard_name = "land_surface_type" 
outf.variables['z'].axis = 'z'
outf.variables['z'][:]=np.arange(1,nzout+1)

outf.createVariable('y','f',('y',),fill_value=mv)
outf.variables['y'].units='m'
outf.variables['y'].long_name = "northing - OSGB36 grid reference" 
outf.variables['y'].standard_name = "projection_y_coordinate" 
outf.variables['y'].units='m'
outf.variables['y'].point_spacing = "even" 
outf.variables['y'][:]=np.arange(0,nyout*1000,1000)+500.0

outf.createVariable('x','f',('x',),fill_value=mv)
outf.variables['x'].units='m'
outf.variables['x'].long_name = "easting - OSGB36 grid reference" 
outf.variables['x'].standard_name = "projection_x_coordinate" 
outf.variables['x'].units='m'
outf.variables['x'].point_spacing = "even" 
outf.variables['x'][:]=np.arange(0,nxout*1000,1000)+500.0

outf.createVariable('lat','f',('y','x'),fill_value=mv)
outf.variables['lat'].long_name = "latitude of grid box centre" 
outf.variables['lat'].standard_name = "latitude" 
outf.variables['lat'].units = "degrees_north" 

outf.createVariable('lon','f',('y','x'),fill_value=mv)
outf.variables['lon'].long_name = "longitude of grid box centre" 
outf.variables['lon'].standard_name = "longitude" 
outf.variables['lon'].units = "degrees_east" 

xg,yg=np.meshgrid(np.arange(0,nxout*1000,1000),np.arange(0,nyout*1000,1000))
xg+=500
yg+=500

glon,glat = bng_to_ll(xg.flatten(),yg.flatten())
glon=glon.reshape([nyout,nxout])
glat=glat.reshape([nyout,nxout])

outf.variables['lat'][:]=glat[:]
outf.variables['lon'][:]=glon[:]

# Coord variable
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

# data variable
outf.createVariable('frac','f',('z','y','x'),fill_value=mv,zlib=True,complevel=6,shuffle=True)
outf.variables['frac'].standard_name='fractional_coverage'
outf.variables['frac'].long_name='Fractional coverage'
outf.variables['frac'].units='none'
outf.variables['frac'].grid_mapping='crs'
outf.variables['frac'].coordinates='z lat lon'
outf.variables['frac'][:]=indata[:]

outf.title = "Fractional coverage" 
outf.description = "Fractional coverage of different surface types at 1km resolution over Great Britain" 
outf.institution = "CEH Wallingford - NERC" 
now=dt.datetime.now()
outf.history = "Created "+now.strftime('%Y-%M-%d %h:%m:%s')
outf.date_created = now.strftime('%Y-%M-%d')


outf.geospatial_lat_min = glat.min()
outf.geospatial_lat_max = glat.max()
outf.geospatial_lon_min = glon.min()
outf.geospatial_lon_max = glon.max()
outf.version = "beta_version" 
outf.Conventions = "CF-1.6" 
outf.creator_name = "Emma L. Robinson" 
outf.creator_email = "emrobi@ceh.ac.uk" 

outf.close()



