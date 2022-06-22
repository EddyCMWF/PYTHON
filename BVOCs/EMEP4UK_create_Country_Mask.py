#!/usr/bin/python2.7

################################################################
# Create a Country Mask for the EMEP4UK data
#   Using the cartopy shapefiles
################################################################

import numpy as np
import netCDF4 as nc
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from cartopy.io.shapereader import natural_earth, Reader

country_names=['England','Scotland','Wales','Northern Ireland','Ireland',\
               'France','Flemish Region','Netherlands','Denmark','Norway']

# Read in the lat and lon grid from an EMEP4UK file
EMEP4UK_latlon_file = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_output/rv4.3/'+\
                      'EMEP4UK_UK_webrun_emep_4.3_Area_Grid_km2.nc'

inf = nc.Dataset(EMEP4UK_latlon_file,'r')
lats=inf.variables['lat'][:]
lons=inf.variables['lon'][:]
inf.close()
nx,ny=lats.shape

# Read in the shape files from the 
shape_records = Reader(natural_earth(resolution='10m',
                                      category='cultural',
                                      name='admin_0_map_units')).records()
geoms=[]
country_pos=[]
for country in shape_records:
    if country.attributes['NAME_LONG'] in country_names:
        print country.attributes['NAME_LONG']
        country_pos.append(country_names.index(country.attributes['NAME_LONG']))
        geoms.append(country.geometry)

# Now loop over lat and long points and see which polygon (if any) there are in
COUNTRY_MASK=np.zeros_like(lats)
for i in range(nx):
    for j in range(ny):
        coordPoint = Point(lons[i,j],lats[i,j])
        for poly,num in zip(geoms,country_pos):
            if poly.contains(coordPoint):
                COUNTRY_MASK[i,j]=num+1
                continue

# Write mask array out to file
EMEP4UK_countryfile='/users/eow/edwcom/EMEP/EMEP4UK/'+ \
                     'EMEP4UK_CountryFile.nc'
outf = nc.Dataset(EMEP4UK_countryfile,'w')

outf.createDimension('x',lats.shape[0])
outf.createDimension('y',lats.shape[1])

outvar=outf.createVariable('lat','float32',('x','y'))
outvar[:]=lats

outvar=outf.createVariable('lon','float32',('x','y'))
outvar[:]=lons

outvar=outf.createVariable('Country','float32',('x','y'))
outvar.note='Water=0, England=1, Scotland=2, Wales=3, Northern Ireland=4, Ireland=5, '\
             + 'France=6, Beglium=7, Netherlands=8, Denmark=9, Norway=10'

outvar[:]=COUNTRY_MASK

outf.title='Country Mask for the EMEP4UK 5km grid'
outf.owner='Edward Comyn-Platt, edwcom@ceh.ac.uk'

outf.close()






