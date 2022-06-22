#!/bin/python2.7

import csv
import netCDF4 as nc
import numpy as np
import sys
import EMEP_tools as EMEPts
import mpl_toolkits.basemap.pyproj as pyproj


def xy_to_lonlat(x,y,EMEP_proj='50'):
    Re    = 6370.         # radius of Earth
    lat_0 = 60.           # defining latitude
    lon_0 = -32.          # rotation angle
    #
    if (EMEP_proj == '50'):
        xpol=8.           # x coord of pole
        ypol=110.         # y coord of pole
        d=50.             # grid length at defining latitude
    elif (EMEP_proj == '150'):
        xpol=3.            # x coord of pole
        ypol=37.           # y coord of pole
        d=150.             # grid length at defining latitude
    else:
        raise Exception('Unrecognised EMEP projection')
    
    M = (Re/d) * ( 1. + np.sin( np.deg2rad(lat_0) ) )
    r = np.sqrt( ((x-xpol)**2.) + ((y-ypol)**2.) )
    
    lat = 90. - ( (360./np.pi) * np.arctan( r/M ) )
    lon = lon_0 + ( (180./np.pi) * np.arctan( (x-xpol)/(ypol-y) ) )
    
    # Correct for points below the pole coord where the rotation angle should be lon_0+180
    lon[y>ypol]   = (lon[y>ypol]+180.)
    # Correct for points which are gt 180., i.e. change to negative
    lon[ lon>(180.) ] = lon[ lon>(180.) ] - 360.
    # correct for pole point
    lon[ (x==xpol) & (y==ypol) ] = 0.
    
    return lon,lat


# Parse arguments, infile, outfile, geofile
infile=sys.argv[1]
outfile=sys.argv[2]

infile='/users/eow/edwcom/EMEP/EMEP_DepositionData2009v2015.csv'
outfile='/users/eow/edwcom/EMEP/EMEP_DepositionData2009v2015_WSG84.nc'

# Create arrays of output lons and lats



# Read infile as a csv file:
indata=list(csv.reader(open(infile,'r')))
# Get indata tile and headers from firt two lines
intitle=indata.pop(0)[0].replace('# ','')
inhdrs=indata.pop(0)
# remove the # from the first header
inhdrs[0]=inhdrs[0].replace('# ','')

# Create a dictionary to store this horrible data in
# Dictionary contains an empty list for each of the headers
indict = {}
for hdr in inhdrs:
    indict[hdr]=[]

# now put the data into the appropriate list in the dictionary
for dat in indata:
    for i in range(len(inhdrs)):
        indict[ inhdrs[i] ].append(dat[i])

# now we can put the usefule bits into np.arrays
data_all=np.array(indict['VALUE'],dtype='float32')
i_coord_all=np.array(indict['i'],dtype='float32')
j_coord_all=np.array(indict['j'],dtype='float32')
component_all=np.array(indict['COMPONENT'])

# create a list of the different components
unique_components=list(set(indict['COMPONENT']))


for comp in unique_components:
    index= component_all==comp
    
    data=data_all[index]
    i_coord=i_coord_all[index]
    j_coord=j_coord_all[index]
    
    # Now use my function defined at the top to convert EMEP coords into lat lon
    lon,lat=xy_to_lonlat(i_coord,j_coord)
    
    #We can now use these lon/lat pairs to create our output index points
    x = int(np.round(lon*2)+360)
    y = int(np.round(lat*2)+180)












    



