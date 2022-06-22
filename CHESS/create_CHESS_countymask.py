#!/bin/env python

import netCDF4 as nc
import numpy as np
from osgeo import gdal, gdalnumeric, ogr, osr
#import Image, ImageDraw
import os, sys
gdal.UseExceptions()
#from shapely.geometry import Point
#from fiona import collection

CHESS_DIR='/users/eow/edwcom/CHESS/'
CHESS_latlon_file=CHESS_DIR+'chess_jules_land_index.nc'
shapefile=CHESS_DIR+'uk_admin_map_shapefile'
outfile=CHESS_DIR+'CHESS_countymask.nc'

#read in lat and lon data
inf=nc.Dataset(CHESS_latlon_file,'r')
lats=inf.variables['lats_2D'][:]
lons=inf.variables['lons_2D'][:]
inf.close()
ny,nx=lats.shape

# Get the driver
driver = ogr.GetDriverByName('ESRI Shapefile')
# Open a shapefile
dataset = driver.Open(shapefile, 0)
layer=dataset.GetLayer()


geo_ref = layer.GetSpatialRef()
point_ref=ogr.osr.SpatialReference()
point_ref.ImportFromEPSG(4326)
ctran=ogr.osr.CoordinateTransformation(point_ref,geo_ref)

COUNTY_MASK=np.zeros_like(lats,dtype='int')

for i in range(nx):
    for j in range(ny):
        coordPoint = ctran.TransformPoint(float(lons[j,i]),float(lats[j,i]))
        pt=ogr.Geometry(ogr.wkbPoint)
        pt.SetPoint_2D(int(coordPoint[2]),coordPoint[0],coordPoint[1])
        for index in xrange(layer.GetFeatureCount()):
            feat=layer.GetFeature(index)
            georef=feat.GetGeometryRef()
#            print  georef.Contains(pt)
            if georef.Contains(pt):
                COUNTY_MASK[j,i]=int(index)


outf=nc.Dataset(outfile,'w')
outf.createDimension('x',nx)
outf.createDimension('y',ny)
outvar=outf.createVariable('County_Mask','int',('y','x'),fill_value=-1)
outvar[:]=COUNTY_MASK
outf.note='CHESS grid mask for UK counties based on uk_admin_map_shapefile'
outf.close()

