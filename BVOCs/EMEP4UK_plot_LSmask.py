#!/usr/bin/python
#
# Python module to plot the EMEP4UK LS mask
#
# Edard Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
#
#import data_info
#import gdal_functions
import plot_map_ECP as plotmap
#import plot_contour
#import plot_functions
#
# Input arguments:
if (len(sys.argv)>1):
    INFILE         = sys.argv[1]
else:
    INFILE         = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_gridfile.nc'

# INFILE = '/users/eow/edwcom/EMEP/EMEPgrid.nc'
#INTER          = sys.argv[2]
#MAP_TYPE       = sys.argv[3]
#iDISPLAY       = sys.argv[4]
#PLOT_OPT       = sys.argv[5]
#

inf=nc.Dataset(INFILE,'r')

grid_dims =  inf.variables['grid_dims'][:]

LSmask = inf.variables['grid_imask'][:]
LSmask = LSmask.reshape(grid_dims)

lon_cens = inf.variables['grid_center_lon'][:]
lon_cens = lon_cens.reshape(grid_dims)
lat_cens = inf.variables['grid_center_lat'][:]
lat_cens = lat_cens.reshape(grid_dims)

lon_cors = inf.variables['grid_corner_lon'][:]
lon_cors = lon_cors.reshape(np.append(grid_dims,4))
lat_cors = inf.variables['grid_corner_lat'][:]
lat_cors = lat_cors.reshape(np.append(grid_dims,4))


plotmap.plot_map(LSmask[0],lon_cens,lat_cens)

