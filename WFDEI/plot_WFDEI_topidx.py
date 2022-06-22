#!/usr/bin/python
#
# Python module to plot differnces between topographic index data (initially on NCEP-CRU grid)
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import sys
import numpy as np
#import argparse
import netCDF4 as nc
#import netcdftime as ncdt
#import matplotlib.pyplot as plt
import plot_tools as PT
#
#

WFDEI_dir = '/users/eow/edwcom/WFD_EI/'
PLOT_DIR=WFDEI_dir+'plots/'
grid_file = 'wfdei-land-mask.nc'

if '-hydro1k' in sys.argv:
    infile='hydro1k_topo_watch_0p5deg.nc' 
    plot_tag='hyrdro1k'
else:
    infile='topoidx_WFDEI_0p5_lp_global_filled.nc'
    plot_tag='Marthews'

vars=['ti_mean','ti_std']
nvars=len(vars)
var_titles=[plot_tag+' Topographic Index Mean', plot_tag+' Topographic Index Standard Deviation']
var_cbar_titles=['Mean TI', 'Standard Deviation TI']
var_dataranges = [ (0,10), (1,2.5) ]

var_colours=[ ('aliceblue','powderblue','lightseagreen','midnightblue') , \
              ('honeydew','lightgreen','darkseagreen','olive')   ]
fill_value=-999.

# read in gridfile
grinf=nc.Dataset(WFDEI_dir+grid_file,'r')
grindex=grinf.variables['land_index'][:]-1
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
lons_2d,lats_2d=np.meshgrid(lons,lats)

# read in data
inf=nc.Dataset(WFDEI_dir+infile,'r')
for ivar in range(nvars):
    indata=inf.variables[vars[ivar]][:]
    indata=np.ma.masked_array(indata[grindex],mask=grindex.mask,fill_value=fill_value)
    indata.data[indata.mask==True]=fill_value
    PT.plot_map(indata,lons_2d,lats_2d, \
                PLOT_TITLE=var_titles[ivar], \
                CBAR_LABEL=var_cbar_titles[ivar], NTICKS=11, \
                COLOURS=var_colours[ivar],INTERPOLATE_COLOURS=True,NLEVELS=100, \
                DATA_RANGE=var_dataranges[ivar], \
                FILE_PLOT=PLOT_DIR+plot_tag+'_'+vars[ivar]+'.png', \
                FONTSIZES=[12,12,14,18], \
                )

inf.close()

