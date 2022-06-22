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

JOEY_dir = '/users/eow/edwcom/JOEY/'
PLOT_DIR = JOEY_dir+'plots/'

infile='julesb.nc'
plot_tag='joey_dynamic_run'

vars=['carb']
nvars=len(vars)
var_titles=['Soil Carbon'] 
var_cbar_titles=['$kg m^-2$']
var_dataranges = [ (0,30) ]

var_colours=[ ('white','khaki','darkgreen')  ]

# read in data
inf=nc.Dataset(JOEY_dir+infile,'r')
lats=inf.variables['lat'][:]-0.5
lons=inf.variables['lon'][:]
lons[lons>=180]-=360.
lons=np.append(lons[360:],lons[:360])
print lons
lons_2d,lats_2d=np.meshgrid(lons,lats)
for ivar in range(nvars):
    indata=inf.variables[vars[ivar]][:]
    indata=np.sum(indata,axis=1)
    indata=np.mean(indata,axis=0)
    plotdata=np.append(indata[:,360:],indata[:,:360],axis=1)
    PT.plot_map(plotdata,lons_2d,lats_2d, \
                PLOT_TITLE=var_titles[ivar], \
                CBAR_LABEL=var_cbar_titles[ivar], \
                MPL_CBAR='YlGn', NLEVELS=11, \
                DATA_RANGE=var_dataranges[ivar], \
                FILE_PLOT=PLOT_DIR+plot_tag+'_'+vars[ivar]+'.png', \
                FONTSIZES=[12,12,14,18], \
                )

inf.close()

