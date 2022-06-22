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
#import os, sys
#import numpy as np
#import argparse
import netCDF4 as nc
#import netcdftime as ncdt
import matplotlib.pyplot as plot
#import plot_tools as PT
#
#
hydro_file='/users/eow/edwcom/CRUNCEP/hydro1k_topo_cruncep_0p5deg.nc'
marthews_file='/users/eow/edwcom/CRUNCEP/topoidx_CRUNCEP_0p5_lp_global.nc'

inf=nc.Dataset(hydro_file,'r')
lons=inf.variables['longitude'][:]
lats=inf.variables['latitude'][:]
ti_mean_H=inf.variables['timean'][:]
ti_std_H=inf.variables['tisd'][:]
inf.close()

inf=nc.Dataset(marthews_file,'r')
ti_mean_M=inf.variables['ti_mean'][:]
ti_std_M=inf.variables['ti_std'][:]
inf.close()

plot.subplot(1,2,1)
plot.plot(ti_mean_H,ti_mean_M,'r.')
plot.plot([0,20],[0,20],'k-')
plot.title('Mean Topographic Index')
plot.ylabel('Marthews TI mean')
plot.xlabel('Hydro 1k TI mean')
plot.xlim([0,20])
plot.ylim([0,15])


plot.subplot(1,2,2)
plot.plot(ti_std_H,ti_std_M,'b.')
plot.plot([0,20],[0,20],'k-')
plot.title('Std. Dev. of Topographic Index')
plot.ylabel('Marthews TI std')
plot.xlabel('Hydro 1k TI std')
plot.xlim([0,20])
plot.ylim([0,10])

plot.show()




