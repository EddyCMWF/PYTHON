#!/usr/bin/python
#
# Python module to plot the hwsd dat on EMEP grid
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import argparse
import netCDF4 as nc
import matplotlib.pyplot as plot
from mpl_toolkits.basemap import Basemap
import plot_map_ECP as PM
#
###################################################################################################################
# Define class
###################################################################################################################
class plotsoilmaps:
        def parse_input(self):
            #
            parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
            # optional
            parser.add_argument('--indir',type=str,help='Input1 dir',required=False, \
                                 default='/users/eow/edwcom/EMEP/hwsd2emep/input_data/')
            parser.add_argument('outdir',type=str,help='Output directory',required=False, \
                                 default='/users/eow/edwcom/EMEP/hwsd2emep/plots/' )
            #
            # positional
            #
            # Parse the arguments
            args=parser.parse_args()
            #
            return args.indir, args.outdir
        

###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    psm=plotsoilmaps()
    indir,outdir=psm.parse_input()
    
    indir='/users/eow/edwcom/EMEP/hwsd2emep/input_data/'
    outdir='/users/eow/edwcom/EMEP/hwsd2emep/plots/'
    
    # First file for EUROPE as it is plotted underneath
    infile1=indir+'hwsdREBIN_emep4ukEUROPE_grid.nc'
    infile2=indir+'hwsd_emep4uk_grid.nc'
    
    #
    # Open infile1 and extract data
    print 'Reading '+infile1+':'
    inf1=nc.Dataset(infile1,'r')
    lons1=inf1.variables['lon'][:]
    lats1=inf1.variables['lat'][:]
    hwsd1=inf1.variables['mu'][0,:,:]
    inf1.close()
    #
    #
    # Open infile2 and extract data
    print 'Reading '+infile2+':'
    inf2=nc.Dataset(infile2,'r')
    lons2=inf2.variables['lon'][:]
    lats2=inf2.variables['lat'][:]
    hwsd2=inf2.variables['mu'][0,:,:]
    inf2.close()
    
    CU_lonrange=[-15,15]
    CU_latrange=[50.,55.]
    #
    #Plot Hydrological Conductivity of saturation
    data_name = 'HWSD'
    data1     = hwsd1   #np.floor(hwsd1/100.)   #
    data2     = hwsd2   #np.floor(hwsd2/100.)   #
    nlevels   = 254
    #
    cbar      = 'Paired'
    cbar2     = 'Paired'
    #
    data_range1 = [0.,20000.]
    data_range2 = [10000.,11000.]
    TickLEVELS2 = [9000., 9500., 10000., 10500., 11000.]
    #
    cbar_title = 'HWSD Class'
    #
    plot_title = 'Harmonised World Soil Database Classes'
    
    # 1. Plot  Entire EUROPE grid at 50 km
    plotname=outdir+'EMEP_EUROPE_'+data_name+'.png'
    PM.plot_map(data1,lons1,lats1, \
                DATA_RANGE=data_range1, \
                MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                TickLEVELS=data_range1, \
                WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=plotname, \
                LATDEL=10., LONDEL=10., RESOLUTION='h', \
                PROJECTION='stere')
    
    # 2. Plot  Entire EUROPE grid at 50 km with EMEP4UK 5km inset
    plotname=outdir+'EMEP4UK_EUROPE_'+data_name+'.png'
    PM.plot_map(data1,lons1,lats1, \
                DATA2=data2,LONS2=lons2,LATS2=lats2, \
                DATA_RANGE=data_range1, \
                MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                TickLEVELS=data_range1, \
                WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=plotname, \
                LATDEL=10., LONDEL=10., RESOLUTION='h', \
                PROJECTION='stere')
    
    # 3. Plot close up EMEP4UK 5km data with EMEP 50km surround
    plotname=outdir+'EMEP4UK_'+data_name+'.png'
    PM.plot_map(data1,lons1,lats1, \
                DATA2=data2,LONS2=lons2,LATS2=lats2, \
                DATA_RANGE=data_range2, \
                MAP_TYPE='Mesh', MPL_CBAR=cbar2, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                TickLEVELS=data_range2, \
                WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=plotname, \
                LATDEL=2., LONDEL=2., RESOLUTION='h', \
                PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)


