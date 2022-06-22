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
import plot_map_ECP as PM
#
###################################################################################################################
# Define class
###################################################################################################################
class plottopomaps:
        def parse_input(self):
            #
            parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
            # optional
            parser.add_argument('--RIVERS',type=bool,help='Boolean flag to plot rivers',required=False,default=False)
            #
            # positional
            parser.add_argument('indir',help='Input1 dir')
            parser.add_argument('outdir',help='Output directory')
            #
            # Parse the arguments
            args=parser.parse_args()
            #
            return args.indir, args.outdir,args.RIVERS
        

###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    ptm=plottopomaps()
    indir,outdir,RIVERS=ptm.parse_input()
    
    infile1=indir+'topidx_emep4ukEUROPE_grid.nc'
    infile2=indir+'topidx_emep4uk_grid.nc'
    #
    # Open infile1 and extract data
    print 'Reading '+infile1+':'
    inf1=nc.Dataset(infile1,'r')
    lons1=inf1.variables['cen_lon'][:]
    lats1=inf1.variables['cen_lat'][:]
    timean1=inf1.variables['timean'][:]
    tistd1=inf1.variables['tistd'][:]
    inf1.close()
    # Open infile2 and extract data
    print 'Reading '+infile2+':'
    inf2=nc.Dataset(infile2,'r')
    lons2=inf2.variables['cen_lon'][:]
    lats2=inf2.variables['cen_lat'][:]
    timean2=inf2.variables['timean'][:]
    tistd2=inf2.variables['tistd'][:]
    inf2.close()
    
    # Close up lat/lon ranges
    CU_lonrange=[-15,15]
    CU_latrange=[50.,55.]
    
    data_names  = ['timean','tistd']
    data1s      = [ timean1, tistd1]
    data2s      = [ timean2, tistd2]
    NLEVELSs    = [    9   ,   9   ]
    data_ranges = [  [0,8] , [0,4] ]
    cbars       = [ 'jet'  , 'jet' ]
    
    cbar_titles = [ 'topographic index mean' , \
                    'topographic index standard deviaiton' ]   
    
    plot_titles = [ 'Topographic Index Mean' , \
                    'Topographic Index Standard Deviaiton' ]
    
    
    for varnum in range(len(data_names)):
            data_name   = data_names[varnum]
            print 'plotting maps for '+data_name
            data1       = data1s[varnum]
            data2       = data2s[varnum]
            nlevels     = NLEVELSs[varnum]
            cbar        = cbars[varnum]
            data_range  = data_ranges[varnum]
            cbar_title = cbar_titles[varnum]
            plot_title = plot_titles[varnum]
            
            # 1. Plot  Entire EUROPE grid at 50 km
            plotname=outdir+'EMEP_EUROPE_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=10., LONDEL=10., RESOLUTION='h', \
                        PROJECTION='stere')
            
            # 2. Plot  Entire EUROPE grid at 50 km with EMEP4UK 5km inset
            plotname=outdir+'EMEP4UK_EUROPE_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA2=data2,LONS2=lons2,LATS2=lats2, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=10., LONDEL=10., RESOLUTION='h', \
                        PROJECTION='stere')
            
            # 3. Plot close up EMEP4UK 5km data with EMEP 50km surround
            plotname=outdir+'EMEP4UK_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA2=data2,LONS2=lons2,LATS2=lats2, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=2., LONDEL=2., RESOLUTION='h', \
                        PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)
            
