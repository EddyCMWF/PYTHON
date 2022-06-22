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
class plotsoilmaps:
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
        # optional
        parser.add_argument('--indir',type=str,help='Input1 dir',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_drive/')
        parser.add_argument('--outdir',type=str,help='Output directory',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_drive/plots/' )
        parser.add_argument('--date_index',type=int,help='date index',required=False, \
                             default=4)
        #
        # positional
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.indir, args.outdir, args.date_index


###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    psm=plotsoilmaps()
    indir,outdir,date_index =psm.parse_input()
    #
    #indir='/users/eow/edwcom/EMEP/EMEP4UK/JULES_drive/'
    #indir='/prj/ALANIS/deposition/'
    #outdir='/users/eow/edwcom/EMEP/EMEP4UK/JULES_drive/plots/'
    #
    # First file for EUROPE as it is plotted underneath
    infile1=indir+'wrfout_d01_2001-07.nc'   ###'wrfout_d01_2001-07-01_00:00:00'  ###'wrfout_d01_2001-07.nc'   ###
    infile2=indir+'wrfout_d03_2001-07.nc'   ###'wrfout_d01_2001-07-01_00:00:00'  ###'wrfout_d01_2001-07.nc'   ###
    #
    #
    WRF_varnames=['SWDOWN','GLW','RAINNC','SNOWNC','U10','V10','PSFC','Q2','T2']
    #
    #
    # Open infile1 and extract data
    print 'Reading '+infile1+':'
    inf1=nc.Dataset(infile1,'r')
    lats1=inf1.variables['XLAT'][:]
    lons1=inf1.variables['XLONG'][:]
    data1s=[]
    longnames1=[]
    units1=[]
    # extract date_index from data 
    for var in WRF_varnames:
        data1s.append(inf1.variables[var][date_index,:,:])
        longnames1.append(inf1.variables[var].description)
        units1.append(inf1.variables[var].units)
    inf1.close()
    #
    #
    # Open infile2 and extract data
    print 'Reading '+infile2+':'
    inf2=nc.Dataset(infile2,'r')
    lats2=inf2.variables['XLAT'][:]
    lons2=inf2.variables['XLONG'][:]
    data2s=[]
    longnames2=[]
    units2=[]
    # extract date_index from data 
    for var in WRF_varnames:
        data2s.append(inf2.variables[var][date_index,:,:])
    inf2.close()
    #
    
    CU_lonrange=[-15,15]
    CU_latrange=[50.,55.]
    #
    #    
    NLEVELSs   = [ 11 for counter in range(len(WRF_varnames))]
    #
    cbars      = ['RdYlBu_r' for counter in range(len(WRF_varnames))]
    #
    data_ranges = [  [0.,1000.], \
                     [0.,500.], \
                     [0.0,0.05], \
                     [0.0,0.05], \
                     [-10.0,10.0], \
                     [-10.0,10.0], \
                     [90e3,105e3], \
                     [0.,0.01], \
                     [270.,310.]  ]
    #
    #
    units1[2] = 'mm s-1'
    #
    #
    for varnum in range(len(WRF_varnames)):
            data_name   = WRF_varnames[varnum]
            print 'plotting maps for '+data_name
            data1       = data1s[varnum]
            data2       = data2s[varnum]
            nlevels     = NLEVELSs[varnum]
            cbar        = cbars[varnum]
            #data_range  = [ min(np.amin(data1),np.amin(data2)), \
            #                max(np.amax(data1),np.amax(data2))  ]
            data_range = data_ranges[varnum]
            cbar_title = data_name+'  ('+units1[varnum]+')'
            plot_title = longnames1[varnum]
            
            print data_range
            
            # 1. Plot  Entire EUROPE grid at 50 km
            plotname=outdir+'WRF_d01_20010701_1200_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=10., LONDEL=10., RESOLUTION='h', \
                        PROJECTION='stere')
            
            # 2. Plot  Entire EUROPE grid at 50 km with EMEP4UK 5km inset
            plotname=outdir+'WRF_d01d03_20010701_1200_'+data_name+'.png'
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
            plotname=outdir+'WRF_d03_20010701_1200_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA2=data2,LONS2=lons2,LATS2=lats2, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=2., LONDEL=2., RESOLUTION='h', \
                        PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)


