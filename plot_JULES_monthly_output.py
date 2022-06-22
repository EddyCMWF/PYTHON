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
class plotJULESmonthly:
        def parse_input(self):
            #
            parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
            # optional
            parser.add_argument('--LATDEL',type=float,help='Latitude gridline spacing',required=False,default=2.5)
            parser.add_argument('--LONDEL',type=float,help='Longitude gridline spacing',required=False,default=2.5)
            #
            # positional
            parser.add_argument('infile',help='Input file')
            parser.add_argument('outdir',help='Output directory')
            #
            # Parse the arguments
            args=parser.parse_args()
            #
            return args.indir, args.outdir, args.LATDEL, args.LONDEL
        

###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    pJm=plotJULESmonthly()
    infile,outdir,LATDEL,LONDEL=pJm.parse_input()
    
    infile='/prj/wetlands_africa/jules/JASMIN/WFD_EI_global/MPI_WFD_EI_global.monthly.nc'
    outdir='/users/eow/edwcom/test_plots/'
    
    
    print 'Reading '+infile
    inf=nc.Dataset(infile,'r')
    lons=inf.variables['lon'][:]
    lats=inf.variables['lat'][:]
    
    
    
    inf.close()
    #
    #
    # Open infile2 and extract data
    print 'Reading '+infile2+':'
    inf2=nc.Dataset(infile2,'r')
    lons2=inf2.variables['lon'][:]
    lats2=inf2.variables['lat'][:]
    # Read top soil layer info
    if (bc):
            sathh2=inf2.variables['sathh'][0,:,:]
            b2=inf2.variables['b'][0,:,:]
    else:
            oneoveralpha2=inf2.variables['oneoveralpha'][0,:,:]
            oneovernminusone2=inf2.variables['oneovernminusone'][0,:,:]
    hcap2=inf2.variables['hcap'][0,:,:]
    hcon2=inf2.variables['hcon'][0,:,:]
    satcon2=inf2.variables['satcon'][0,:,:]
    vcrit2=inf2.variables['vcrit'][0,:,:]
    vsat2=inf2.variables['vsat'][0,:,:]
    vwilt2=inf2.variables['vwilt'][0,:,:]
    cs2=inf2.variables['cs'][0,:,:]
    inf2.close()
    
    CU_lonrange=[-15,15]
    CU_latrange=[50.,55.]
    #
    #Plot Hydrological Conductivity of saturation
    data_names = ['hcap',     'hcon','satcon','vcrit','vsat','vwilt','soilC' ]
    data1s     = [ hcap1*1e-6, hcon1, satcon1, vcrit1, vsat1, vwilt1, cs1    ] 
    data2s     = [ hcap2*1e-6, hcon2, satcon2, vcrit2, vsat2, vwilt2, cs2    ]
    NLEVELSs   = [ 16,         9,     15,      11,     13,    13,     8      ]
    #
    cbars      = ['RdYlBu_r','RdYlBu_r','RdYlBu_r','RdYlBu_r','RdYlBu_r','RdYlBu_r','YlGn_r']
    #
    data_ranges = [ [0.0,1.5],   \
                    [0.0,0.4],   \
                    [0.0,0.07], \
                    [0.0,0.5],   \
                    [0.3,0.8],   \
                    [0.0,0.3],  \
                    [0.0,35]   ]
    #
    cbar_titles = ['Heat Capacity ($J.m^{-3}.K^{-1}$ x$10^{6}$)',\
                   'Thermal Conductivity ($W.m^{-1}.K^{-1}$)',   \
                   'Hydraulic Conductivity ($kg.m^{-2}.s^{-1}$)',\
                   'Volumetric Water Content ($m^{3}.m^{3}$)',   \
                   'Volumetric Water Content ($m^{3}.m^{3}$)',   \
                   'Volumetric Water Content ($m^{3}.m^{3}$)',   \
                   'Carbon Content ($kg.m^{-2}$)'                ]
    #
    plot_titles = ['Heat Capacity of Dry Soil',                 \
                   'Thermal Conductivity of Dry Soil',          \
                   'Hydraulic Conductivity at Saturation',      \
                   'Volumetric Water Content at Critical Point',\
                   'Volumetric Water Content at Saturation',    \
                   'Volumetric Water Content at Wilting Point', \
                   'Soil Carbon Content to 1m Depth'            ]
    
    if (bc):
            data_names.append( 'sathh' )
            data_names.append( 'bexp' )
            data1s.append( sathh1 )
            data1s.append( b1 )
            data2s.append( sathh2 )
            data2s.append( b2 )
            NLEVELSs.append( 12 )
            NLEVELSs.append( 13 )
            cbars.append('RdYlBu_r')
            cbars.append('RdYlBu_r')
            data_ranges.append( [0.0,0.055] )
            data_ranges.append( [0.0,12] )
            cbar_titles.append( 'Hydraulic Conductivity ($m^{-1}$)' )
            cbar_titles.append( 'Brooks & Corey Exponent' )
            plot_titles.append( 'Hydraulic Conductivity at Saturation')
            plot_titles.append( 'Brooks & Corey Exponent' )
    else:
            data_names.append( 'oneoveralpha' )
            data_names.append( 'oneovernminusone' )
            data1s.append( oneoveralpha1 )
            data1s.append( oneovernminusone1 )
            data2s.append( oneoveralpha2 )
            data2s.append( oneovernminusone2 )
            NLEVELSs.append( 11 )
            NLEVELSs.append( 11 )
            cbars.append('RdYlBu_r')
            cbars.append('RdYlBu_r')
            data_ranges.append([0.0,1.0])
            data_ranges.append([0.0,10.0])
            cbar_titles.append('Reciprocal of van G alpha parameter',)
            cbar_titles.append('Reciprocal of van G n parameter minus one')
            plot_titles.append('Reciprocal of van G alpha parameter')
            plot_titles.append('Reciprocal of van G n parameter minus one')
    
            


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
            plotname=outdir+'EMEP_EUROPE_'+tag+'_soilparams_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=10., LONDEL=10., RESOLUTION='h', \
                        PROJECTION='stere')
            
            # 2. Plot  Entire EUROPE grid at 50 km with EMEP4UK 5km inset
            plotname=outdir+'EMEP4UK_EUROPE_'+tag+'_soilparams_'+data_name+'.png'
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
            plotname=outdir+'EMEP4UK_'+tag+'_soilparams_'+data_name+'.png'
            PM.plot_map(data1,lons1,lats1, \
                        DATA2=data2,LONS2=lons2,LATS2=lats2, \
                        DATA_RANGE=data_range, \
                        MAP_TYPE='Mesh', MPL_CBAR=cbar, NLEVELS=nlevels, CBAR_ORIENTATION='vertical', \
                        WIDTH=8, HEIGHT=8, PLOT_LABEL=cbar_title, \
                        PLOT_TITLE=plot_title, FONTSIZES=[12,12,12,18], \
                        iDISPLAY='N', FILE_PLOT=plotname, \
                        LATDEL=2., LONDEL=2., RESOLUTION='h', \
                        PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)
