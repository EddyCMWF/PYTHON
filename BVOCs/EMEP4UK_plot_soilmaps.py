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
            parser.add_argument('--plot_tag',type=str,help='Tag for plot naming convention',required=False,default='')
            parser.add_argument('--LATDEL',type=float,help='Latitude gridline spacing',required=False,default=2.5)
            parser.add_argument('--LONDEL',type=float,help='Longitude gridline spacing',required=False,default=2.5)
            parser.add_argument('--RIVERS',type=bool,help='Boolean flag to plot rivers',required=False,default=False)
            #
            # positional
            parser.add_argument('infile',help='Input file')
            parser.add_argument('outdir',help='Output directory')
            #
            # Parse the arguments
            args=parser.parse_args()
            #
            return args.infile, args.outdir, args.plot_tag, args.LATDEL, args.LONDEL, args.RIVERS
        

###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    psm=plotsoilmaps()
    infile,outdir,plot_tag,LATDEL,LONDEL,RIVERS=psm.parse_input()
    #
    # Open infile and extract data
    inf=nc.Dataset(infile,'r')
    lons=inf.variables['lon'][:]
    lats=inf.variables['lat'][:]
    # Read top soil layer info
    sathh=inf.variables['sathh'][0,:,:]
    b=inf.variables['b'][0,:,:]
    hcap=inf.variables['hcap'][0,:,:]
    hcon=inf.variables['hcon'][0,:,:]
    satcon=inf.variables['satcon'][0,:,:]
    vcrit=inf.variables['vcrit'][0,:,:]
    vsat=inf.variables['vsat'][0,:,:]
    vwilt=inf.variables['vwilt'][0,:,:]
    cs=inf.variables['cs'][0,:,:]
    
    #Plot Hydrological Conductivity of saturation
    PM.plot_map(sathh,lons,lats, \
                DATA_RANGE=[0.0,0.055], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=12, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Hydraulic Conductivity ($m^{-1}$)',\
                PLOT_TITLE='Hydraulic Conductivity at Saturation', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'sathh.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(b,lons,lats, \
                DATA_RANGE=[0.0,12], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=13, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Brooks & Corey Exponent',\
                PLOT_TITLE='Brooks & Corey Exponent', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'bexp.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(hcap*1e-6,lons,lats, \
                DATA_RANGE=[1.0,1.5], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=11, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Heat Capacity ($J.m^{-3}.K^{-1}$ x$10^{6}$)',\
                PLOT_TITLE='Heat Capacity of Dry Soil', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'hcap.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(hcon,lons,lats, \
                DATA_RANGE=[0.0,0.5], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=11, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Thermal Conductivity ($W.m^{-1}.K^{-1}$)',\
                PLOT_TITLE='Thermal Conductivity of Dry Soil', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'hcon.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(satcon,lons,lats, \
                DATA_RANGE=[0.0,0.014], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=8, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Hydraulic Conductivity ($kg.m^{-2}.s^{-1}$)',\
                PLOT_TITLE='Hydraulic Conductivity at Saturation', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'satcon.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(vcrit,lons,lats, \
                DATA_RANGE=[0.0,0.5], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=11, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Volumetric Water Content ($m^{3}.m^{3}$)',\
                PLOT_TITLE='Volumetric Water Content at Critical Point', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'vcrit.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(vsat,lons,lats, \
                DATA_RANGE=[0.3,0.6], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=7, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Volumetric Water Content ($m^{3}.m^{3}$)',\
                PLOT_TITLE='Volumetric Water Content at Saturation', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'vsat.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)

    PM.plot_map(vwilt,lons,lats, \
                DATA_RANGE=[0.0,0.25], \
                MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=6, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Volumetric Water Content ($m^{3}.m^{3}$)',\
                PLOT_TITLE='Volumetric Water Content at Wilting Point', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'vwilt.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)
    
    PM.plot_map(cs,lons,lats, \
                DATA_RANGE=[0.0,0.25], \
                MAP_TYPE='Mesh', MPL_CBAR='YlGn_r', NLEVELS=6, \
                WIDTH=12, HEIGHT=8, PLOT_LABEL='Carbon Content ($kg.m^{-2}$)',\
                PLOT_TITLE='Soil Carbon Content to 1m Depth', FONTSIZES=[12,12,12,18], \
                iDISPLAY='N', FILE_PLOT=outdir+plot_tag+'vwilt.png', \
                LATDEL=LATDEL, LONDEL=LONDEL, RESOLUTION='h', RIVERS=RIVERS)
