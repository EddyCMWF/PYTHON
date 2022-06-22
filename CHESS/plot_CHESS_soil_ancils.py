#!/bin/env python


import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import sys
from utils import io
import PlotTools.plot_tools as PT
import maths_tools.MathsTools as MT

infile=sys.argv[1]

outdir=io.optional_argparse('-outdir','./')
plot_vars=io.optional_argparse('-plot_vars',None)

colours=io.optional_argparse('-colours','khaki,goldenrod,darkred')
colours=colours.split(',')
print(colours)

if plot_vars==None:
    if 'vg' in infile:
        plot_vars=[ 'oneoveralpha','oneovernminusone','hcap','hcon','satcon',\
                     'vcrit','vsat','vwilt','cs'] 
        max_vals =[0.8,10,1.2,0.24,0.07,0.8,0.8,0.8,30]
        plot_tag='VG_'
    elif 'bc' in infile:
        plot_vars=[ 'sathh','b','hcap','hcon','satcon',\
                     'vcrit','vsat','vwilt','cs'] 
        max_vals =[0.05,11.,1.3,0.4,0.012,0.3,0.5,0.2,30]
        plot_tag='BC_'
    else:
        print('Not a BC or VG soil file, please define parameters')
else:
    plot_vars=plot_vars.split(',')
    max_vals =[None for var in plot_vars]
    plot_tag=io.optional_argparse('-plot_tag','')

lat_name=io.optional_argparse('-lat_name','lat')
lon_name=io.optional_argparse('-lon_name','lon')
z_dim=int(io.optional_argparse('-z_dim','0'))

lat_range=[48,60]
lon_range=[-8,4]

print(infile)
print(plot_vars)
print(lat_name,lon_name)

inf=nc.Dataset(infile,'r')
lats=inf.variables[lat_name][:]
lons=inf.variables[lon_name][:]

print(lats.shape, lons.shape)
for var,max_val in zip(plot_vars,max_vals):
    print(var)
    plot_data=inf.variables[var][:].squeeze()
    plot_title=inf.variables[var].long_name.replace('_',' ')
    cbar_title=inf.variables[var].units
    print(plot_title,cbar_title)
    if len(plot_data.shape)>2:
        plot_data=plot_data[0,:]
    
    if var=='hcap':
        plot_data/=1e6
        cbar_title='M'+cbar_title
    
    if max_val!=None:
        data_range=[0,max_val]
    else:
        data_range=[0,MT.round2SignifFigs(np.max(plot_data),2)]

    print(data_range)
    print(plot_data.shape)
    PT.plot_map(plot_data,lons,lats,RESOLUTION='i',PROJECTION='chess',
                DATA_RANGE=data_range,
                COLOURS=colours,INTERPOLATE_COLOURS=True,
                NLEVELS=100,iDISPLAY='N',NTICKS=6,TICK_FORMAT='%0.2f',
                FILE_PLOT=outdir+plot_tag+var+'.png',
                PLOT_TITLE=plot_title,FONTSIZES=[12,12,14,18],
                LATDEL=2,LONDEL=2,CBAR_LABEL=cbar_title,
                LON_RANGE=lon_range)
                #LAT_RANGE=lat_range,LON_RANGE=lon_range)


inf.close()


