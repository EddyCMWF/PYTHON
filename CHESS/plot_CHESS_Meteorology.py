#!/bin/env python
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import sys
from utils import io
import PlotTools.plot_tools as PT
import maths_tools.MathsTools as MT

indir=sys.argv[1]
if (indir[-1]!='/'): indir+='/'

outdir=io.optional_argparse('-outdir','./')
plot_tag=io.optional_argparse('-plot_tag','')
YYYYMM=io.optional_argparse('-YYYYMM','201206')

colours=io.optional_argparse('-colours','blue,white,red')
colours=colours.split(',')
print(colours)

plot_vars=io.optional_argparse('-plot_vars','dtr,precip,rlds,sfcWind,huss,psurf,rsds,tas')
plot_vars=plot_vars.split(',')
nVARs=len(plot_vars)

VARinf = { 'tas': {'range':[0,15],'tickformat':'%0.1f','method':'mean',\
                   'scalefactor':1,'offset':-273.15,'units':'$^o$C' }, \
           'dtr': {'range':[0,10],'tickformat':'%0.1f','method':'mean',\
                   'scalefactor':1,'offset':0,'units':'$^o$C' }, \
           'huss':{'range':[0,10],'tickformat':'%0.1f','method':'mean',\
                   'scalefactor':1e3,'offset':0,'units':'g kg$^{-1}$'}, \
           'precip':{'range':[0,10],'tickformat':'%0.12f','method':'mean',\
                   'scalefactor':86400.,'offset':0,'units':'mm day$^{-1}$'}, \
           'psurf':{'range':[900,1050],'tickformat':'%0.1f','method':'mean',\
                   'scalefactor':1e-2,'offset':0,'units':'hPa'}, \
           'rlds':{'range':[350,400],'tickformat':'%0.1f','method':'max',\
                   'scalefactor':1,'offset':0,'units':'W m$^{-2}$'}, \
           'rsds':{'range':[290,320],'tickformat':'%0.1f','method':'max',\
                   'scalefactor':1,'offset':0,'units':'W m$^{-2}$'}, \
           'sfcWind':{'range':[4,12],'tickformat':'%0.1f','method':'max',\
                      'scalefactor':1,'offset':0,'units':'m s$^{-1}$'}, \
        }

lat_name=io.optional_argparse('-lat_name','lat')
lon_name=io.optional_argparse('-lon_name','lon')

lat_range=[50,60]
lon_range=[-8,4]

print(indir,outdir)
print(plot_vars)
print(lat_name,lon_name)

#get lat and lon from first variable file
infile=indir+'chess_'+plot_vars[0]+'_'+YYYYMM+'.nc'
print(infile)
inf=nc.Dataset(infile,'r')
lats=inf.variables[lat_name][:]
lons=inf.variables[lon_name][:]
inf.close()

print(lats.shape, lons.shape)
#for var,max_val in zip(plot_vars,max_vals):
for var in plot_vars:
    print(var)
    infile=indir+'chess_'+var+'_'+YYYYMM+'.nc'
    print(infile)
    inf=nc.Dataset(infile,'r')
    plotvar=inf.variables[var]
    method=VARinf[var]['method']
    if method=='mean':
        plot_data=np.mean( plotvar[:],axis=plotvar.dimensions.index('time') )
    elif method=='max':
        plot_data=np.max( plotvar[:],axis=plotvar.dimensions.index('time') )
    
    plot_data=plot_data*VARinf[var]['scalefactor']
    plot_data=plot_data+VARinf[var]['offset']
    plot_title=plotvar.long_name.replace('_',' ')+' - '+YYYYMM+' - '+method
    cbar_title=plotvar.units
    filename=outdir+plot_tag+var+'.png'
    print(plot_title,cbar_title)
    print(filename)

    data_range=VARinf[var]['range']

    print(data_range)
    PT.plot_map(plot_data,lons,lats,RESOLUTION='i',PROJECTION='chess',
                DATA_RANGE=data_range,
                COLOURS=colours,INTERPOLATE_COLOURS=True,
                NLEVELS=100,iDISPLAY='N',NTICKS=6,TICK_FORMAT='%0.1f',
                FILE_PLOT=filename,
                PLOT_TITLE=plot_title,FONTSIZES=[12,12,12,14],
                LATDEL=2,LONDEL=2,CBAR_LABEL=cbar_title,
                LON_RANGE=lon_range)


inf.close()


