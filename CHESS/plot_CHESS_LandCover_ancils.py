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
frac_name=io.optional_argparse('-frac_name','frac')
plot_tag=io.optional_argparse('-plot_tag','LC_')

colours=io.optional_argparse('-colours','indigo,gold,red')
colours=colours.split(',')
print(colours)

LCnames=io.optional_argparse('-LCnames','Broadleaf_Tree,Needleleaf_Tree,C3_Grass,Shrub,Crop,Urban,Lake,Bare_Soil').split(',')
nLCs=len(LCnames)

LCindex=io.optional_argparse('-LCindex',None)
if LCindex==None:
    LCindex=range(nLCs)
else:
    LCindex=[int(index) for index in LCindex.split(',')]

lat_name=io.optional_argparse('-lat_name','lat')
lon_name=io.optional_argparse('-lon_name','lon')

lat_range=[50,60]
lon_range=[-8,4]
minLCvalue=float(io.optional_argparse('-minLCvalue',0.01))

print(infile,outdir,frac_name)
print(LCnames)
print(lat_name,lon_name)

inf=nc.Dataset(infile,'r')
lats=inf.variables[lat_name][:]
lons=inf.variables[lon_name][:]

LCdata=inf.variables[frac_name][:]
inf.close()

# plot maximum landcover factor
maxLC=np.argmax(LCdata,axis=0)

maxLC=np.ma.masked_array(maxLC,mask=LCdata[0,:].mask)

LCcolours=['springgreen','darkgreen','yellowgreen','peru','khaki','red','blue','saddlebrown']
filename=outdir+plot_tag+'maxLC.png'
plot_title='Maximum Land Cover'
CBar_names=['BL Tree','NL Tree','C3 Grass','Shrub','Crop','Urban','Lake','Bare Soil']
PT.plot_map(maxLC,lons,lats,RESOLUTION='i',PROJECTION='chess',
            COLOURS=LCcolours,CLEVELS=[0,1,2,3,4,5,6,7,8],
            FILE_PLOT=filename,TickLEVELS=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5],
            TickLABELS=CBar_names,
            PLOT_TITLE=plot_title,FONTSIZES=[12,12,14,18],
            CBAR_ORIENTATION='vertical',CBAR_TICK_LENGTH=0,
            LATDEL=2.5,LONDEL=2.5,
            LAT_RANGE=lat_range,LON_RANGE=lon_range)

quit()


#plot cover of each tile
for index in LCindex:
    plot_data=np.ma.masked_less(LCdata[index,:].data,minLCvalue)
    plot_title=LCnames[index].replace('_',' ')
    
    plot_data=np.ma.masked_array(plot_data,mask=(plot_data.mask)|(plot_data.data<=0))

    cbar_title='Fraction of Grid Cell'
    print(plot_title,cbar_title)
    data_range=[0,1]
    filename=outdir+plot_tag+LCnames[index]+'.png'

    PT.plot_map(plot_data,lons,lats,RESOLUTION='i',PROJECTION='chess',
                DATA_RANGE=data_range,
                COLOURS=colours,INTERPOLATE_COLOURS=True,
                NLEVELS=100,iDISPLAY='N',NTICKS=6,TICK_FORMAT='%0.1f',
                FILE_PLOT=filename, 
                PLOT_TITLE=plot_title,FONTSIZES=[12,12,14,18],
                LATDEL=90,LONDEL=360,CBAR_LABEL=cbar_title,
                LON_RANGE=lon_range)





