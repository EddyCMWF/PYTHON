#!/bin/env python2.7

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon
import os
import datetime as dt
import PlotTools.plot_tools as PTs
from calendar import monthrange
from copy import deepcopy

import brewer2mpl

L_overwrite = True
plot_dir='/prj/GREENHOUSE/CHESS_runs/PERFECT-TREAT/'

JULES_dir='/prj/GREENHOUSE/CHESS_runs/jules_output/Jvn4.5/'
#JULES_filetag='chess_v1.0_ECP.monthly_vegcarb.YYYY.nc'
#JULES_vars=['npp_gb','gpp_gb','fch4_wetl'] #,'resp_s_gb','cs','resp_p_gb',\
   #            'fch4_wetl','smcl', 'tstar_gb']
JULES_filetag='chess_v1.0_ECP.monthly_soil.YYYY.nc'
JULES_vars=['fwetl']

JULES_CHESS_gridfile='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

CHESS_landcover_file='/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'

kgC_to_molesC = ( 1e3/12.01 )

start_year=2001
end_year=2009  #14
nYEARS=end_year-start_year+1
month_names=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
days_in_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
nMONTHs=12
hist_month_objects=[ dt.date(start_year,imonth+1,1) for imonth in range(nMONTHs) ]
line_month_objects=[ dt.date(start_year,imonth+1,15) for imonth in range(nMONTHs) ]

grinf=nc.Dataset(JULES_CHESS_gridfile,'r')
grindex=grinf.variables['index_2D'][:]
lats_2D=grinf.variables['lats_2D'][:]
lons_2D=grinf.variables['lons_2D'][:]
grinf.close()

LCinf=nc.Dataset(CHESS_landcover_file,'r')
LCfrac=LCinf.variables['frac'][:]
x=LCinf.variables['x'][:]
y=LCinf.variables['y'][:]
LCinf.close()

# Read in JULES data and convert units:
JULES_DATA={}
print('Reading Data')
for year in range(start_year,end_year+1):
    print(year)
    print(JULES_filetag.replace('YYYY',str(year)))
    Jinf=nc.Dataset(JULES_dir+JULES_filetag.replace('YYYY',str(year)),'r')
    temp_time=nc.num2date(Jinf.variables['time_bounds'][:,0].squeeze(), 
                          units=Jinf.variables['time'].units,
                          calendar=Jinf.variables['time'].calendar)
    if year==start_year:
        JULES_time=temp_time
        JULES_monthdays = np.array([ monthrange(dtime.year,dtime.month)[1] for dtime in temp_time])
    else:
        JULES_time=np.append(JULES_time,temp_time)
        JULES_monthdays = np.append(JULES_monthdays,
                                np.array([ monthrange(dtime.year,dtime.month)[1] for dtime in temp_time]) )

    for var in JULES_vars:
        invar = Jinf.variables[var]
        indata = invar[:].squeeze()[...,grindex]
        indata = indata*np.ones_like(grindex)
        indata.data[indata.mask==True]=indata.fill_value
        if year==start_year:
            JULES_DATA[var]={'data': indata.data, 
                             'units': invar.units, 
                             'long_name': invar.long_name, 
                             'fill_value':indata.fill_value }
        else:
            JULES_DATA[var]['data']=np.append(JULES_DATA[var]['data'],
                                              indata.data,axis=0)

    Jinf.close()   

nTSTEPS,nY,nX=JULES_DATA['fwetl']['data'].shape
print('Converting Data')
var='fwetl'
JULES_DATA[var]['long_name']='Fraction of ground water fed wetland'
JULES_DATA[var]['File_tag']='Wetland_Area'
JULES_DATA[var]['data']=np.ma.masked_equal(JULES_DATA[var]['data'],
                                                   JULES_DATA[var]['fill_value'] )
JULES_DATA[var]['units']='km $^{2}$'
JULES_DATA[var]['data_range']=[0,0.25]
JULES_DATA[var]['data_range_MN']=[0,0.25]
JULES_DATA[var]['hist_range']=[0,1e4]
JULES_DATA[var]['Annual-Mean']=np.mean(JULES_DATA[var]['data'].
                               reshape(end_year-start_year+1,-1,nY,nX),axis=1) 
JULES_DATA[var]['climatology']=np.mean(JULES_DATA[var]['data'].
                               reshape(end_year-start_year+1,-1,nY,nX),axis=0)
JULES_DATA[var]['Monthly-TS']=np.sum(JULES_DATA[var]['data'].
                                     reshape(nTSTEPS,-1),axis=1)
JULES_DATA[var]['Annual-TS']=np.sum(JULES_DATA[var]['Annual-Mean'].
                                 reshape(nYEARS,-1),axis=1) 

half_window=6
JULES_DATA[var]['Monthly-TS-rolling'] = np.array(
                                          [np.mean(JULES_DATA[var]['Monthly-TS'][i-half_window:i+half_window])
                                            for i in range(half_window,nTSTEPS-half_window)])

PLOT_FWETL=True
if PLOT_FWETL:
    print('plotting fwetl')
    var='fwetl'
    # Methane emission plots:
    os.system('mkdir -p '+plot_dir+var+'/')
    CMAP=brewer2mpl.get_map('YlGnBu','Sequential','9',reverse=True)
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=0.8)
    PLOT_DICT=JULES_DATA[var] 
    Ann_plot_units='km$^{2}$ km$^{-2}$'
    Mon_plot_units='km$^{2}$ km$^{-2}$'
    
    ANNUAL_FWETL=True
    if ANNUAL_FWETL:
      print('Annual Wetland Area plot:')
      for i in range(nYEARS):
        print(str(i+start_year))
        FILE_NAME = plot_dir+var+'/Annual_'+PLOT_DICT['File_tag']+str(i+start_year)+'.png'
        if ( (not os.path.isfile(FILE_NAME)) or L_overwrite):
            print(FILE_NAME)
            PLOT_TITLE='Annual Mean '+PLOT_DICT['long_name']+'\n'+str(i+start_year)
            PLOT_DATA=PLOT_DICT['Annual-Mean'][i,:]
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                    DATA_RANGE=PLOT_DICT['data_range'],
                    LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                    CMAP=MPL_CMAP,NLEVELS=100,NTICKS=6,TICK_FORMAT='%0.2f',
                    CBAR_LABEL=Ann_plot_units,
                    CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                    PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                    )

    MONTHLY_CH4=False
    if MONTHLY_CH4:
      print('Monthly emission plots:')
      for i in range(nTSTEPS):
        datestring=JULES_time[i].strftime('%Y-%m')
        print(datestring)
        FILE_NAME = plot_dir+var+'/'+PLOT_DICT['File_tag']+datestring+'.png'
        if ( (not os.path.isfile(FILE_NAME)) or L_overwrite):
            print(FILE_NAME)
            PLOT_TITLE='Month Mean '+PLOT_DICT['long_name']+'\n '+datestring
            PLOT_DATA=PLOT_DICT['data'][i,:]
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                    DATA_RANGE=PLOT_DICT['data_range_MN'],
                    LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                    CMAP=MPL_CMAP,NLEVELS=100,NTICKS=6,TICK_FORMAT='%0.2f',
                    CBAR_LABEL=Mon_plot_units,
                    CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                    PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                    )
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,4))
    ax2=ax.twinx()
    ax.bar(JULES_time,PLOT_DICT['Monthly-TS'],width=JULES_monthdays,
            color=CMAP.hex_colors[4],lw=1,edgecolor=CMAP.hex_colors[1])
    #ax2.plot(JULES_time[5::12],PLOT_DICT['Annual-TS'],
    #         color=CMAP.hex_colors[0],lw=3)
    ax2.plot(JULES_time[half_window:-half_window],PLOT_DICT['Monthly-TS-rolling'],
             color=CMAP.hex_colors[0],lw=3)
    ax.set_xlim([dt.date(start_year,1,1),dt.date(end_year+1,1,1)]) 
    #ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel(PLOT_DICT['units'],color=CMAP.hex_colors[4],
                    fontsize=16,fontweight='bold')
    ax2.set_ylabel(PLOT_DICT['units'],color=CMAP.hex_colors[0],
                    fontsize=16,fontweight='bold')
    ax2.set_ylim(PLOT_DICT['hist_range'])
    ax.set_title('Total '+PLOT_DICT['File_tag'].replace('_',' ')+' for Great Britain',
                    fontsize=20,fontweight='bold')
    fig.savefig(plot_dir+var+'/Histogram_'+str(start_year)+str(end_year)+'.png')
    plt.close()
   




