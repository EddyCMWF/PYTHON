#!/bin/env python3

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import os
import datetime as dt
import PlotTools.plot_tools as PTs

import brewer2mpl


plot_dir='/prj/GREENHOUSE/CHESS_runs/RoyalSocPres/WW_diffs/'

JULES_dir='/prj/GREENHOUSE/CHESS_runs/jules_output/'

std_subdir='albmar_Jvn4.5/'
std_filetag='chess_v1.1.monthly.YYYY.nc'

WW_subdir='WinterWheat_Jvn4.5/'
WW_filetag='chess_wheat.month.YYYY.nc'

comp_vars=['gpp_gb','fqw_gb','ftl_gb']

JULES_CHESS_gridfile='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

CHESS_landcover_file='/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'

kgC_to_molesC = ( 1e3/12.01 )

year=2014
month_names=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
days_in_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
nMONTHs=12
hist_month_objects=[ dt.date(year,imonth+1,1) for imonth in range(nMONTHs) ]
line_month_objects=[ dt.date(year,imonth+1,15) for imonth in range(nMONTHs) ]

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
STDinf=nc.Dataset(JULES_dir+std_subdir+std_filetag.replace('YYYY',str(year)),'r')
WWinf=nc.Dataset(JULES_dir+WW_subdir+WW_filetag.replace('YYYY',str(year)),'r')
std_DATA={}
WW_DATA={}
for var in comp_vars:
    invar = STDinf.variables[var]
    indata = invar[:].squeeze()[...,grindex]
    indata = indata*np.ones_like(grindex)
    indata.data[indata.mask==True]=indata.fill_value
    std_DATA[var]={'data': indata, \
           'units': invar.units, \
           'long_name': invar.long_name }
    invar = WWinf.variables[var]
    indata = invar[:].squeeze()[...,grindex]
    indata = indata*np.ones_like(grindex)
    indata.data[indata.mask==True]=indata.fill_value
    WW_DATA[var]={'data': indata, \
          'units': invar.units, \
          'long_name': invar.long_name }

STDinf.close()   
WWinf.close()   


WW_DATA['fqw_gb']['long_name']='Latent Heat Flux'
WW_DATA['fqw_gb']['data_range']=[0,10]
std_DATA['fqw_gb']['data']=np.sum(std_DATA['fqw_gb']['data'].transpose(1,2,0)\
                  *days_in_month,axis=2)/365.
std_DATA['fqw_gb']['long_name']='Latent Heat Flux'
std_DATA['fqw_gb']['data_range']=[0,10]

WW_DATA['ftl_gb']['long_name']='Sensible Heat Flux'
WW_DATA['ftl_gb']['data_range']=[0,10]
std_DATA['ftl_gb']['data']=np.sum(std_DATA['ftl_gb']['data'].transpose(1,2,0)\
                  *days_in_month,axis=2)/365.
std_DATA['ftl_gb']['long_name']='Sensible Heat Flux'
std_DATA['ftl_gb']['data_range']=[0,10]

WW_DATA['gpp_gb']['long_name']='Gross Primary Productivity'
WW_DATA['gpp_gb']['data_range']=[0,4]
# change units to Gg per year
WW_DATA['gpp_gb']['data']*= 86400.*365.*1e-6
std_DATA['gpp_gb']['data']=np.sum(std_DATA['gpp_gb']['data'].transpose(1,2,0)\
                  *days_in_month,axis=2)/365.
std_DATA['gpp_gb']['long_name']='Gross Primary Productivity'
std_DATA['gpp_gb']['data_range']=[0,4]

DIFF_DATA=WW_DATA.copy()

for var in DIFF_DATA:
    DIFF_DATA[var]['data']-=std_DATA[var]['data']

os.system('mkdir -p '+plot_dir)
CMAP=brewer2mpl.get_map('RdBr','Diverging','9')
MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=1.0)

for var in comp_vars:
    var_plotdir=plot_dir
#    os.system('mkdir -p '+var_plotdir)
    
    FILE_NAME=var_plotdir+'Annual_'+var+'_mass.png'
    PTs.plot_map(ANNUAL_FLUX_mass,lons_2D,lats_2D,
          DATA_RANGE=JULES_DATA[var]['data_range'],TICK_FORMAT='%.1f',
          LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
          CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Ann_plot_units_mass,
          CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
          PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
          )
    
    PLOT_TITLE='Annual Mean '+JULES_DATA[var]['long_name']
    FILE_NAME=var_plotdir+'Annual_'+var+'_mols.png'
    PTs.plot_map(ANNUAL_FLUX_mols,lons_2D,lats_2D,
             DATA_RANGE=JULES_DATA[var]['data_range'],TICK_FORMAT='%.1f',
             LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
             CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Ann_plot_units_mols,
             CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
             PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
             )
    
    for imonth in range(nMONTHs):
        PLOT_TITLE=month_names[imonth]+' Mean '+JULES_DATA[var]['long_name']
        FILE_NAME=var_plotdir+'%02i_'%(imonth+1)+month_names[imonth]+'_'+var+'_mols.png'
        PTs.plot_map(MONTH_FLUX_mols[imonth],lons_2D,lats_2D,
             DATA_RANGE=[-1.,3.],TICK_FORMAT='%.1f',
             LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
             CMAP=MPL_CMAP_month,NLEVELS=100,NTICKS=11,CBAR_LABEL=Mon_plot_units_mols,
             CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
             PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
             )
        PLOT_TITLE=month_names[imonth]+' Mean '+JULES_DATA[var]['long_name']
        FILE_NAME=var_plotdir+'%02i_'%(imonth+1)+month_names[imonth]+'_'+var+'_mass.png'
        PTs.plot_map(MONTH_FLUX_mass[imonth],lons_2D,lats_2D,
             DATA_RANGE=[-1.,5.],TICK_FORMAT='%.1f',
             LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
             CMAP=MPL_CMAP_month,NLEVELS=100,NTICKS=11,CBAR_LABEL=Mon_plot_units_mols,
             CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
             PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
             )
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    ax.bar(hist_month_objects,(DAILY_FLUX_byMONTH*days_in_month)+0.2,days_in_month-0.4,
        color=CMAP.hex_colors[-4],lw=1,edgecolor=CMAP.hex_colors[-1]) 
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel('Tg CO$_2$')
    ax.set_title('Total '+JULES_DATA[var]['long_name']+' Great Britain')
    fig.savefig(var_plotdir+'Histogram_True.png')
    plt.close()
    
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    ax.bar(hist_month_objects,(DAILY_FLUX_byMONTH*365./12.)+0.2,days_in_month-0.4,\
        color=CMAP.hex_colors[-4],lw=1,edgecolor=CMAP.hex_colors[-1]) 
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel('Tg CO$_2$')
    ax.set_title('Total '+JULES_DATA[var]['long_name']+' Great Britain')
    fig.savefig(var_plotdir+'Histogram_EqualMonth.png')
    plt.close()



PLOT_TSTAR=True
if PLOT_TSTAR:
    # Methane emission plots:
    os.system('mkdir -p '+plot_dir+'tstar/')
    CMAP=brewer2mpl.get_map('RdYlBu','Diverging','9',reverse=True)
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=0.55)
    
    #Unit strings for plots
    text_units='degC'
    plot_units='$^{o}$C'
    data_range=[-5,20]

    JULES_Tstar=JULES_DATA['tstar_gb']['data']-273.15
    JULES_Ann_Tstar=np.mean(JULES_Tstar,axis=0)
    JULES_Tstar_timeseries=np.mean(JULES_Tstar.reshape(nMONTHs,-1),axis=1)
    JULES_Tstar_timeseries_std=np.std(JULES_Tstar.reshape(nMONTHs,-1),axis=1)
    

    PLOT_TITLE='Mean Surface Temperature'
    FILE_NAME=plot_dir+'tstar/Annual_MeanSurfaceTemperature.png'
    PTs.plot_map(JULES_Ann_Tstar,lons_2D,lats_2D,
         DATA_RANGE=data_range,TICK_FORMAT='%0.1f',
         LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
         CMAP=MPL_CMAP,NLEVELS=100,NTICKS=6,CBAR_LABEL=plot_units,
         CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
         PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
         )
    
    for imonth in range(nMONTHs):
    PLOT_TITLE=month_names[imonth]+' Mean Surface Temperature'
    FILE_NAME=plot_dir+'tstar/%02i_'%(imonth+1)+month_names[imonth]+\
            '_MeanSurfaceTemperature.png'
    PTs.plot_map(JULES_Tstar[imonth,:],lons_2D,lats_2D,
              DATA_RANGE=data_range,TICK_FORMAT='%0.1f',
              LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
              CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=plot_units,
              CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
              PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
              )
    

    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    polygon_x=line_month_objects+line_month_objects[::-1]
    polygon_y=np.append(JULES_Tstar_timeseries-JULES_Tstar_timeseries_std,\
            JULES_Tstar_timeseries[::-1]+JULES_Tstar_timeseries_std)
    polygon_y2=np.append(JULES_Tstar_timeseries-(JULES_Tstar_timeseries_std*2),\
             JULES_Tstar_timeseries[::-1]+(JULES_Tstar_timeseries_std*2))
    ax.fill(polygon_x,polygon_y2,c=CMAP.hex_colors[3])
    ax.fill(polygon_x,polygon_y,c=CMAP.hex_colors[2])
    ax.plot(line_month_objects,JULES_Tstar_timeseries,\
        color=CMAP.hex_colors[-1],lw=1)
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel(plot_units)
    ax.set_title('Mean Temperature for Great Britain')
    fig.savefig(plot_dir+'tstar/TimeSeries.png')
    plt.close()
    
  
