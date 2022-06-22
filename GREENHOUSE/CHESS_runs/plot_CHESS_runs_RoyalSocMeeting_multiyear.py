#!/bin/env python3

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import os
import datetime as dt
import PlotTools.plot_tools as PTs

import brewer2mpl


plot_dir='/prj/GREENHOUSE/CHESS_runs/RoyalSocPres/'

JULES_dir='/prj/GREENHOUSE/CHESS_runs/jules_output/Jvn4.5/'
JULES_filetag='chess_v1.0_ECP.monthly_vegcarb.YYYY.nc'
JULES_vars=['npp_gb','gpp_gb','fch4_wetl'] #,'resp_s_gb','cs','resp_p_gb',\
#            'fch4_wetl','smcl', 'tstar_gb']

JULES_CHESS_gridfile='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

CHESS_landcover_file='/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'

kgC_to_molesC = ( 1e3/12.01 )

start_year=2001
end_year=2014
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
    Jinf=nc.Dataset(JULES_dir+JULES_filetag.replace('YYYY',str(year)),'r')
    temp_time=nc.num2date(Jinf.variables['time_bounds'][:,0].squeeze(), 
                          units=Jinf.variables['time'].units )
    if year==start_year:
        JULES_time=temp_time
    else:
        JULES_time=np.append(JULES_time,temp_time)
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

nTSTEPS,nY,nX=JULES_DATA['fch4_wetl']['data'].shape
print('Converting Data')
var='fch4_wetl'
JULES_DATA[var]['long_name']='Methane flux from natural wetlands'
JULES_DATA[var]['data']=np.ma.masked_equal(JULES_DATA[var]['data'],
                                                   JULES_DATA[var]['fill_value'] )
JULES_DATA[var]['data']*=1e-9 
JULES_DATA[var]['units']='kg m$^{-2}$ s$^{-1}$'
JULES_DATA[var]['data_range']=[0,2000]
JULES_DATA[var]['data_range_MN']=[0,200]
JULES_DATA[var]['hist_range']=[0,0.11]
JULES_DATA[var]['Annual-Mean']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=1) 
JULES_DATA[var]['climatology']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=0)
JULES_DATA[var]['Monthly-TS']=np.sum(JULES_DATA[var]['data'].\
                                     reshape(nTSTEPS,-1),axis=1)  \
                                 *86400.*(365./12.)*1e6
JULES_DATA[var]['Annual-TS']=np.sum(JULES_DATA[var]['Annual-Mean'].
                                    reshape(nYEARS,-1),axis=1)  \
                                 *86400.*365.*1e6

var='npp_gb'
JULES_DATA[var]['data']=np.ma.masked_equal(JULES_DATA[var]['data'],
                                           JULES_DATA[var]['fill_value'] )
JULES_DATA[var]['long_name']='Net Primary Productivity'
JULES_DATA[var]['data_range']=[0,3]
JULES_DATA[var]['data_range_MN']=[0,5]
JULES_DATA[var]['hist_range']=[0,200]
JULES_DATA[var]['Annual-Mean']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=1) 
JULES_DATA[var]['climatology']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=0) 
JULES_DATA[var]['Monthly-TS']=np.sum(JULES_DATA[var]['data'].\
                                     reshape(nTSTEPS,-1),axis=1)  \
                                 *86400.*(365./12.)*1e6
JULES_DATA[var]['Annual-TS']=np.sum(JULES_DATA[var]['Annual-Mean'].
                                    reshape(nYEARS,-1),axis=1)  \
                                 *86400.*365.*1e6

var='gpp_gb'
JULES_DATA[var]['data']=np.ma.masked_equal(JULES_DATA[var]['data'],
                                           JULES_DATA[var]['fill_value'] )
JULES_DATA[var]['long_name']='Gross Primary Productivity'
JULES_DATA[var]['data_range']=[1,4]
JULES_DATA[var]['data_range_MN']=[0,6]
JULES_DATA[var]['hist_range']=[0,300]
JULES_DATA[var]['data_range_monthly']=[2,4]
JULES_DATA[var]['Annual-Mean']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=1) 
JULES_DATA[var]['climatology']=np.mean(JULES_DATA[var]['data'].\
                               reshape(end_year-start_year+1,12,nY,nX),axis=0) 
JULES_DATA[var]['Monthly-TS']=np.sum(JULES_DATA[var]['data'].\
                                     reshape(nTSTEPS,-1),axis=1)  \
                                 *86400.*(365./12.)*1e6
JULES_DATA[var]['Annual-TS']=np.sum(JULES_DATA[var]['Annual-Mean'].
                                    reshape(nYEARS,-1),axis=1)  \
                                 *86400.*365.*1e6


PLOT_FCH4=True
if PLOT_FCH4:
    print('plotting fch4_wetl')
    # Methane emission plots:
    os.system('mkdir -p '+plot_dir+'fch4_wetl/')
    CMAP=brewer2mpl.get_map('YlGnBu','Sequential','9',reverse=True)
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=0.8)
    FCH4_DICT=JULES_DATA['fch4_wetl'] 
    #Unit strings for plots
    #text_units='GgCH4'
    #Ann_plot_units_mols='mmolCH4 m$^{-2}$ s$^{-1}$'
    Ann_plot_units_mass='kgCH4 km$^{-2}$ yr$^{-1}$'
    #Mon_plot_units_mols='mmolCH4 m$^{-2}$ s$^{-1}$'
    Mon_plot_units_mass='kgCH4 km$^{-2}$ mo$^{-1}$'
    
    ANNUAL_CH4=False
    if ANNUAL_CH4:
      print('Annual emission plotss:')
      for i in range(nYEARS):
        print(str(i+start_year))
        FILE_NAME = plot_dir+'fch4_wetl/Annual_Wetland_Emission_mass_'+\
                            str(i+start_year)+'.png'
        PLOT_TITLE='Annual Mean '+FCH4_DICT['long_name']+'\n'+str(i+start_year)
        PLOT_DATA=FCH4_DICT['Annual-Mean'][i,:]*(16.04/12.01)*1e6*(86400.*365.)
        PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                     DATA_RANGE=FCH4_DICT['data_range'],TICK_FORMAT='%0.1f',
                     LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                     CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                     CBAR_LABEL=Ann_plot_units_mass,
                     CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                     PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                     )

    MONTHLY_CH4=False
    if MONTHLY_CH4:
      print('Monthly emission plots:')
      for i in range(nTSTEPS):
        datestring=JULES_time[i].strftime('%Y-%m')
        print(datestring)
        FILE_NAME = plot_dir+'fch4_wetl/Wetland_Emission_mass_'+\
                            datestring+'.png'
        PLOT_TITLE='Annual Mean '+FCH4_DICT['long_name']+'\n '+datestring
        PLOT_DATA=FCH4_DICT['data'][i,:]*(16.04/12.01)*1e6*(86400.*365./12)
        PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                     DATA_RANGE=FCH4_DICT['data_range_MN'],TICK_FORMAT='%0.1f',
                     LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                     CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                     CBAR_LABEL=Mon_plot_units_mass,
                     CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                     PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                     )
    

    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    ax2=ax.twinx()
    ax.bar(JULES_time,FCH4_DICT['Monthly-TS']*1e-6,365./12.,\
            color=CMAP.hex_colors[4],lw=1,edgecolor=CMAP.hex_colors[1])
    ax2.plot(JULES_time[5::12],FCH4_DICT['Annual-TS']*1e-9,\
             color=CMAP.hex_colors[0],lw=3)
    ax.set_xlim([dt.date(start_year,1,1),dt.date(end_year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel('GgCH$_4$ per Month',color=CMAP.hex_colors[4],
                    fontsize=16,fontweight='bold')
    ax2.set_ylabel('GgCH$_4$ per Annum',color=CMAP.hex_colors[0],
                    fontsize=16,fontweight='bold')
    ax2.set_ylim(FCH4_DICT['hist_range'])
    ax.set_title('Total Methane Emission for Great Britain',
                    fontsize=20,fontweight='bold')
    fig.savefig(plot_dir+'fch4_wetl/Histogram_'+str(start_year)+str(end_year)+'.png')
    plt.close()
    
PLOT_CO2=True 
if PLOT_CO2:
    # CO2 flux plots:
    os.system('mkdir -p '+plot_dir+'co2_flux/')
    CMAP=brewer2mpl.get_map('YlOrBr','Sequential','9')
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=1.0)
    CMAP_month=brewer2mpl.get_map('PuOr','Diverging','9',reverse=True)
    MPL_CMAP_month=CMAP_month.get_mpl_colormap(N=100,gamma=0.5)
   
    text_units='TgCO2'
    Ann_plot_units_mols='umolCO2 m$^{-2}$ s$^{-1}$'
    Ann_plot_units_mass='gC m$^{-2}$ day$^{-1}$'
    Mon_plot_units_mols='umolCO2 m$^{-2}$ s$^{-1}$'
    Mon_plot_units_mass='gC m$^{-2}$ day$^{-1}$'
   
    for var in ['npp_gb','gpp_gb']:
        print('Plotting '+var)
        VARDICT=JULES_DATA[var]
        varplotdir=plot_dir+'co2_flux'+'/'+var+'/'
        ANNUAL_CO2=False
        if ANNUAL_CO2:
          print('Annual uptake: ')
          for i in range(nYEARS):
            print(str(i+start_year))
            FILE_NAME = varplotdir+'Annual_'+var+'_mass_'+\
                                str(i+start_year)+'.png'
            PLOT_TITLE='Annual Mean '+VARDICT['long_name']+'\n'+str(i+start_year)
            PLOT_DATA=VARDICT['Annual-Mean'][i,:]*1e3*(86400.)
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                         DATA_RANGE=JULES_DATA[var]['data_range'],TICK_FORMAT='%0.1f',
                         LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                         CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                         CBAR_LABEL=Ann_plot_units_mass,
                         CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                         PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                         )
            
            FILE_NAME = varplotdir+'Annual_'+var+'_mols_'+\
                                str(i+start_year)+'.png'
            PLOT_TITLE='Annual Mean '+VARDICT['long_name']
            PLOT_DATA=VARDICT['Annual-Mean'][i,:]*kgC_to_molesC*1e6
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                         DATA_RANGE=JULES_DATA[var]['data_range'],TICK_FORMAT='%0.1f',
                         LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                         CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                         CBAR_LABEL=Ann_plot_units_mols,
                         CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                         PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                         )

        MONTHLY_CO2=False
        if MONTHLY_CO2:
          print('Monthly uptake: ')
          for i in range(nTSTEPS):
            datestring=JULES_time[i].strftime('%Y-%m')
            print(datestring)
            FILE_NAME = varplotdir+var+'_mass_'+\
                                datestring+'.png'
            PLOT_TITLE='Annual Mean '+VARDICT['long_name']+'\n '+datestring
            PLOT_DATA=VARDICT['data'][i,:]*(1e3*(86400.))
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                   DATA_RANGE=JULES_DATA[var]['data_range_MN'],TICK_FORMAT='%0.1f',
                   LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                   CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                   CBAR_LABEL=Mon_plot_units_mass,
                   CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                   PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                   )
            
            FILE_NAME = varplotdir+var+'_mols_'+\
                                datestring+'.png'
            PLOT_TITLE='Annual Mean '+VARDICT['long_name']+'\n '+datestring
            PLOT_DATA=VARDICT['data'][i,:]*kgC_to_molesC*1e6
            PTs.plot_map(PLOT_DATA,lons_2D,lats_2D,
                   DATA_RANGE=JULES_DATA[var]['data_range_MN'],TICK_FORMAT='%0.1f',
                   LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                   CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,
                   CBAR_LABEL=Mon_plot_units_mols,
                   CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                   PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                   )
    

        fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
        ax2=ax.twinx()
        ax.bar(JULES_time,VARDICT['Monthly-TS']*1e-9,365./12.,\
                color=CMAP.hex_colors[-5],lw=1,edgecolor=CMAP.hex_colors[-2])
        ax2.plot(JULES_time[5::12],VARDICT['Annual-TS']*1e-9,\
                 color=CMAP.hex_colors[-1],lw=3)
        ax.set_xlim([dt.date(start_year,1,1),dt.date(end_year+1,1,1)]) 
        ax.set_xlabel('Date',fontsize=13)
        ax.set_ylabel('TgCO$_2$ per Month',color=CMAP.hex_colors[-5],
                    fontsize=16,fontweight='bold')
        ax2.set_ylabel('TgCO$_2$ per Annum',color=CMAP.hex_colors[-1],
                    fontsize=16,fontweight='bold')
        ax2.set_ylim(VARDICT['hist_range'])
        ax.set_title('Total '+VARDICT['long_name']+' for Great Britain',
                    fontsize=20,fontweight='bold')
        fig.savefig(varplotdir+'Histogram_'+str(start_year)+str(end_year)+'.png')
        plt.close()
    

