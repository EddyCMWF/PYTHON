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
JULES_vars=['npp_gb','gpp_gb','resp_s_gb','cs','resp_p_gb',\
            'fch4_wetl','smcl', 'tstar_gb']

JULES_CHESS_gridfile='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

CHESS_landcover_file='/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'

kgC_to_molesC = ( 1e3/12.01 )

year=2014
month_names=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
days_in_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
nMONTHs=12
hist_month_objects=[ dt.date(year,imonth+1,1) for imonth in range(nMONTHs) ]
line_month_objects=[ dt.date(year,imonth+1,15) for imonth in range(nMONTHs) ]
print(days_in_month.shape)

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
Jinf=nc.Dataset(JULES_dir+JULES_filetag.replace('YYYY',str(year)),'r')
JULES_DATA={}
for var in JULES_vars:
    invar = Jinf.variables[var]
    indata = invar[:].squeeze()[...,grindex]
    indata = indata*np.ones_like(grindex)
    indata.data[indata.mask==True]=indata.fill_value
    JULES_DATA[var]={'data': indata, \
                     'units': invar.units, \
                     'long_name': invar.long_name }

Jinf.close()   

# convert units accordingly
JULES_DATA['cs']['units']='kg m$^{-2}$'

JULES_DATA['fch4_wetl']['data']*=1e-9 #
JULES_DATA['fch4_wetl']['units']='kg m$^{-2}$ s$^{-1}$'

JULES_DATA['resp_s_gb']['data']*=0.162
# scales to P.Levy natural co2 flux (0.45 kg per year)
JULES_DATA['resp_s_gb']['units']='kg m$^{-2}$ s^{-1}$'
JULES_DATA['resp_s_gb']['long_name']='Soil Respiration'
JULES_DATA['resp_s_gb']['data_range']=[0,3]

#JULES_DATA['npp_gb']['data']*=86400.
#JULES_DATA['npp_gb']['units']='kg m$^{-2}$ s$^{-1}$'
JULES_DATA['npp_gb']['long_name']='Net Primary Productivity'
JULES_DATA['npp_gb']['data_range']=[2,3]

#JULES_DATA['gpp_gb']['data']*=86400.
#JULES_DATA['gpp_gb']['units']='kg m$^{-2}$ s$^{-1}$'
JULES_DATA['gpp_gb']['long_name']='Gross Primary Productivity'
JULES_DATA['gpp_gb']['data_range']=[0,4]

#JULES_DATA['resp_p_gb']['data']*=86400.
#JULES_DATA['resp_p_gb']['units']='kg m$^{-2}$ day$^{-1}$'
JULES_DATA['resp_p_gb']['long_name']='Plant Respiration'
JULES_DATA['resp_p_gb']['data_range']=[1,2]

JULES_DATA['smcl']['data']=np.sum(JULES_DATA['smcl']['data'][:,:3,:],axis=1)*1e-3
JULES_DATA['smcl']['units']='m^3 m^-3'
JULES_DATA['smcl']['long_name']='Soil moisture content of top 1m'

JULES_DATA['nee_gb']={'data':JULES_DATA['npp_gb']['data']-JULES_DATA['resp_s_gb']['data'],
                      'units':'kg km$^{-2}$ day$^{-1}$', 
                      'long_name':'Net Ecosystem Exchange', 
                      'data_range':[0,4] }


PLOT_FCH4=False
if PLOT_FCH4:
    # Methane emission plots:
    os.system('mkdir -p '+plot_dir+'fch4_wetl/')
    CMAP=brewer2mpl.get_map('YlGnBu','Sequential','9',reverse=True)
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=0.8)
    
    #Unit strings for plots
    text_units='GgCH4'
    Ann_plot_units_mols='mmolCH4 m$^{-2}$ s$^{-1}$'
    Ann_plot_units_mass='kgCH4 km$^{-2}$ yr$^{-1}$'
    Mon_plot_units_mols='mmolCH4 m$^{-2}$ s$^{-1}$'
    Mon_plot_units_mass='kgCH4 km$^{-2}$ mo$^{-1}$'
    
    #Emission maps by month in millimolesC m^-2 s^-1
    JULES_ch4_emis_mols=JULES_DATA['fch4_wetl']['data']*kgC_to_molesC*1e3
    # Annual mean emission rate
    JULES_Ann_ch4_emis_mols=np.mean(JULES_ch4_emis_mols,axis=0)

    #Emission maps by month in kgC km^-2 day^-1
    JULES_ch4_emis_mass=JULES_DATA['fch4_wetl']['data']*1e6*(16.04/12.01)
    #Annual Emission map of Total emission kgC km^-2 yr^-1
    JULES_Ann_ch4_emis_mass=np.mean(JULES_ch4_emis_mass,axis=0)*365.

    # Total Emission scaler GgCH4 yr^-1
    TOTAL_Ann_ch4_emis=np.sum(JULES_Ann_ch4_emis_mass)*1e-6
    
    # monthly CH4 emission GgC day^-1
    DAILY_ch4_emis_byMONTH=np.array([np.sum(JULES_ch4_emis_mass[imonth,:])*1e-6 \
                                        for imonth in range(nMONTHs)])
    
    outf=open(plot_dir+'fch4_wetl/fch4_wetl_info_'+str(year)+'.txt','w')
    outf.write('National Wetland Methane Emission Inventory '+ text_units+'\n')
    outf.write('Annual =%8.2f\n'%(TOTAL_Ann_ch4_emis))
    for imonth in range(nMONTHs):
        outf.write('%7a=%8.2f\n'%(month_names[imonth],\
                np.sum(DAILY_ch4_emis_byMONTH[imonth]*days_in_month[imonth]) ) )
    outf.close()

    PLOT_TITLE='Annual Wetland Methane Emission'
    FILE_NAME=plot_dir+'fch4_wetl/Annual_Wetland_Emission_mass.png'
    PTs.plot_map(JULES_Ann_ch4_emis_mass,lons_2D,lats_2D,
            DATA_RANGE=[0,2000],TICK_FORMAT='%0.1f',
            LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
            CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Ann_plot_units_mass,
            CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
            PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
            )
    
    FILE_NAME=plot_dir+'fch4_wetl/Annual_Wetland_Emission_mols.png'
    PTs.plot_map(JULES_Ann_ch4_emis_mols,lons_2D,lats_2D,
                 DATA_RANGE=[0,1.0],TICK_FORMAT='%0.1f',
                 LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                 CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Ann_plot_units_mols,
                 CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                 PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                 )
    
    for imonth in range(nMONTHs):
        PLOT_TITLE=month_names[imonth]+' Wetland Methane Emission'
        FILE_NAME=plot_dir+'fch4_wetl/%02i_'%(imonth+1)+month_names[imonth]+\
                    '_Wetland_Emission_mass.png'
        PTs.plot_map(JULES_ch4_emis_mass[imonth,:]*days_in_month[imonth],lons_2D,lats_2D,
                DATA_RANGE=[0,200],TICK_FORMAT='%0.1f',
                LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Mon_plot_units_mass,
                CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                )
        
        FILE_NAME=plot_dir+'fch4_wetl/%02i_'%(imonth+1)+month_names[imonth]+\
                    '_Wetland_Emission_mols.png'
        PTs.plot_map(JULES_ch4_emis_mols[imonth,:],lons_2D,lats_2D,
                DATA_RANGE=[0,1.0],TICK_FORMAT='%0.1f',
                LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=Mon_plot_units_mols,
                CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                )

    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    ax.bar(hist_month_objects,(DAILY_ch4_emis_byMONTH*days_in_month),days_in_month-0.4,\
            color=CMAP.hex_colors[5],lw=1,edgecolor=CMAP.hex_colors[1]) 
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel('GgCH4')
    ax.set_title('Total Methane Emission for Great Britain')
    fig.savefig(plot_dir+'fch4_wetl/Histogram_True.png')
    plt.close()
    
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    ax.bar(hist_month_objects,(DAILY_ch4_emis_byMONTH*365./12.),days_in_month-0.4,\
            color=CMAP.hex_colors[5],lw=1,edgecolor=CMAP.hex_colors[1]) 
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)]) 
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel('GgCH4')
    ax.set_title('Total Methane Emission for Great Britain')
    fig.savefig(plot_dir+'fch4_wetl/Histogram_EqualMonth.png')
    plt.close()


PLOT_CO2=False
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
   
    for var in ['npp_gb','gpp_gb','resp_s_gb','resp_p_gb','nee_gb']:
        var_plotdir=plot_dir+'co2_flux/'+var+'/'
        os.system('mkdir -p '+var_plotdir)
        # maps by month in molsC m^-2 s^-1
        MONTH_FLUX_mols=JULES_DATA[var]['data']*kgC_to_molesC*1e6
        #Annual NPP map imolsC m^-2 s^-1
        ANNUAL_FLUX_mols=np.mean(MONTH_FLUX_mols,axis=0)
        
        # maps by month in gC m^-2 day^-1
        MONTH_FLUX_mass=JULES_DATA[var]['data']*86400.*1e3
        #Annual NPP map gC m^-2 day^-1
        ANNUAL_FLUX_mass=np.mean(MONTH_FLUX_mass,axis=0)

        # Total Emission scaler kgC yr^-1
        TOTAL_Ann_FLUX=np.sum(ANNUAL_FLUX_mass) * 1e6 * 1e-12 *365
        #                               * m^2 to km^2 * kg to Tg
        
        # monthly CH4 emission kgC day^-1
        DAILY_FLUX_byMONTH=np.array([np.sum(MONTH_FLUX_mass[imonth,:])*1e6*1e-12 \
                                    for imonth in range(nMONTHs)])
    
        outf=open(var_plotdir+var+'_info_'+str(year)+'.txt','w')
        outf.write('National '+JULES_DATA[var]['long_name']+' Inventory ('+ text_units+')\n')
        outf.write('Annual =%8.2f\n'%(TOTAL_Ann_FLUX))
        for imonth in range(nMONTHs):
            outf.write('%7a=%8.2f\n'%(month_names[imonth],\
                        np.sum(DAILY_FLUX_byMONTH[imonth]*days_in_month[imonth]) ) )
        outf.close()
        
        PLOT_TITLE='Annual Mean '+JULES_DATA[var]['long_name']
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
    
   
PLOT_SMCL=True
if PLOT_SMCL:
    # Methane emission plots:
    os.system('mkdir -p '+plot_dir+'smcl/')
    CMAP=brewer2mpl.get_map('GnBu','Sequential','9')
    MPL_CMAP=CMAP.get_mpl_colormap(N=100,gamma=1.0)
    
    #Unit strings for plots
    plot_units='m$^3$ m$^{-3}$'
    data_range=[0,0.5]

    JULES_SMCL=JULES_DATA['smcl']['data']
    JULES_Ann_SMCL=np.mean(JULES_SMCL,axis=0)
    JULES_SMCL_timeseries=np.mean(JULES_SMCL.reshape(nMONTHs,-1),axis=1)
    JULES_SMCL_timeseries_std=np.std(JULES_SMCL.reshape(nMONTHs,-1),axis=1)
    

    PLOT_TITLE='Mean Annual Soil Moisture'
    FILE_NAME=plot_dir+'smcl/Annual_MeanSoilMoisture.png'
    PTs.plot_map(JULES_Ann_SMCL,lons_2D,lats_2D,
                 DATA_RANGE=data_range,TICK_FORMAT='%0.1f',
                 LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                 CMAP=MPL_CMAP,NLEVELS=100,NTICKS=6,CBAR_LABEL=plot_units,
                 CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                 PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                 )
    
    for imonth in range(nMONTHs):
        PLOT_TITLE=month_names[imonth]+' Mean Soil Temperature'
        FILE_NAME=plot_dir+'smcl/%02i_'%(imonth+1)+month_names[imonth]+\
                    '_MeanSoilMoisture.png'
        PTs.plot_map(JULES_SMCL[imonth,:],lons_2D,lats_2D,
                      DATA_RANGE=data_range,TICK_FORMAT='%0.1f',
                      LATDEL=2.5,LONDEL=2.5,PROJECTION='chess',LON_RANGE=[-8,4],
                      CMAP=MPL_CMAP,NLEVELS=100,NTICKS=11,CBAR_LABEL=plot_units,
                      CBAR_ORIENTATION='vertical',iDISPLAY='N',iCLOSE='Y',
                      PLOT_TITLE=PLOT_TITLE,FILE_PLOT=FILE_NAME,
                      )

    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(12,5))
    polygon_x=line_month_objects+line_month_objects[::-1]
    polygon_y=np.append(JULES_SMCL_timeseries-JULES_SMCL_timeseries_std,\
                        JULES_SMCL_timeseries[::-1]+JULES_SMCL_timeseries_std)
    polygon_y2=np.append(JULES_SMCL_timeseries-(JULES_SMCL_timeseries_std*2),\
                         JULES_SMCL_timeseries[::-1]+(JULES_SMCL_timeseries_std*2))
    ax.fill(polygon_x,polygon_y2,c=CMAP.hex_colors[1])
    ax.fill(polygon_x,polygon_y,c=CMAP.hex_colors[3])
    ax.plot(line_month_objects,JULES_SMCL_timeseries,\
            color=CMAP.hex_colors[-1],lw=2)
    ax.set_xlim([dt.date(year,1,1),dt.date(year+1,1,1)])
    ax.set_ylim([0,1])
    ax.set_xlabel('Date',fontsize=13)
    ax.set_ylabel(plot_units)
    ax.set_title('Mean Soil Moisture for Great Britain')
    fig.savefig(plot_dir+'smcl/TimeSeries.png')
    plt.close()
    
   




