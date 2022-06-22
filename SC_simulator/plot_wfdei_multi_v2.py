#!/usr/bin/env python
################################################################################
# 
# Program: plot_OPERATIONALoutput.py 
# 
# Python Script to plot the test OPERATIONAL output provided by the 
#     JASMIN simulation
#  Data located: /users/eow/edwcom/SC_simulator/OPERATIONAL_output/work_GLOBAL/
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
# PLOTS[0]:  GLOBAL Timeseries 
# PLOTS[1]:  TRANSCOM Timeseries 
# PLOTS[2]:  GLOBAL Map
# PLOTS[3]:  Zonal time-series
################################################################################
#

import pylab as plt
import numpy as np
import netCDF4 as nc
import sys, os
import netcdftime as nctime

import plot_tools as PT
import data_info_SC 
###################################################################################
# Define functions
###################################################################################
# Round to given significant figures
def round2SignifFigs(vals,n):
    mags = 10.0**np.floor(np.log10(np.abs(vals))) # order of mag's
    outvals = np.around(vals/mags,n-1)*mags       # round(val/omag)*omag
    try:
        outvals[np.where(np.isnan(vals))] = 0.0   # where order of mag = 0, set to zero
    except:
        if np.isnan(outvals):
            outvals=0.0
    #
    return outvals
#
###################################################################################


INTERACTIVE=sys.argv[1]
iDISPLAY=sys.argv[2]

#INTERACTIVE='N'
#iDISPLAY='N'
#PLOTS='YYYYYYYYYYNN'

fill_value=-999.
CS_pools = ['DPM', 'RPM', 'BIO', 'HUM']
CS_units = '$kg$C $m^{-2}$'
if '-CS_maxes' in sys.argv:
    argloc = sys.argv.index('-CS_maxes')
    temp=sys.argv.pop(argloc)
    temp_CSmaxes = sys.argv.pop(argloc)
    CS_maxes=temp_CSmaxes.replace('[','').replace(']','').split(',')
    CS_maxes=[float(CSmax) for CSmax in CS_maxes]
    print temp[1:], ' = ',CS_maxes,type(CS_maxes),len(CS_maxes)
    del temp, temp_CSmaxes
else:
    CS_maxes = [0.1,5,0.4,12]
    #CS_maxes = [10,10,10,10,]

nPOOLs = len(CS_pools)


BASE_DIR      = '/users/eow/edwcom/SC_simulator/'
OUTPUT_DIR    = BASE_DIR+'plots/'
GRID_FILE     = '/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
TRANSCOM_FILE = '/users/eow/edwcom/WFD_EI/TRANSCOM_Regions_WFDEI_landpoints.nc'

JULES_sources= data_info_SC.jules_WFDEI_sources()
if INTERACTIVE=='Y':
    print 'Available Sources: '
    for i,fn in zip(range(len(JULES_sources)),JULES_sources):
        print i,':',fn[0]
    temp_sources = raw_input('Enter sources to plot seperated by commas: ')
    SOURCES = temp_sources.split(',')
    del temp_sources
else:
    temp_sources=sys.argv[3]
    SOURCES=temp_sources.replace('[','').replace(']','').split(',')
    del temp_sources

if INTERACTIVE=='Y':
    PLOTS=''
    PLOTS+=raw_input('Produce GLOBAL Timeseries? (Y/N) ')
    PLOTS+=raw_input('Produce TRANSCOM Timeseries? (Y/N) ')
    PLOTS+=raw_input('Produce Global Soil Carbon Maps of final simulation year? (Y/N) ')
    PLOTS+=raw_input('Produce Zonal Time-series of first CS dataset? (Y/N) ') 
    PLOTS+=raw_input('Produce Regional Soil Carbon Maps of final simulation year? (Y/N) ')
else:
    PLOTS=sys.argv[4]

if PLOTS[4]=='Y':
    # Set sub-region as UK
    SUBREGION_NAME='UK'
    SUBREGION_LIMITS=np.array([49,-11,61,2.5])
    if INTERACTIVE=='Y':
        change_temp = raw_input('Default region is '+SUBREGION_NAME+\
                                ' with limits='+str(SUBREGION_LIMITS)+\
                                '\nDo you want to define a different region? (Y/N) ' )
        if change_temp=='Y':
            SUBREGION_NAME=raw_input('Give a name for the subregion: ')
            SUBREGION_LIMITS=input('Enter region limits as min_lat,min_lon,max_lat,max_lon: ')
    else:
        if '-region_name' in sys.argv:
            argloc = sys.argv.index('-region_name')
            temp=sys.argv.pop(argloc)
            SUBREGION_NAME= sys.argv.pop(argloc)
            del temp
        print 'Region = ',SUBREGION_NAME
        if '-region_limits' in sys.argv:
            argloc = sys.argv.index('-region_limits')
            temp=sys.argv.pop(argloc)
            temp_limits= sys.argv.pop(argloc)
            SUBREGION_LIMITS=np.array([float(limit) for limit in temp_limits])
            del temp
            del temp_limits
        print 'Region Limits = ',SUBREGION_LIMITS


plot_tag=''
if INTERACTIVE=='Y':
    plot_tag+=raw_input('Enter a tag for the plot filenames: ')
elif '-plot_tag' in sys.argv:
    argloc = sys.argv.index('-plot_tag')
    temp=sys.argv.pop(argloc)
    plot_tag+= sys.argv.pop(argloc)
    print temp[1:],'=',plot_tag
    del temp

OUTPUT_DIR+=plot_tag

# Read grid file
print 'Reading gridfile: '+GRID_FILE
grinf=nc.Dataset(GRID_FILE,'r')
grindex=grinf.variables['land_index'][:]-1
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
grimask=np.ones_like(grindex)
lons_2d,lats_2d = np.meshgrid(lons,lats)

# Read TRANSCOM file 
print 'Reading TRANSCOM file: '+TRANSCOM_FILE
trinf = nc.Dataset(TRANSCOM_FILE,'r')
TRANSCOM_REGIONS = trinf.variables['transcom_regions'][:]
trinf.close()
#TRAindex = range(360,720)
#for i in range(360):
#    TRAindex.append(i)
#TRANSCOM_REGIONS = TRANSCOM_REGIONS[:,TRAindex]
TRANSCOM_Names = [ 'North American Boreal',    \
                   'North American Temperate', \
                   'South American Tropical',  \
                   'South American Temperate', \
                   'Northern Africa',          \
                   'Southern Africa',          \
                   'Eurasia Boreal',           \
                   'Eurasia Temperate',        \
                   'Tropical Asia',            \
                   'Australia',                \
                   'Europe',                   ]


CS_DATA   = []
TIME_DATA = []
FRAC_DATA = []
print 'Reading in CS data'
for iSOURCE in SOURCES:
    iS = int(iSOURCE)
    print JULES_sources[iS][0]
    inf = nc.Dataset(JULES_sources[iS][1],'r')
    
    if (iSOURCE==SOURCES[0]):
        CSlats=inf.variables['latitude'][:].squeeze()
        CSlons=inf.variables['longitude'][:].squeeze()

    cs    = inf.variables['cs'][:].squeeze()
    frac  = inf.variables['frac'][:].squeeze()
    time  = nctime.num2date(inf.variables['time'][:],\
            units=inf.variables['time'].units,       \
            calendar='standard'                      )

    temp_Imask= frac[:,8,:]>0.1
    ICEmask = np.array( [ temp_Imask for i in range(nPOOLs) ] ).transpose(1,0,2)
    del temp_Imask

    cs[ ICEmask==True ] = fill_value
    cs=np.ma.masked_equal(cs,fill_value)
    
    CS_DATA.append( cs ) 
    TIME_DATA.append( time )
    FRAC_DATA.append( frac ) 
                       


# Global time series plots
if (PLOTS[0]=='Y'):
    
    # Create figure and plot axis if plotting global timeseries
    FIG = plt.figure(figsize=[18,12])
    nplts_wdth   = 2.
    nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
    for iPOOL in range(nPOOLs):
        #print CS_pools[iPOOL]
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
        AX.set_title(CS_pools[iPOOL])
        for iSOURCE in SOURCES:
            iS  = int(iSOURCE)
            iCS = SOURCES.index(iSOURCE)

            CS_timeseries \
                = np.array( [ np.mean(CS_DATA[iCS][i,iPOOL,:]) \
                               for i in range(CS_DATA[iCS].shape[0]) ] )
            
            AX.plot(TIME_DATA[iCS],CS_timeseries,label=JULES_sources[iS][0],
                    color=JULES_sources[iS][2],lw=2)
        
        
        AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])
        #AX.set_yscale('log')

    AX.legend( bbox_to_anchor=(-0.1,-0.11),loc=10,borderaxespad=0.,ncol=min(len(SOURCES),6) )
    FIG.text(0.03,0.5,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
    FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])
    FIG.suptitle('Global Time-series of the 4 Soil Carbon Pools',fontsize=24)
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_Global_Timeseries.png')
        plt.close()


# Timeseries by Transcom region 
if (PLOTS[1]=='Y'):
    # Loop round each CS pool
    #iPOOL=3
    for iPOOL in range(nPOOLs):
        FIG=plt.figure(figsize=[20,14])
        nplts_wdth = 4
        nplts_hght = 3
        
        #PLOT the Global time-series in top right corner
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,1)
        AX.set_title('GLOBAL')
        for iSOURCE in SOURCES:
            iS  = int(iSOURCE)
            iCS = SOURCES.index(iSOURCE)
            CS_timeseries \
                = np.array( [ np.mean(CS_DATA[iCS][i,iPOOL,:]) \
                              for i in range(CS_DATA[iCS].shape[0]) ] )
            
            AX.plot(TIME_DATA[iCS],CS_timeseries,
                    color=JULES_sources[iS][2],lw=2)

        AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])

        for iTRANS in range(len(TRANSCOM_Names)):        
            # Create axis for transom region
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iTRANS+2)
            AX.set_title(TRANSCOM_Names[iTRANS])
            
            # Get index of TRANSCOM regions
            TRindex = np.where(TRANSCOM_REGIONS==iTRANS+1)
            
            for iSOURCE in SOURCES:
                iS  = int(iSOURCE)
                iCS = SOURCES.index(iSOURCE)
                CS_timeseries \
                    = np.array( [ np.mean(CS_DATA[iCS][i,iPOOL,TRindex[0]]) \
                                  for i in range(CS_DATA[iCS].shape[0]) ] )
                
                AX.plot(TIME_DATA[iCS],CS_timeseries,label=JULES_sources[iS][0], \
                        color=JULES_sources[iS][2],lw=2)
        
            AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])
    
        AX.legend( bbox_to_anchor=(-1.1,-0.15),loc=10,borderaxespad=0.,ncol=min(len(SOURCES),6) )
        FIG.text(0.03,0.5,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
        FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])
        FIG.suptitle('Time-series of the '+CS_pools[iPOOL]+' soil carbon pool by TRANSCOM region',\
                     fontsize=26)
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+CS_pools[iPOOL]+'_SoilCarbon_TRANSCOM_Timeseries.png')
            plt.close()


# Plot global maps of mean soil carbon for final year
if (PLOTS[2]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
        FIG = plt.figure(figsize=[18,12])
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 2.
        nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
        for iPOOL in range(nPOOLs):
            #print CS_pools[iPOOL]
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
            AX.set_title(CS_pools[iPOOL])

            CS_mapdata =  np.mean(CS_DATA[iCS][-12:,iPOOL,:],axis=0) 
            # convert to 2d
            CS_mapdata = CS_mapdata[grindex]*grimask
            maxval = CS_maxes[iPOOL] 
        
            PT.plot_map(CS_mapdata,lons_2d,lats_2d,            \
                        DATA_RANGE=[0,maxval],                 \
                        PLOT_TITLE=CS_pools[iPOOL],            \
                        COLOURS=['white','orange','brown'],    \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,   \
                        CBAR_LABEL='Soil Carbon '+CS_units,    \
                        RESOLUTION='c',FONTSIZES=[12,12,15,18],\
                        AXIS=AX,     )

        
        FIG.tight_layout(rect=[0.0,0.05,1.0,0.92],h_pad=6)
        FIG.suptitle( 'Global Soil Carbon Map, '+JULES_sources[iS][0].replace('_','-'),fontsize=32)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'Global_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


#Plot timeseries of zonal bands (Global, NH,Tropical and SH) of WFDEI-BigSpin (or CS_DATA[0] for example plots)
if (PLOTS[3]=='Y'):
    # Create figure and plot axis if plotting global timeseries
    FIG = plt.figure(figsize=[18,12])
    nplts_wdth   = 2.
    nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
    
    index_names = ['Global','Northern Lats','Tropical','Southern Lats']
    indexes     = [ CSlats>-9999., CSlats>30 , (CSlats<15)&(CSlats>-15), CSlats<-30 ]
    colours     = [ 'red' , 'blue', 'orange', 'darkgreen' ]

    for iPOOL in range(nPOOLs):
        #print CS_pools[iPOOL]
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
        AX.set_title(CS_pools[iPOOL])

        for iZON in range(len(indexes)):
            index=indexes[iZON]
            # print index.shape
            CS_timeseries \
                = np.array( [ np.mean(CS_DATA[0][i,iPOOL,index]) \
                              for i in range(CS_DATA[0].shape[0]) ] )
            
            AX.plot(TIME_DATA[0],CS_timeseries,label=index_names[iZON],
                    color=colours[iZON],lw=2)
        
        
        AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])

    AX.legend( bbox_to_anchor=(-0.1,-0.11),loc=10,borderaxespad=0.,ncol=4 )
    FIG.text(0.03,0.5,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
    FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])
    FIG.suptitle('Zonal Time-series of the 4 Soil Carbon Pools',fontsize=24)
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_Zonal_Timeseries.png')
        plt.close()



# Plot sub-region maps of mean soil carbon for final year
if (PLOTS[4]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
        # caluculate figure size based on lat_range/lon_range ratio
        lon_range=(SUBREGION_LIMITS[1]-SUBREGION_LIMITS[3])
        lat_range=(SUBREGION_LIMITS[0]-SUBREGION_LIMITS[2])
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(9*xy_ratio),12])
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 2.
        nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
        for iPOOL in range(nPOOLs):
            #print CS_pools[iPOOL]
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
            AX.set_title(CS_pools[iPOOL])

            CS_mapdata =  np.mean(CS_DATA[iCS][-12:,iPOOL,:],axis=0) 
            # convert to 2d
            CS_mapdata = CS_mapdata[grindex]*grimask
            maxval = CS_maxes[iPOOL] 
        
            PT.plot_map(CS_mapdata,lons_2d,lats_2d,                            \
                        DATA_RANGE=[0,maxval],                                 \
                        PLOT_TITLE=CS_pools[iPOOL],                            \
                        COLOURS=['white','orange','brown'],                    \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,                   \
                        CBAR_LABEL='Soil Carbon '+CS_units,TICK_FORMAT='%0.2f',\
                        RESOLUTION='i',FONTSIZES=[12,12,15,18],                \
                        LON_RANGE=[SUBREGION_LIMITS[1],SUBREGION_LIMITS[3]],   \
                        LAT_RANGE=[SUBREGION_LIMITS[0],SUBREGION_LIMITS[2]],   \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),       \
                        AXIS=AX,     )

        
        FIG.tight_layout(rect=[0.0,0.05,1.0,0.92],h_pad=6)
        FIG.suptitle( SUBREGION_NAME+' Soil Carbon Map, '+JULES_sources[iS][0].replace('_','-'),fontsize=30*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+SUBREGION_NAME+'_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


