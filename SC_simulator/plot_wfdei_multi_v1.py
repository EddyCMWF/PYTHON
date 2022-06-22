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
CS_maxes = [0.1,5,0.3,10]
nPOOLs = len(CS_pools)



BASE_DIR      = '/users/eow/edwcom/SC_simulator/'
OUTPUT_DIR    = BASE_DIR+'plots/'
GRID_FILE     = '/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
TRANSCOM_FILE = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5.nc'

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
    PLOTS+=raw_input('Produce GLOBAL Timeseries? (Y/N)')
    PLOTS+=raw_input('Produce TRANSCOM Timeseries? (Y/N)')
    PLOTS+=raw_input('Produce Soil Carbon Maps of final simulation year? (Y/N)')
else:
    PLOTS=sys.argv[4]

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
TRAindex = range(360,720)
for i in range(360):
    TRAindex.append(i)
TRANSCOM_REGIONS = TRANSCOM_REGIONS[:,TRAindex]
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

for iSOURCE in SOURCES:
    iS = int(iSOURCE)
    inf = nc.Dataset(JULES_sources[iS][1],'r')
    
    cs    = inf.variables['cs'][:].squeeze()
    frac  = inf.variables['frac'][:].squeeze()
    time  = nctime.num2date(inf.variables['time'][:],\
            units=inf.variables['time'].units,       \
            calendar='standard'                      )

    frac_2d = frac[:,:,grindex]*grimask
    frac_2d.data[frac_2d.mask==True]=fill_value
    temp_Imask= frac_2d[:,8,:]>0.1
    ICEmask = np.array( [ temp_Imask for i in range(nPOOLs) ] ).transpose(1,0,2,3)
    del temp_Imask

    cs_2d = cs[:,:,grindex]*grimask
    cs_2d.mask[ ICEmask ] = True
    cs_2d.data[cs_2d.mask==True]=fill_value
    
    CS_DATA.append( cs_2d ) 
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

    AX.legend( bbox_to_anchor=(-0.1,-0.1),loc=10,borderaxespad=0.,ncol=min(len(SOURCES),6) )
    FIG.tight_layout(rect=[0,0.05,1.0,0.96])
    FIG.suptitle('Global Time-series of the 4 Soil Carbon Pools',fontsize=20)
    
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
                    = np.array( [ np.mean(CS_DATA[iCS][i,iPOOL,TRindex[0],TRindex[1]]) \
                                  for i in range(CS_DATA[iCS].shape[0]) ] )
                
                AX.plot(TIME_DATA[iCS],CS_timeseries,label=JULES_sources[iS][0], \
                        color=JULES_sources[iS][2],lw=2)
        
            AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])
    
        AX.legend( bbox_to_anchor=(-1.1,-0.15),loc=10,borderaxespad=0.,ncol=min(len(SOURCES),6) )
        FIG.tight_layout(rect=[0,0.05,1.0,0.95])
        FIG.suptitle('Time-series of the '+CS_pools[iPOOL]+' soil carbon pool by TRANSCOMM region',\
                     fontsize=20)
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+CS_pools[iPOOL]+'_SoilCarbon_TRANSCOM_Timeseries.png')
            plt.close()


# Plot global maps of mean soil carbon for final year
if (PLOTS[2]=='Y'):
    
    for iSPINUP in range(3):
        spindata=CS_DATA[iSPINUP]

        # Create figure to plot on
        FIGURE=plt.figure(figsize=(24,14))
        for iPOOL in range(4):
            AXIS = FIGURE.add_subplot(2,2,iPOOL+1)
            plotdata=np.mean(spindata[-12:,iPOOL,:],axis=0)
            # maximum colour bar value = mean * 2 to 1SF
            maxval = CS_maxes[iPOOL]  #round2SignifFigs(np.mean(plotdata)*2,1)
            
            PT.plot_map(plotdata,lons_2d,lats_2d,              \
                        DATA_RANGE=[0,maxval],                 \
                        PLOT_TITLE=CS_pools[iPOOL],            \
                        COLOURS=['white','orange','brown'],    \
                        INTERPOLATE_COLOURS=True,NLEVELS=15,   \
                        CBAR_LABEL='Soil Carbon '+CS_units,    \
                        RESOLUTION='c',FONTSIZES=[12,12,15,18],\
                        AXIS=AXIS,\
                        )
         
        FIGURE.suptitle( 'Global Soil Carbon Map, '+ \
                         'Spinup-'+str(iSPINUP+1),fontsize=32)
        
        FILE_PLOT =  'Global_Soil_Carbon_Map_Spinup-'+str(iSPINUP+1)+'.png'   
        FIGURE.savefig(OUTPUT_DIR+FILE_PLOT, bbox_inches='tight',pad_inches=0.3)
        plt.close()

    # Big Spin run
    # Create figure to plot on
    FIGURE=plt.figure(figsize=(24,14))
    for iPOOL in range(4):
        AXIS = FIGURE.add_subplot(2,2,iPOOL+1)
        # get mean of final year of data for each pool
        plotdata=np.mean(BigSpin_CS_DATA_2D[-12:,iPOOL,:],axis=0)
        # maximum colour bar value = mean * 2 to 1SF
        maxval = CS_maxes[iPOOL] 
        
        PT.plot_map(plotdata,lons_2d,lats_2d,              \
                    DATA_RANGE=[0,maxval],                 \
                    PLOT_TITLE=CS_pools[iPOOL],            \
                    COLOURS=['white','orange','brown'],    \
                    INTERPOLATE_COLOURS=True,NLEVELS=15,   \
                    CBAR_LABEL='Soil Carbon '+CS_units,    \
                    RESOLUTION='c',FONTSIZES=[12,12,15,18],\
                    AXIS=AXIS,        )

    FIGURE.suptitle( 'Global Soil Carbon Map, '+ \
                     'Big Spin',fontsize=32)
        
    FILE_PLOT =  'Global_Soil_Carbon_Map_BigSpin.png'   
    FIGURE.savefig(OUTPUT_DIR+FILE_PLOT, bbox_inches='tight',pad_inches=0.3)
    plt.close()



