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
PLOTS=sys.argv[2]
iDISPLAY=sys.argv[3]

#INTERACTIVE='N'
#iDISPLAY='N'
#PLOTS='YYYYYYYYYYNN'

fill_value=-999.
CS_pools = ['DPM', 'RPM', 'BIO', 'HUM']
CS_units = '$kg$C $m^{-2}$'
CS_maxes = [0.1,10,1.0,25]

if INTERACTIVE=='Y':
    PLOTS=''
    PLOTS+=raw_input('Produce GLOBAL Timeseries? (Y/N)')
    PLOTS+=raw_input('Produce TRANSCOM Timeseries? (Y/N)')
    PLOTS+=raw_input('Produce ')


DATA_DIR      = '/users/eow/edwcom/SC_simulator/OPERATIONAL_output/work_GLOBAL/'
OUTPUT_DIR    = DATA_DIR+'plots/'
GRID_FILE     = '/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
TRANSCOM_FILE = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5.nc'

spin_files = [ 'J4.3_WFDEI_GLOBAL_spin1.monthly_mean.nc', \
               'J4.3_WFDEI_GLOBAL_spin2.monthly_mean.nc', \
               'J4.3_WFDEI_GLOBAL_spin3.monthly_mean.nc'  ]

BigSpin_file='/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_nocomp/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin_nocomp.monthly_mean.nc'

Koven_file='/users/eow/edwcom/SC_simulator/JULES_output/KovenMethod/JULES_v4.3_TRIFFID_RsQ10_WFDEI_GLOBAL_KovenMethod_spin2.monthly_mean.nc'

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
TRANSCOM_REGIONS = nc.Dataset(TRANSCOM_FILE,'r').\
                   variables['transcom_regions'][:]
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


# Read Big Spin data for comparison
print 'Reading Big Spin CS data'
inf = nc.Dataset(BigSpin_file,'r')
BigSpin_CS_DATA    = inf.variables['cs'][:].squeeze()
BigSpin_TIME_DATA = nctime.num2date(inf.variables['time'][:],          \
                                    units=inf.variables['time'].units, \
                                    calendar='standard'                )
BS_frac = inf.variables['frac'][0,:].squeeze()
BS_frac_2D = BS_frac[:,grindex]*grimask
BS_frac_2D.data[BS_frac_2D.mask==True]=fill_value
BS_frac_2D.fill_value=fill_value
inf.close()
ICEmask = BS_frac_2D[8,:]>=0.1

# Convert data to 2D and apply ocean mask
BigSpin_CS_DATA_2D = BigSpin_CS_DATA[:,:,grindex]*grimask
BigSpin_CS_DATA_2D.data[BigSpin_CS_DATA_2D.mask==True]=fill_value

# Apply ICE mask
BigSpin_CS_DATA_2D[:,:,ICEmask]=fill_value
BigSpin_CS_DATA_2D = np.ma.masked_equal(BigSpin_CS_DATA_2D.data,fill_value)


# Read Big Spin data for comparison
print 'Reading Koven CS data'
inf = nc.Dataset(Koven_file,'r')
Koven_CS_DATA    = inf.variables['cs'][:].squeeze()
Koven_TIME_DATA = nctime.num2date(inf.variables['time'][:],          \
                                    units=inf.variables['time'].units, \
                                    calendar='standard'                )
inf.close()

# Convert data to 2D and apply ocean mask
Koven_CS_DATA_2D = Koven_CS_DATA[:,:,grindex]*grimask
Koven_CS_DATA_2D.data[Koven_CS_DATA_2D.mask==True]=fill_value

# Apply ICE mask
Koven_CS_DATA_2D[:,:,ICEmask]=fill_value
Koven_CS_DATA_2D = np.ma.masked_equal(Koven_CS_DATA_2D.data,fill_value)


# Read SC simulator cs data for each of the 3 spin ups
print 'Reading SC simulator data'
CS_DATA = []
CS_time_DATA = []
for file in spin_files:
    inf = nc.Dataset(DATA_DIR+file,'r')
    cs_data = inf.variables['cs'][:].squeeze()
    cs_data_2D = cs_data[:,:,grindex]*grimask
    cs_data_2D.data[cs_data_2D.mask==True]=fill_value
    cs_data_2D[:,:,ICEmask]=fill_value
    cs_data_2D=np.ma.masked_equal(cs_data_2D,fill_value)
    CS_DATA.append(cs_data_2D)
    
    cs_time_data = nctime.num2date(inf.variables['time'][:],          \
                                   units=inf.variables['time'].units, \
                                   calendar='standard'                )
    
    CS_time_DATA.append(cs_time_data)

    inf.close()


# Global time series plots
if (PLOTS[0]=='Y'):
    
    # Create figure and plot axis if plotting global timeseries
    FIG = plt.figure(figsize=[18,12])
    nplts_wdth   = 2.
    nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
    for iPOOL in range(len(CS_pools)):
        #print CS_pools[iPOOL]
        BS_CS_Global_timeseries \
            = np.array( [ np.mean(BigSpin_CS_DATA_2D[i,iPOOL,:]) \
                          for i in range(BigSpin_CS_DATA_2D.shape[0]) ] )
        
        CS_Global_timeseries = []
        for data in CS_DATA:
            CS_Global_timeseries.append( \
                    np.array( [ np.mean(data[i,iPOOL,:]) \
                                for i in range(data.shape[0]) ] ) \
                                         )
        
        # Plot data if plotting global timeseries
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
        AX.set_title(CS_pools[iPOOL])
        # Plot Big Spin data
        AX.plot(BigSpin_TIME_DATA,BS_CS_Global_timeseries,
                label='BigSpin Run')
        # Plot the three spin ups
        for i in range(3):
            AX.plot(CS_time_DATA[i],CS_Global_timeseries[i],\
                    label='SC sim - spinup '+str(i+1))
        
        AX.set_ybound(lower=0,upper=np.max(BS_CS_Global_timeseries)*2)
    
    AX.legend()
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_Global_Timeseries.png')


# Timeseries by Transcom region 
# Only plot the Hummus
if (PLOTS[1]=='Y'):
    FIG=plt.figure(figsize=[24,18])
    nplts_wdth = 3
    nplts_hght = 4
    iPOOL=3
        
    # Calculate Global Time series
    BS_CS_Global_timeseries \
        = np.array( [ np.mean(BigSpin_CS_DATA_2D[i,iPOOL,:]) \
                      for i in range(BigSpin_CS_DATA_2D.shape[0]) ] )
    CS_Global_timeseries = []
    for data in CS_DATA:
        CS_Global_timeseries.append( \
                np.array( [ np.mean(data[i,iPOOL,:]) \
                for i in range(data.shape[0]) ] ) \
                                     )

    # Plot the Global data in the top left panel
    AX  = FIG.add_subplot(nplts_hght,nplts_wdth,1)
    AX.set_title('Global')
    # Plot Big Spin data
    BS_line=AX.plot(BigSpin_TIME_DATA,BS_CS_Global_timeseries,\
                    label='Big Spin Run')
    CS1_line = AX.plot(CS_time_DATA[0],CS_Global_timeseries[0],\
                       label='SC sim - spinup 1')
    CS2_line = AX.plot(CS_time_DATA[1],CS_Global_timeseries[1],\
                       label='SC sim - spinup 2')
    CS3_line = AX.plot(CS_time_DATA[2],CS_Global_timeseries[2],\
                       label='SC sim - spinup 3')
    AX.legend()
    AX.set_ybound(lower=0,upper=np.max(BS_CS_Global_timeseries)*2)
    
    for iTRANS in range(len(TRANSCOM_Names)):        
        # Get index of TRANSCOM regions
        TRindex = np.where(TRANSCOM_REGIONS==iTRANS+1)
        BS_CS_TRANS_timeseries \
          = np.array( [ np.mean(BigSpin_CS_DATA_2D[i,iPOOL,TRindex[0],TRindex[1]]) \
                        for i in range(BigSpin_CS_DATA_2D.shape[0]) ] )
        CS_TRANS_timeseries = []
        for data in CS_DATA:
            CS_TRANS_timeseries.append( \
                    np.array( [ np.mean(data[i,iPOOL,:]) \
                                for i in range(data.shape[0]) ] ) \
                                         )
        
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iTRANS+2)
        AX.set_title(TRANSCOM_Names[iTRANS])
        # Plot Big Spin data
        AX.plot(BigSpin_TIME_DATA,BS_CS_TRANS_timeseries,
                label='BigSpin Run')
        # Plot the three spin ups
        for i in range(3):
            AX.plot(CS_time_DATA[i],CS_TRANS_timeseries[i],\
                    label='SC sim - spinup '+str(i+1))
    
        AX.set_ybound(lower=0,upper=np.max(BS_CS_TRANS_timeseries)*1.5)
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_TRANSCOM_Timeseries.png')


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




# Global time series with Koven data
if (PLOTS[3]=='Y'):
    
    # Create figure and plot axis if plotting global timeseries
    FIG = plt.figure(figsize=[18,12])
    nplts_wdth   = 2.
    nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
    for iPOOL in range(len(CS_pools)):
        #print CS_pools[iPOOL]
        BS_CS_Global_timeseries \
            = np.array( [ np.mean(BigSpin_CS_DATA_2D[i,iPOOL,:]) \
                          for i in range(BigSpin_CS_DATA_2D.shape[0]) ] )
       
        # Just take final simulator run
        CS_Global_timeseries \
            = np.array( [ np.mean(CS_DATA[2][i,iPOOL,:]) \
                            for i in range(CS_DATA[2].shape[0]) ] ) 

        # Koven Data
        Koven_Global_timeseries \
            = np.array( [ np.mean(Koven_CS_DATA_2D[i,iPOOL,:]) \
                          for i in range(Koven_CS_DATA_2D.shape[0]) ] )
       

        # Plot data if plotting global timeseries
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
        AX.set_title(CS_pools[iPOOL])
        # Plot Big Spin data
        AX.plot(BigSpin_TIME_DATA,BS_CS_Global_timeseries,
                label='BigSpin Run')
        # Plot the third spin ups
        AX.plot(CS_time_DATA[2],CS_Global_timeseries,\
                label='SC simulator')
        # Plot the Koven data
        AX.plot(Koven_TIME_DATA,BS_CS_Global_timeseries,
                label='Koven Run')
        
        #AX.set_ybound(lower=0,upper=np.max(BS_CS_Global_timeseries)*2)
    
    AX.legend()
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_Global_Timeseries_withKoven.png')

