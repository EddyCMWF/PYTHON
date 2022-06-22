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
# PLOTS[0] - Soil Carbon Time-series for full CHESS grid
# PLOTS[1] - Time-series by landcover
# PLOTS[2] - Maps of Mean Soil Carbon Pool Maps for final simulation year
# PLOTS[3] - Maps of Mean Total Soil Carbon for final simulation year 
# PLOTS[4] - Maps of Mean Total Soil Carbon difference (source[0]-source[1])
# PLOTS[5] - Regional Maps of Mean Soil Carbon Pools for final simulation year
# PLOTS[6] - Regional Maps of Mean Total Soil Carbon for final simulation year
# PLOTS[7] - Regional Maps of Mean Total Soil Carbon difference (source[0]-source[1])
# PLOTS[8] - 
# PLOTS[9] -
# PLOTS[10] - Maps of Dominent Land Cover
# PLOTS[11] - Maps of Mean Litterfall for full simulation
# PLOTS[12] - Maps of Mean Litterfall by PFT for full simulation
# PLOTS[13] - Maps of Mean Soil Respiration for full simulation
# PLOTS[14] - Maps of Mean Soil Respiration by SCpool for full simulation
# PLOTS[15] - Maps of Mean Vegetation Carbon for full simulation
################################################################################
#

import pylab as plt
import numpy as np
import netCDF4 as nc
import sys, os

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
kappa_s_RothC= np.array([ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10 ])
CS_tickformats = ['%0.2f','%0.1f','%0.2f','%0.1f','%0.1f']
CS_units = '$kg$C $m^{-2}$'
CS_DIFF_RANGE=[-2,2]
DIFF_COL=['darkbrown','lawngreen']
nPOOLs = len(CS_pools)

if '-CS_maxes' in sys.argv:
    argloc = sys.argv.index('-CS_maxes')
    temp=sys.argv.pop(argloc)
    temp_CSmaxes = sys.argv.pop(argloc)
    CS_maxes=temp_CSmaxes.replace('[','').replace(']','').split(',')
    CS_maxes=[float(CSmax) for CSmax in CS_maxes]
    print temp[1:], ' = ',CS_maxes,type(CS_maxes),len(CS_maxes)
    del temp, temp_CSmaxes
else:
    CS_maxes = [0.16,8,0.4,16,20]
    #CS_maxes = [10,10,10,10,]

RESP_S_units = '$kg$C $m^{-2}$ $yr^{-1}$'
RESP_S_tickformats = ['%0.2f','%0.1f','%0.2f','%0.1f','%0.1f']
RESP_S_conv_factor = 3600.*24*360. # Seconds in a 360 day year
if '-RESP_S_maxes' in sys.argv:
    argloc = sys.argv.index('-RESP_S_maxes')
    temp=sys.argv.pop(argloc)
    temp_RESPSmaxes = sys.argv.pop(argloc)
    RESP_S_maxes=temp_RESPSmaxes.replace('[','').replace(']','').split(',')
    RESP_S_maxes=[float(RESPSmax) for RESPSmax in RESP_S_maxes]
    print temp[1:], ' = ',RESP_S_maxes,type(RESP_S_maxes),len(RESP_S_maxes)
    del temp, temp_RESPSmaxes
else:
    RESP_S_maxes = [0.4,0.8,0.1,0.1,1.4]

RESP_FACT_units = '$yr^{-1}$'
RESP_FACT_tickformats =  ['%0.1f','%0.2f','%0.2f','%0.3f','%0.1f']
RESP_FACT_NLEVELS     = [ 11, 7, 7, 6 , 11 ]
if '-RESP_FRAC_maxes' in sys.argv:
    argloc = sys.argv.index('-RESP_FRAC_maxes')
    temp=sys.argv.pop(argloc)
    temp_RESPFmaxes = sys.argv.pop(argloc)
    RESP_FRAC_maxes=temp_RESPFmaxes.replace('[','').replace(']','').split(',')
    RESP_FRAC_maxes=[float(RESPFmax) for RESPFmax in RESP_FRAC_maxes]
    print temp[1:], ' = ',RESP_FRAC_maxes,type(RESP_FRAC_maxes),len(RESP_FRAC_maxes)
    del temp, temp_RESPFmaxes
else:
    RESP_FRAC_maxes = [10,0.3,0.6,0.02,5]

PFT_names=['Broadleaf','Needleleaf','C3 Grass','Shrub','Crop']
nPFTs=5
LIT_C_units = '$kg$C $m^{-2}$ $yr^{-1}$'
LIT_C_tickformat = '%0.1f'
if '-LIT_C_max' in sys.argv:
    argloc = sys.argv.index('-LIT_C_max')
    temp=sys.argv.pop(argloc)
    LIT_C_max = float(sys.argv.pop(argloc))
    del temp
else:
    LIT_C_max = 1.4

C_VEG_units = '$kg$C $m^{-2}$'
C_VEG_tickformat = '%0.2f'
if '-C_VEG_max' in sys.argv:
    argloc = sys.argv.index('-C_VEG_maxes')
    temp=sys.argv.pop(argloc)
    temp_CVEGmaxes = sys.argv.pop(argloc)
    C_VEG_maxes=temp_CVEGmaxes.replace('[','').replace(']','').split(',')
    C_VEG_maxes=[float(C_VEGmax) for C_VEGmax in C_VEG_maxes]
    print temp[1:], ' = ',C_VEG_maxes,type(C_VEG_maxes),len(C_VEG_maxes)
    del temp, temp_CVEGmaxes
else:
    C_VEG_maxes = [4.,8.,0.4,0.8,0.4,8.]


BASE_DIR      = '/users/eow/edwcom/SC_simulator/'
OUTPUT_DIR    = BASE_DIR+'plots/'
GRID_FILE     = '/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

JULES_sources= data_info_SC.jules_CHESS_sources()
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
    # PLOTS[0]
    PLOTS+=raw_input('Produce time-series for full CHESS grid? (Y/N) ')
    # PLOTS[1]
    PLOTS+=raw_input('Produce Timeseries by dominent landcover? (Y/N) ')
    # PLOTS[2]
    PLOTS+=raw_input('Produce Maps of Mean Soil Carbon Pool Maps for final simulation year? (Y/N) ')
    # PLOTS[3]
    PLOTS+=raw_input('Produce Maps of Mean Total Soil Carbon for final simulation year? (Y/N) ')
    # PLOTS[4]
    PLOTS+=raw_input('Produce Maps of Mean Total Soil Carbon difference (source[0]-source[1])? (Y/N) ')
    # PLOTS[5]
    PLOTS+=raw_input('Produce Regional Maps of Mean Soil Carbon Pools for final simulation year? (Y/N) ')
    # PLOTS[6]
    PLOTS+=raw_input('Produce Regional Maps of Mean Total Soil Carbon for final simulation year? (Y/N) ')
    # PLOTS[7]
    PLOTS+=raw_input('Produce Regional Maps of Mean Total Soil Carbon difference (source[0]-source[1])? (Y/N) ')
    # PLOTS[8]
    PLOTS+='N'
    # PLOTS[9]
    PLOTS+='N'
    # PLOTS[10]
    PLOTS+=raw_input('Produce Maps of Dominent Land Cover? (Y/N) ')
    # PLOTS[11]
    PLOTS+=raw_input('Produce Maps of Mean Litterfall for full simulation? (Y/N) ')
    # PLOTS[12]
    PLOTS+=raw_input('Produce Maps of Mean Litterfall by PFT for full simulation? (Y/N) ')
    # PLOTS[13]
    PLOTS+=raw_input('Produce Maps of Mean Soil Respiration for full simulation? (Y/N) ')
    # PLOTS[14]
    PLOTS+=raw_input('Produce Maps of Mean Soil Respiration by SCpool for full simulation? (Y/N) ')
    # PLOTS[15]
    PLOTS+=raw_input('Produce Maps of Mean Vegetation Carbon for full simulation? (Y/N) ')
    # PLOTS[16]
    PLOTS+=raw_input('Produce Maps of Mean Vegetation Carbon by PFT for full simulation? (Y/N) ')
    # PLOTS[17]
    PLOTS+=raw_input('Produce Maps of Mean Soil Resp Factor for full simulation? (Y/N) ')
    # PLOTS[18]
    PLOTS+=raw_input('Produce Maps of Mean Soil Resp Factor by SCpool for full simulation? (Y/N) ')
    # PLOTS[19]
    PLOTS+=raw_input('Produce Maps of (kappa_dpm*C_dpm+kappa_rpm*C_rpm)/(kappa_bio*C_bio+kappa_hum*C_hum)? (Y/N) ')
else:
    PLOTS=sys.argv[4]

if PLOTS[1]=='Y':
    LC_names     = [ 'Broadleaf Tree','Needleleaf Tree','C3 Grass','shrub','Crops','Baresoil']
    LC_indexes   = [ 0,1,2,3,4,7 ]
    LC_threshold = 0.6
    print 'Default land covers to plot are:'
    print LC_names
    print 'With a minimum fraction for dominent cover of ', LC_threshold
    if INTERACTIVE=='Y':
        temp_change=raw_input('Do you want to change these options? (Y/N) ')
        if temp_change=='Y':
            temp_LC_names = raw_input('Enter lnadcover names seperated by commas: ')
            LC_names = temp_LC_names.split(',')
            LC_indexes = input('Enter corresponding indexes seperated by commas (from 0 index)')
            LC_threshold = input('Enter minimum land fraction for dominent cover')

if (PLOTS[5]=='Y') or (PLOTS[6]=='Y') or (PLOTS[7]=='Y'):
    # Set sub-region as UK
    SUBREGION_NAME='Scotland'
    SUBREGION_LIMITS=np.array([55,-6,59,-2])
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
if (OUTPUT_DIR[-1]=='/')&(not os.path.isdir(OUTPUT_DIR)):
      os.mkdir(OUTPUT_DIR)

#Print out option for repeat run in non-interactive mode:
print INTERACTIVE+' '+iDISPLAY+' '+str(SOURCES)+' '+PLOTS+' '+'-plot_tag '+plot_tag 

# Read grid file
print 'Reading gridfile: '+GRID_FILE
grinf=nc.Dataset(GRID_FILE,'r')
lats_2D   =grinf.variables['lats_2D'][:]
lons_2D   =grinf.variables['lons_2D'][:]
grindex   =grinf.variables['index_2D'][:]   # minus 1 to convert from Fortran to C index
grimask   =np.ones_like(grindex)
#flandex   =inf.variables['index_to1D'][:] 
grinf.close()


CS_DATA   = []
TIME_DATA = []
FRAC_DATA = []
LIT_C_DATA = []
RESP_S_DATA = []
C_VEG_DATA = []
print 'Reading in CS data'
for iSOURCE in SOURCES:
    iS = int(iSOURCE)
    print JULES_sources[iS][0]
    print 'file: '+JULES_sources[iS][1]
    inf = nc.Dataset(JULES_sources[iS][1],'r')
    
    if (iSOURCE==SOURCES[0]):
        CSlats=inf.variables['latitude'][:].squeeze()
        CSlons=inf.variables['longitude'][:].squeeze()
        # Calculate limits and range
        lon_limits= [np.floor(CSlons.min()),np.ceil(CSlons.max())]
        lon_range = lon_limits[1]-lon_limits[0]
        lat_limits= [np.floor(CSlats.min()),np.ceil(CSlats.max())]
        lat_range = lat_limits[1]-lat_limits[0]


    cs    = inf.variables['cs'][:].squeeze()
    frac  = inf.variables['frac'][:].squeeze()
    time  = nc.num2date(inf.variables['time'][:],\
            units=inf.variables['time'].units,       \
            calendar='standard'                      )

    CS_DATA.append( cs ) 
    TIME_DATA.append( time )
    FRAC_DATA.append( frac )
    
    # if plotting litterfall extract data
    if (PLOTS[11]=='Y') | (PLOTS[12]=='Y'): 
        lit_c = inf.variables['lit_c'][:].squeeze()
        # Apply fractional cover weight
        lit_c = lit_c * frac[:,:5,:]
        LIT_C_DATA.append(lit_c)
   
    # if plotting soil respiration extract data
    if (PLOTS[13]=='Y') | (PLOTS[14]=='Y') | (PLOTS[17]=='Y') | (PLOTS[18]=='Y'):
        resp_s = inf.variables['resp_s'][:].squeeze()
        RESP_S_DATA.append(resp_s* RESP_S_conv_factor)

    # if plotting c_veg extract data
    if (PLOTS[15]=='Y') | (PLOTS[16]=='Y'):
        c_veg = inf.variables['c_veg'][:].squeeze()
        c_veg = c_veg * frac[:,:5,:]
        C_VEG_DATA.append(c_veg)
    inf.close()

# full CHESS grid time series plots
if (PLOTS[0]=='Y'):
    print 'Producing time-series for full CHESS grid'
    # Create figure and plot axis if plotting CHESS timeseries
    FIG = plt.figure(figsize=[18,12])
    nplts_wdth   = 2.
    nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)
        
    for iPOOL in range(nPOOLs):
        print CS_pools[iPOOL]+':'
        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
        AX.set_title(CS_pools[iPOOL])
        for iSOURCE in SOURCES:

            iS  = int(iSOURCE)
            iCS = SOURCES.index(iSOURCE)
            loop_data=np.ma.masked_invalid(CS_DATA[iCS][:,iPOOL,:])
            
            CS_timeseries =  np.mean(loop_data,axis=1)

            print JULES_sources[iS][0]+' mean = ',np.mean(CS_timeseries)        
            AX.plot(TIME_DATA[iCS],CS_timeseries,label=JULES_sources[iS][0],
                    color=JULES_sources[iS][2],lw=2)
        
        AX.set_ybound(lower=0,upper=CS_maxes[iPOOL])
        #AX.set_yscale('log')

    AX.legend( bbox_to_anchor=(-0.1,-0.11),loc=10,borderaxespad=0.,ncol=min(len(SOURCES),6) )
    FIG.text(0.03,0.5,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
    FIG.tight_layout(rect=[0.05,0.05,1.0,0.95])
    FIG.suptitle('Time-series of the 4 Soil Carbon Pools',fontsize=30)
    
    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'SoilCarbon_Timeseries.png')
        plt.close()

# Timeseries by Land Class 
if (PLOTS[1]=='Y'):
    # Loop round each CS pool
    #iPOOL=3
    for iPOOL in range(nPOOLs):
        FIG=plt.figure(figsize=[18,12])
        nplts_wdth = 3
        nplts_hght = 2
        
        for iLC in range(len(LC_names)):        
            LC_name = LC_names[iLC]
            # Create axis for transom region
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iLC+1)
            AX.set_title(LC_name)
            

            for iSOURCE in SOURCES:
                iS  = int(iSOURCE)
                iCS = SOURCES.index(iSOURCE)

                LCmask=FRAC_DATA[iCS][:,LC_indexes[iLC],:]<LC_threshold
                data_iLC=np.ma.masked_array(CS_DATA[iCS][:,iPOOL,:],mask=LCmask)
                CS_timeseries =  np.mean(data_iLC,axis=1)
                
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


# Plot maps of mean soil carbon pools for final year for full CHESS grid
if (PLOTS[2]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
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
        
            PT.plot_map(CS_mapdata,lons_2D,lats_2D,                      \
                        DATA_RANGE=[0,maxval],                           \
                        PLOT_TITLE=CS_pools[iPOOL],                      \
                        COLOURS=['white','orange','brown'],              \
                        INTERPOLATE_COLOURS=True,NLEVELS=9,              \
                        CBAR_LABEL='Soil Carbon '+CS_units,              \
                        TICK_FORMAT=CS_tickformats[iPOOL],               \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5), \
                        LON_RANGE=lon_limits,LAT_RANGE=lat_limits,       \
                        RESOLUTION='i',FONTSIZES=[12,12,15,18],          \
                        AXIS=AX,                                         )

        
        FIG.tight_layout(rect=[0.02,0.05,1.0,0.92],h_pad=6)
        FIG.suptitle( 'UK-CHESS Soil Carbon Map, '+JULES_sources[iS][0].replace('_','-'),fontsize=25*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'CHESS-UK_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of mean total soil carbon for final year for full CHESS grid
if (PLOTS[3]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        #AX  = FIG.add_subplot(1,1,1)
        CS_mapdata = np.sum(CS_DATA[iCS][-12:,:,:],axis=1)
        CS_mapdata =  np.mean(CS_mapdata,axis=0) 
        # convert to 2d
        CS_mapdata = CS_mapdata[grindex]*grimask
        maxval = CS_maxes[4] 
        
        PT.plot_map(CS_mapdata,lons_2D,lats_2D,                       \
                    DATA_RANGE=[0,maxval],                            \
                    COLOURS=['white','orange','brown'],               \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,               \
                    CBAR_LABEL='Soil Carbon '+CS_units,               \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=CS_tickformats[4],CBAR_TICK_LENGTH=20,\
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Total Soil Carbon, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Total_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of mean total soil carbon difference for first 2 sources
if (PLOTS[4]=='Y'):
    
    iS1  = int(SOURCES[0])
    iCS1 = 0
    
    iS2  = int(SOURCES[1])
    iCS2 = 1

    xy_ratio = (lon_range/lat_range) 
    FIG = plt.figure(figsize=[int(12*xy_ratio),12])
        
    # subtract source-2 from source-1 (For last 12 months)
    #denom = (CS_DATA[iCS1][-12:,:,:]+CS_DATA[iCS2][-12:,:,:])/2
    CS_mapdata = (CS_DATA[iCS1][-12:,:,:]-CS_DATA[iCS2][-12:,:,:])
    # Sum over the 4 pools (axis=1)
    CS_mapdata = np.sum(CS_mapdata,axis=1)
    # Mean ove the months (axis=0)
    CS_mapdata = np.mean(CS_mapdata,axis=0)

    
    # convert to 2d 
    CS_mapdata = CS_mapdata[grindex]*grimask

    urban_mask = FRAC_DATA[iCS1][-1,5,grindex]*grimask
    urban_mask = np.ma.masked_less(urban_mask,0.5)

    maxval = CS_maxes[4] 
    
    #subtitle = JULES_sources[iS1][0].replace('_','-')+' minus '+\
    #           JULES_sources[iS2][0].replace('_','-')
    
    COLOURS  = [JULES_sources[iS2][2],'white',JULES_sources[iS1][2]]
    PT.plot_map(CS_mapdata,lons_2D,lats_2D,                       \
                GREYMASK=urban_mask,\
                DATA_RANGE=CS_DIFF_RANGE, COLOURS=COLOURS,        \
                INTERPOLATE_COLOURS=True,NLEVELS=200,NTICKS=11,   \
                CBAR_LABEL='Soil Carbon '+CS_units,               \
                CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                TICK_FORMAT=CS_tickformats[4],CBAR_TICK_LENGTH=20,\
                LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                RIVERS=True,FIGURE=FIG                            )
      
    FIG.tight_layout(rect=[0.01,0.02,1.0,0.94])
    FIG.suptitle( 'UK-CHESS Total Soil Carbon Difference',\
                   fontsize=28*xy_ratio)
    FIG.text(0.03,0.03,JULES_sources[iS2][3],ha='left',fontsize=17)
    FIG.text(0.97,0.03,JULES_sources[iS1][3],ha='right',fontsize=17)

    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'UK_Total_Soil_Carbon_Difference_Map_'+JULES_sources[iS1][0]+\
                                '_'+JULES_sources[iS2][0]+'.png')
        plt.close()


# Plot sub-region maps of mean soil carbon for final year
if (PLOTS[5]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # caluculate figure size based on lat_range/lon_range ratio
        sub_lon_range=(SUBREGION_LIMITS[1]-SUBREGION_LIMITS[3])
        sub_lat_range=(SUBREGION_LIMITS[0]-SUBREGION_LIMITS[2])
        xy_ratio = (sub_lon_range/sub_lat_range) 
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
        
            PT.plot_map(CS_mapdata,lons_2D,lats_2D,                            \
                        DATA_RANGE=[0,maxval],                                 \
                        PLOT_TITLE=CS_pools[iPOOL],                            \
                        COLOURS=['white','orange','brown'],                    \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,                   \
                        CBAR_LABEL='Soil Carbon '+CS_units,                    \
                        TICK_FORMAT=CS_tickformats[iPOOL],                     \
                        RESOLUTION='i',FONTSIZES=[12,12,15,18],                \
                        LON_RANGE=[SUBREGION_LIMITS[1],SUBREGION_LIMITS[3]],   \
                        LAT_RANGE=[SUBREGION_LIMITS[0],SUBREGION_LIMITS[2]],   \
                        LONDEL=int(sub_lon_range/5),LATDEL=int(sub_lat_range/5),\
                        AXIS=AX,     )

        
        
        FIG.tight_layout(rect=[0.02,0.05,1.0,0.92],h_pad=6)
        FIG.suptitle( SUBREGION_NAME+' Soil Carbon Map, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=30*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+SUBREGION_NAME+'_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot sub-region maps of mean total soil carbon for final year
if (PLOTS[6]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # caluculate figure size based on lat_range/lon_range ratio
        sub_lon_range=(SUBREGION_LIMITS[1]-SUBREGION_LIMITS[3])
        sub_lat_range=(SUBREGION_LIMITS[0]-SUBREGION_LIMITS[2])
        xy_ratio = (sub_lon_range/sub_lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
        # Create figure and plot axis if plotting global timeseries

        CS_mapdata = np.sum(CS_DATA[iCS][-12:,:,:],axis=1)
        CS_mapdata = np.mean(CS_mapdata,axis=0) 
        # convert to 2d
        CS_mapdata = CS_mapdata[grindex]*grimask
        maxval = CS_maxes[4] 
        
        PT.plot_map(CS_mapdata,lons_2D,lats_2D,                            \
                    DATA_RANGE=[0,maxval],                                 \
                    COLOURS=['white','orange','brown'],                    \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,                   \
                    CBAR_LABEL='Soil Carbon '+CS_units,                    \
                    TICK_FORMAT=CS_tickformats[4],                     \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],                \
                    LON_RANGE=[SUBREGION_LIMITS[1],SUBREGION_LIMITS[3]],   \
                    LAT_RANGE=[SUBREGION_LIMITS[0],SUBREGION_LIMITS[2]],   \
                    LONDEL=max(int(lon_range/5.),0.5),LATDEL=max(int(lat_range/5.),0.5),       \
                    RIVERS=True,FIGURE=FIG                            )

        
        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Total Soil Carbon, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+SUBREGION_NAME+'_Soil_Carbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of land cover
if (PLOTS[10]=='Y'):
    
    iS  = int(SOURCES[0])
    iCS = 0
    
    frac_data = FRAC_DATA[iCS][-1,:]
    dom_frac  = np.argmax(frac_data,axis=0)
    dom_frac = dom_frac[grindex]*grimask
    
    xy_ratio = (lon_range/lat_range) 
    FIG = plt.figure(figsize=[int(12*xy_ratio),12])
    PFT_names = ['Broadleaf','Needleleaf','C3 Grass','Shrub','Crop','Urban','Lake','Bare Soil']
    COLOURS  = ['olive','darkgreen','gold','orange','palegreen','red','blue','saddlebrown']
    PT.plot_map(dom_frac,lons_2D,lats_2D,                       \
                COLOURS=COLOURS,CLEVELS=range(0,9),DATA_RANGE=[0,9],  \
                CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                TICK_FORMAT='',CBAR_TICK_LENGTH=20,\
                LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                RIVERS=True,FIGURE=FIG                            )
      
    FIG.tight_layout(rect=[0.02,0.02,1.0,0.95])
    FIG.suptitle( 'UK-CHESS Dominant Land Cover Map',\
                   fontsize=28*xy_ratio)
    for i in range(len(PFT_names)):
        FIG.text(0.06+(i*0.112),0.015,PFT_names[i])

    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'UK_Dominant_LandCover_Map.png')
        plt.close()


# Plot maps of land cover
if (PLOTS[10]=='Y'):
    
    iS  = int(SOURCES[0])
    iCS = 0
    
    frac_data = FRAC_DATA[iCS][-1,:]
    dom_frac  = np.argmax(frac_data,axis=0)
    dom_frac = dom_frac[grindex]*grimask
    #plt.imshow(dom_frac,origin='bottom')
    #plt.colorbar()
    #plt.show()
    
    xy_ratio = (lon_range/lat_range) 
    FIG = plt.figure(figsize=[int(12*xy_ratio),12])
    PFT_names = ['Broadleaf','Needleleaf','C3 Grass','Shrub','Crop','Urban','Lake','Bare Soil']
    COLOURS  = ['olive','darkgreen','gold','orange','palegreen','red','blue','saddlebrown']
    PT.plot_map(dom_frac,lons_2D,lats_2D,                       \
                COLOURS=COLOURS,CLEVELS=range(0,9),DATA_RANGE=[0,9],  \
                CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                TICK_FORMAT='',CBAR_TICK_LENGTH=20,\
                LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                RIVERS=True,FIGURE=FIG                            )
      
    FIG.tight_layout(rect=[0.02,0.02,1.0,0.95])
    FIG.suptitle( 'UK-CHESS Dominant Land Cover Map',\
                   fontsize=28*xy_ratio)
    for i in range(len(PFT_names)):
        FIG.text(0.06+(i*0.112),0.015,PFT_names[i])

    if (iDISPLAY=='Y'):
        plt.show()
    else:
        FIG.savefig(OUTPUT_DIR+'UK_Dominant_LandCover_Map.png')
        plt.close()


# Plot maps of total litter fall
if (PLOTS[11]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        #AX  = FIG.add_subplot(1,1,1)
        lit_c_mapdata = np.sum(LIT_C_DATA[iCS][:,:,:],axis=1)
        lit_c_mapdata = np.mean(lit_c_mapdata,axis=0) 
        # convert to 2d
        lit_c_mapdata = np.ma.masked_array(lit_c_mapdata[grindex],mask=grindex.mask)
        maxval = LIT_C_max
        
        PT.plot_map(lit_c_mapdata,lons_2D,lats_2D,                    \
                    DATA_RANGE=[0,maxval],                            \
                    COLOURS=['white','moccasin','green'],             \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                    CBAR_LABEL='Litterfall '+LIT_C_units,             \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=LIT_C_tickformat,CBAR_TICK_LENGTH=20, \
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Litterfall Rate, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_Litterfall_Rate_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of total litter fall by PFT
if (PLOTS[12]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 3.
        nplts_hght = np.ceil(nPFTs/nplts_wdth)
        cr_ratio = (nplts_wdth/nplts_hght)
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(15*xy_ratio*cr_ratio),15])
        
        for iPFT in range(nPFTs):
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPFT+1)
            lit_c_mapdata = LIT_C_DATA[iCS][:,iPFT,:]
            lit_c_mapdata = np.mean(lit_c_mapdata,axis=0) 
            # convert to 2d
            lit_c_mapdata = np.ma.masked_array(lit_c_mapdata[grindex],mask=grindex.mask)
            maxval = LIT_C_max/2.
            
            PT.plot_map(lit_c_mapdata,lons_2D,lats_2D,                    \
                        DATA_RANGE=[0,maxval],                            \
                        PLOT_TITLE=PFT_names[iPFT],                       \
                        COLOURS=['white','moccasin','green'],             \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                        CBAR_LABEL='Litterfall '+LIT_C_units,             \
                        CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                        TICK_FORMAT=LIT_C_tickformat,                     \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                        LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                        RESOLUTION='i',FONTSIZES=[12,12,15,18],           \
                        AXIS=AX                                            )


        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Litterfall Rate by PFT, '+JULES_sources[iS][0].replace('_','-'),\
                       fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_Litterfall_Rate_Map_ByPFT_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of total soil respiration fall
if (PLOTS[13]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        #AX  = FIG.add_subplot(1,1,1)
        resp_s_mapdata = np.sum(RESP_S_DATA[iCS][:,:,:],axis=1)
        resp_s_mapdata = np.mean(resp_s_mapdata,axis=0) 
        # convert to 2d
        resp_s_mapdata = np.ma.masked_array(resp_s_mapdata[grindex],mask=grindex.mask)
        maxval = RESP_S_maxes[-1]
        
        PT.plot_map(resp_s_mapdata,lons_2D,lats_2D,                   \
                    DATA_RANGE=[0,maxval],                            \
                    COLOURS=['white','moccasin','orangered'],            \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                    CBAR_LABEL='Soil Respiration ('+RESP_S_units+')', \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=RESP_S_tickformats[-1],CBAR_TICK_LENGTH=20,\
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Soil Respiration Rate, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_SoilResp_Rate_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of soil respiration by pool
if (PLOTS[14]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 2.
        nplts_hght = np.ceil(nPOOLs/nplts_wdth)
        cr_ratio = (nplts_wdth/nplts_hght)
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(15*xy_ratio*cr_ratio),15])
        
        for iPOOL in range(nPOOLs):
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
            resp_s_mapdata = RESP_S_DATA[iCS][:,iPOOL,:]
            resp_s_mapdata = np.mean(resp_s_mapdata,axis=0) 
            # convert to 2d
            resp_s_mapdata = np.ma.masked_array(resp_s_mapdata[grindex],mask=grindex.mask)
            maxval = RESP_S_maxes[iPOOL]/2
        
            
            PT.plot_map(resp_s_mapdata,lons_2D,lats_2D,                   \
                        DATA_RANGE=[0,maxval],                            \
                        PLOT_TITLE=CS_pools[iPOOL],                       \
                        COLOURS=['white','moccasin','orangered'],         \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                        CBAR_LABEL='Soil Respiration ('+RESP_S_units+')', \
                        CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                        TICK_FORMAT=RESP_S_tickformats[iPOOL],            \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                        LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                        RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                        AXIS=AX                                           )

        
        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Soil Respiration Rate by Pool, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_SoilResp_Rate_Map_ByPool_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of total vegetation carbon
if (PLOTS[15]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        #AX  = FIG.add_subplot(1,1,1)
        c_veg_mapdata = np.sum(C_VEG_DATA[iCS][:,:,:],axis=1)
        c_veg_mapdata = np.mean(c_veg_mapdata,axis=0) 
        # convert to 2d
        c_veg_mapdata = np.ma.masked_array(c_veg_mapdata[grindex],mask=grindex.mask)
        maxval = C_VEG_maxes[5]
        
        PT.plot_map(c_veg_mapdata,lons_2D,lats_2D,                    \
                    DATA_RANGE=[0,maxval],                            \
                    COLOURS=['white','moccasin','goldenrod'],         \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                    CBAR_LABEL='Vegetation Carbon ('+C_VEG_units+')', \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=C_VEG_tickformat,CBAR_TICK_LENGTH=20, \
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Vegetation Carbon, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_VegetationCarbon_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of total litter fall by PFT
if (PLOTS[16]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 3.
        nplts_hght = np.ceil(nPFTs/nplts_wdth)
        cr_ratio = (nplts_wdth/nplts_hght)
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(15*xy_ratio*cr_ratio),15])
        
        for iPFT in range(nPFTs):
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPFT+1)
            c_veg_mapdata = C_VEG_DATA[iCS][:,iPFT,:]
            c_veg_mapdata = np.max(c_veg_mapdata,axis=0) 
            # convert to 2d
            c_veg_mapdata = np.ma.masked_array(c_veg_mapdata[grindex],mask=grindex.mask)
            maxval = C_VEG_maxes[iPFT]
            PT.plot_map(c_veg_mapdata,lons_2D,lats_2D,                    \
                        DATA_RANGE=[0,maxval],                            \
                        PLOT_TITLE=PFT_names[iPFT],                       \
                        COLOURS=['white','moccasin','goldenrod'],         \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                        CBAR_LABEL='Vegetation Carbon ('+C_VEG_units+')', \
                        CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                        TICK_FORMAT=C_VEG_tickformat,                     \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                        LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                        RESOLUTION='i',FONTSIZES=[12,12,15,18],           \
                        AXIS=AX                                            )


        AX  = FIG.add_subplot(nplts_hght,nplts_wdth,nplts_hght*nplts_wdth)
        c_veg_mapdata = np.sum(C_VEG_DATA[iCS][:,:,:],axis=1)
        c_veg_mapdata = np.max(c_veg_mapdata,axis=0) 
        # convert to 2d
        c_veg_mapdata = np.ma.masked_array(c_veg_mapdata[grindex],mask=grindex.mask)
        maxval = C_VEG_maxes[5]
        
        PT.plot_map(c_veg_mapdata,lons_2D,lats_2D,                    \
                    DATA_RANGE=[0,maxval],                            \
                    PLOT_TITLE='Total',                               \
                    COLOURS=['white','moccasin','goldenrod'],         \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                    CBAR_LABEL='Vegetation Carbon ('+C_VEG_units+')', \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=C_VEG_tickformat,CBAR_TICK_LENGTH=20, \
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    AXIS=AX                                           )

        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Max Vegetation Carbon by PFT, '+JULES_sources[iS][0].replace('_','-'),\
                       fontsize=28*xy_ratio)
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Max_VegetationCarbon_Map_ByPFT_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of total soil respiration factor
if (PLOTS[17]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        #AX  = FIG.add_subplot(1,1,1)
        resp_fact_mapdata = np.sum(RESP_S_DATA[iCS][:,:,:]/CS_DATA[iCS][:,:,:],axis=1)
        resp_fact_mapdata = np.mean(resp_fact_mapdata,axis=0) 
        # convert to 2d
        resp_fact_mapdata = np.ma.masked_array(resp_fact_mapdata[grindex],mask=grindex.mask)
        maxval = RESP_FRAC_maxes[-1]
        
        urban_mask = FRAC_DATA[iCS][-1,5,grindex]*grimask
        urban_mask = np.ma.masked_less(urban_mask,0.5)
        
        PT.plot_map(resp_fact_mapdata,lons_2D,lats_2D,                \
                    DATA_RANGE=[0,maxval],                            \
                    GREYMASK=urban_mask,\
                    COLOURS=['white','moccasin','orangered'],         \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,              \
                    CBAR_LABEL='Soil Respiration Factor ('+RESP_FACT_units+')', \
                    CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                    TICK_FORMAT=RESP_FACT_tickformats[-1],CBAR_TICK_LENGTH=20,\
                    LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Soil Respiration Factor, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_SoilResp_Factor_Map_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot maps of soil respiration by pool
if (PLOTS[18]=='Y'):
    
    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        # Create figure and plot axis if plotting global timeseries
        nplts_wdth   = 2.
        nplts_hght = np.ceil(nPOOLs/nplts_wdth)
        cr_ratio = (nplts_wdth/nplts_hght)
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(15*xy_ratio*cr_ratio),15])
        
        urban_mask = FRAC_DATA[iCS][-1,5,grindex]*grimask
        urban_mask = np.ma.masked_less(urban_mask,0.5)

        for iPOOL in range(nPOOLs):
            AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
            resp_fact_mapdata = RESP_S_DATA[iCS][:,iPOOL,:]/CS_DATA[iCS][:,iPOOL,:]
            resp_fact_mapdata = np.max(resp_fact_mapdata,axis=0) 
            # convert to 2d
            resp_fact_mapdata = np.ma.masked_array(resp_fact_mapdata[grindex],mask=grindex.mask)
            maxval = RESP_FRAC_maxes[iPOOL] 
        
            PT.plot_map(resp_fact_mapdata,lons_2D,lats_2D,                \
                        DATA_RANGE=[0,maxval],                            \
                        GREYMASK=urban_mask,\
                        PLOT_TITLE=CS_pools[iPOOL],                       \
                        COLOURS=['white','moccasin','orangered'],         \
                        INTERPOLATE_COLOURS=True,NLEVELS=RESP_FACT_NLEVELS[iPOOL],  \
                        CBAR_LABEL='Soil Respiration Factor ('+RESP_FACT_units+')', \
                        CBAR_SIZE='4%',CBAR_PAD=0.4,                      \
                        TICK_FORMAT=RESP_FACT_tickformats[iPOOL],         \
                        LONDEL=int(lon_range/5),LATDEL=int(lat_range/5),  \
                        LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                        RESOLUTION='i',FONTSIZES=[15,15,18,21],           \
                        AXIS=AX                                           )

        
        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Mean Soil Respiration Factor by Pool, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            FIG.savefig(OUTPUT_DIR+'UK_Mean_SoilRespFactor_Map_ByPool_'+JULES_sources[iS][0]+'.png')
            plt.close()


# Plot Maps of (kappa_dpm*C_dpm+kappa_rpm*C_rpm)/(kappa_bio*C_bio+kappa_hum*C_hum)
if (PLOTS[19]=='Y'):

    for iSOURCE in SOURCES:
        iS  = int(iSOURCE)
        iCS = SOURCES.index(iSOURCE)
    
        xy_ratio = (lon_range/lat_range) 
        FIG = plt.figure(figsize=[int(12*xy_ratio),12])
            
        DPM_kapC = CS_DATA[iCS][-12:,0,:]*kappa_s_RothC[0]
        RPM_kapC = CS_DATA[iCS][-12:,1,:]*kappa_s_RothC[1]
        BIO_kapC = CS_DATA[iCS][-12:,2,:]*kappa_s_RothC[2]
        HUM_kapC = CS_DATA[iCS][-12:,3,:]*kappa_s_RothC[3]

        CS_mapdata = (DPM_kapC+RPM_kapC)/(BIO_kapC+HUM_kapC) 
        CS_mapdata = np.mean(CS_mapdata,axis=0) 
        # convert to 2d
        CS_mapdata = CS_mapdata[grindex]*grimask
        maxval = CS_maxes[4]
        print(np.mean(CS_mapdata[CS_mapdata.mask==False]),\
                np.std(CS_mapdata[CS_mapdata.mask==False]),\
                np.median(CS_mapdata[CS_mapdata.mask==False]))
        PT.plot_map(CS_mapdata,lons_2D,lats_2D,                            \
                    DATA_RANGE=[0,10],\
                    COLOURS=['white','yellow','orange','brown'],                    \
                    INTERPOLATE_COLOURS=True,NLEVELS=200,NTICKS=11,        \
                    CBAR_LABEL='Factor ',                    \
                    TICK_FORMAT=CS_tickformats[4],                     \
                    RESOLUTION='i',FONTSIZES=[15,15,18,21],                \
                    LON_RANGE=lon_limits,LAT_RANGE=lat_limits,        \
                    LONDEL=max(int(lon_range/5.),0.5),LATDEL=max(int(lat_range/5.),0.5),       \
                    RIVERS=True,FIGURE=FIG                            )

        
        FIG.tight_layout(rect=[0.01,0.02,1.0,0.95])
        FIG.suptitle( 'UK-CHESS Eleanor Factor, '+JULES_sources[iS][0].replace('_','-'),\
                        fontsize=28*xy_ratio)
        
        if (iDISPLAY=='Y'):
            plt.show()
        else:
            print(OUTPUT_DIR+'Eleanor_Factor_'+JULES_sources[iS][0]+'.png')
            FIG.savefig(OUTPUT_DIR+'Eleanor_Factor_'+JULES_sources[iS][0]+'.png')
            plt.close()


