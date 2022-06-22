#!/usr/bin/env python2.7
################################################################################
# 
# Program: SC_simulator.py
# 
# Python Script to plot simulated soil carbon estimate 
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
################################################################################
#

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import plot_tools as PTs
import sys, os, glob


data_dir='/users/eow/edwcom/SC_simulator/output/'
plots_DIR=data_dir+'plots/'

data_strlen=len(data_dir)
files=glob.glob(data_dir+'SCsim*.nc')

for i in range(len(files)):
    print str(i)+': '+(files[i])[data_strlen:]

plottag=(files[i])[data_strlen:-3]

iFILE=input('Select a file: ')

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'
############################################################
# read in parsed arguments

REGIONS = ['Global','UK']
REG_lat_lims = [ (-90,90), (49,61) ]
REG_lon_lims = [ (-180,180), (-11,2.5) ]
REG_resolutions = [ 'c', 'h' ]
REG_gridspc  = [ (60,30), (2.5,2.5), (2.5,2.5) ]
REG_nlevels  = [ 10, 8, 8 ]

print 'Available Regions: '
for i in range(len(REGIONS)):
    print str(i)+': '+REGIONS[i]

iREGION=input('Select a region: \n')
while iREGION>=len(REGIONS):
    iREGION=input('Selected region does not exist, please select another region: \n')

region=REGIONS[iREGION]
lat_limits=REG_lat_lims[iREGION]
lon_limits=REG_lon_lims[iREGION]
resolution=REG_resolutions[iREGION]
gridspc = REG_gridspc[iREGION]
nlevels = REG_nlevels[iREGION]

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

##################################################
# Set parameters and variables
nYEARs=33
nMONTHs=12
nlandpoints=67209
nPOOLs=4
pool_names     = ['DPM','RPM','BIO','HUM']
pool_longnames = [ 'Decomposable Plant Material', \
                   'Resistant Plant Material',    \
                   'Microbial Biomass',           \
                   'Hummus'                       ]
pool_maxvals = [ 0.08, 6.0, 0.6, 24 ]

pool_CLEVELS = [ [0,0.01,0.02,0.04,0.6,0.08,0.1,0.12], \
                 [0,0.5,1,2,4,6,8,10,12], \
                 [0,0.01,0.02,0.04,0.08,0.12,0.16,0.22,0.30], \
                 [0,1.,2.,4.,8.,12,16.,20,24]  ]

select_trans   = [ 0, 1, 2, 3 ]


##################################################


SoilRespFunctions = ['Q10t','RothCit']
print 'Available Soil Respiration Functions: '
for i in range(len(SoilRespFunctions)):
    print str(i)+': '+SoilRespFunctions[i]

iSRF=input('Select an Soil Respiration Funtion: ')
  
while iSRF>=len(SoilRespFunctions):
    iSource=input('Selected Soil Respiration Funciton does not exist, please select another function: \n')

plottag=plottag+SoilRespFunctions[iSRF]

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

######################################################################
print 'Read in data'
# Open and read the land index file
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'
print 'Grid file: ',WFDEI_gridfile
inf=nc.Dataset(WFDEI_gridfile,'r')
lats   =inf.variables['latitude'][:]
lons   =inf.variables['longitude'][:]
grindex=inf.variables['land_index'][:]-1   # minus 1 to convert from Fortran to C index
grifrac=inf.variables['land_fraction'][:]
grimask=np.ones_like(grindex)
inf.close()
lons_2D,lats_2D=np.meshgrid(lons,lats)

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

# Open and read in Soil Carbon source data
inf=nc.Dataset(files[iFILE],'r')

SC_data_1D = inf.variables['C_4pools_'+SoilRespFunctions[iSRF]][:].squeeze()  

SC_data_2D=SC_data_1D[:,grindex]*grimask
   
print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'
#############################
print 'Plotting Soil Carbon '+region+' Map - '+plottag
FIG = plt.figure(figsize=(18,12))
for ipool in range(nPOOLs):
    AXIS=FIG.add_subplot(2,2,ipool+1)
    PTs.plot_map(SC_data_2D[ipool,:,:], \
                 lons_2D,lats_2D, \
                 AXIS=AXIS, \
                 LAT_RANGE=lat_limits, LON_RANGE=lon_limits, \
                 LONDEL=gridspc[0], LATDEL=gridspc[1], \
                 DATA_RANGE=[0,pool_maxvals[ipool]], \
                 MPL_CBAR='YlOrBr', NLEVELS=nlevels, \
                 CBAR_SIZE='8%',CBAR_PAD=0.2,\
                 PLOT_TITLE=pool_longnames[ipool], \
                 CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                 RESOLUTION=resolution)

#FIG.subplots_adjust(wspace=-0.2,hspace=0.3)
FIG.suptitle('Soil Carbon - ', fontsize=30 )
print 'Writing image to: '+plots_DIR+region+'_CS_map_'+plottag+'.png'
FIG.savefig(plots_DIR+region+'_CS_map_'+plottag+'.png', bbox_inches='tight')
FIG.clear()



