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
import sys, os


############################################################
# read in parsed arguments

REGIONS = ['Global','UK']
REG_lat_lims = [ (-90,90), (49,61) ]
REG_lon_lims = [ (-180,180), (-11,2.5) ]
REG_resolutions = [ 'c', 'h' ]
REG_gridspc  = [ (60,30), (2.5,2.5), (2.5,2.5) ]
REG_nlevels  = [ 10, 8, 8 ]

##################################################
# Select Investigation Region
print 'Available Regions: '
REGIONS = ['Global','UK']
for i in range(len(REGIONS)):
    print str(i)+': '+REGIONS[i]

#iREGION=input('Select a region: \n')
iREGION=1
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
pool_maxvals = [ 0.04, 3.0, 0.1, 5 ]

pool_CLEVELS = [ [0,0.01,0.02,0.04,0.6,0.08,0.1,0.12], \
                 [0,0.5,1,2,4,6,8,10,12], \
                 [0,0.01,0.02,0.04,0.08,0.12,0.16,0.22,0.30], \
                 [0,1.,2.,4.,8.,12,16.,20,24]  ]

select_trans   = [ 0, 1, 2, 3 ]



########################################################################
print 'Setting Filenames and Directories'
# Input files:
#WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
CHESS_dir        = '/users/eow/edwcom/CHESS/'
gridfile         = CHESS_dir+'chess_jules_land_index.nc'
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/'
SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'

plots_DIR          = '/users/eow/edwcom/SC_simulator/output/plots/WFDEI/'+region+'/CS/'
os.system('mkdir -p '+plots_DIR)

plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/CHESS/'+region+'/CS/'
plots_DIR = plots_DIR.replace('-fixed','')
os.system('mkdir -p '+plots_DIR)

SoilCarbon_Sources = [ 'SCsim CHESS JULES 4.1 Run' ,  \
                        'Ancillary Soil Carbon'       ]

SoilCarbon_files  = [ SCsimulator_DIR+'SCsim_J4.1_CHESSrun_FastSoilCarbon.nc',         \
                        '/prj/chess/data/1km/v1.0/ancil/chess_soilparams_hwsd_bc.nc'   ]

SoilCarbon_plottags= [ 'SCsim_CHESSrun', \
                        'Ancillary CS'   ]

print 'Available Soil Carbon Sources: '
for i in range(len(SoilCarbon_Sources)):
    print str(i)+': '+SoilCarbon_Sources[i]

iSource = input('Select a source: ')
while iSource>=len(SoilCarbon_Sources):
    iSource=input('Selected Soil Carbon Source does not exist, please select another source: \n')

plottag = SoilCarbon_plottags[iSource]

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

if (SoilCarbon_Sources[iSource][:5]=='SCsim'):
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
inf=nc.Dataset(gridfile,'r')
lats_2D   =inf.variables['lats_2D'][:]
lons_2D   =inf.variables['lons_2D'][:]
grindex=inf.variables['index_2D'][:]   # minus 1 to convert from Fortran to C index
grimask=np.ones_like(grindex)
flandex=inf.variables['index_to1D'][:] 
inf.close()

plot_lons,plot_lats=lons_2D,lats_2D

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

# Open and read in Soil Carbon source data
print 'Soil Carbon file: '+SoilCarbon_files[iSource]
inf=nc.Dataset(SoilCarbon_files[iSource],'r')

if (SoilCarbon_Sources[iSource][:5]=='SCsim'):
    SC_data_1D = inf.variables['C_4pools_'+SoilRespFunctions[iSRF]][:].squeeze()
    if ('limatology' in SoilCarbon_Sources[iSource]):
        SC_data_1D = np.mean(SC_data_1D,axis=1)
    SC_data_2D=SC_data_1D[:,grindex]*grimask
elif (SoilCarbon_Sources[iSource][:5]=='JULES'):
    SC_data_1D = np.mean(inf.variables['cs'][:].squeeze(),axis=0)
    SC_data_2D=SC_data_1D[:,grindex]*grimask
elif (SoilCarbon_Sources[iSource][:5]=='Ancil'):
    SC_data_2D = inf.variables['cs'][:].squeeze()
    print 'SC_data_2D.shape = ',SC_data_2D.shape
    #SC_data_1D = SC_data_2D[:,flandex[0,:],flandex[1,:] ]


        
print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'
#############################
print 'Plotting Soil Carbon '+region+' Map - '+SoilCarbon_Sources[iSource]
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
FIG.suptitle('Soil Carbon - '+SoilCarbon_Sources[iSource], fontsize=30 )
print 'Writing image to: '+plots_DIR+region+'_CS_map_'+plottag+'.png'
FIG.savefig(plots_DIR+region+'_CS_map_'+plottag+'.png', bbox_inches='tight')
FIG.clear()



