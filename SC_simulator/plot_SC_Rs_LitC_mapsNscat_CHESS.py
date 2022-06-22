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
import scipy.stats as stats
import netCDF4 as nc
import netcdftime as nctime
import plot_tools as PTs
import sys, os

##########################################
# Set Universal parameters and variables
nYEARs=33
nMONTHs=12
nlandpoints=67209
nPOOLs=4
nPFTs=5
PFT_names      = [ 'BL','NL','C3','C4','Sh']
PFT_longnames  = [ 'Broadleaf','Needleleaf','C3 Grass','C4 Grass','Shrub' ]
PFT_colours    = [ 'lime','darkgreen', 'yellow','yellowgreen','orange' ]
pool_names     = ['DPM','RPM','BIO','HUM']
pool_longnames = [ 'Decomposable Plant Material', \
                   'Resistant Plant Material',    \
                   'Microbial Biomass',           \
                   'Hummus'                       ]

##################################################
# Select Investigation Region
print 'Available Regions: '
REGIONS = ['Global','UK','UK-fixed']
for i in range(len(REGIONS)):
    print str(i)+': '+REGIONS[i]

iREGION=input('Select a region: \n')
while iREGION>=len(REGIONS):
    iREGION=input('Selected region does not exist, please select another region: \n')


############################################################
# Set eEgionally dependent Variables
REG_lat_lims = [ (-90,90), (49,61), (49,61) ]
REG_lon_lims = [ (-180,180), (-11,2.5), (-11,2.5) ]
REG_resolutions = [ 'c', 'h', 'h' ]
REG_textpos  = [ (-175,-50,-15), (-10.5,60,-0.8), (-10.5,60,-0.8)  ]
REG_gridspc  = [ (30,15), (2.5,2.5), (2.5,2.5) ]
REG_nlevels  = [ 10, 8, 8 ]
REG_CSlims   = [ (0,15.), None, (0,5) ]
REG_LitClims = [ (0,1.), None, (0,0.5) ]
REG_Rslims   = [ (0,15), None, (0,5) ]


region=REGIONS[iREGION]
lat_limits=REG_lat_lims[iREGION]
lat_indices = [ (ll+90)*2 for ll in lat_limits ]
lon_limits=REG_lon_lims[iREGION]
lon_indices = [ (ll+180)*2 for ll in lon_limits ]

resolution=REG_resolutions[iREGION]
textpos = REG_textpos[iREGION]
gridspc = REG_gridspc[iREGION]
nlevels = REG_nlevels[iREGION]
CSlims  = REG_CSlims[iREGION]
LitClims= REG_LitClims[iREGION]
Rslims  = REG_Rslims[iREGION]

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

########################################################################
print 'Setting Filenames and Directories'
# Input files:
#WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
CHESS_dir        = '/users/eow/edwcom/CHESS/'
gridfile         = CHESS_dir+'chess_jules_land_index.nc'
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/'
SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'

plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/CHESS/'+region+'/CS/'
plots_DIR = plots_DIR.replace('-fixed','')
os.system('mkdir -p '+plots_DIR)

SoilCarbon_Sources = [ 'SCsim CHESS JULES 4.1 Run' ]                          

SoilCarbon_files  = [ SCsimulator_DIR+'SCsim_J4.1_CHESSrun_FastSoilCarbon.nc' ] 

SoilCarbon_plottags= [ 'SCsim_CHESSrun'  ]

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

######################################################################
print 'Read in data'
# Open and read the land index file
print 'Grid file: ',gridfile
inf=nc.Dataset(gridfile,'r')
lats_2D   =inf.variables['lats_2D'][:]
lons_2D   =inf.variables['lons_2D'][:]
grindex=inf.variables['index_2D'][:]   # minus 1 to convert from Fortran to C index
grimask=np.ones_like(grindex)
inf.close()

plot_lons,plot_lats=lons_2D,lats_2D

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

# Get ICE mask from sthu in JULES BS file
#print 'Getting ICEmask from: '+SoilCarbon_files[0]
#ICEmask=np.ones_like( \
#        np.ma.masked_equal( \
#            nc.Dataset(SoilCarbon_files[0],'r').variables['sthu'][0,0,:].squeeze() \
#          , 0.0 ) )

# Open and read in Soil Carbon source data
print 'Soil Carbon file: '+SoilCarbon_files[iSource]
inf=nc.Dataset(SoilCarbon_files[iSource],'r')

if (SoilCarbon_Sources[iSource][:5]=='SCsim'):
    SC_data_1D   = inf.variables['C_4pools_'+SoilRespFunctions[iSRF]][:].squeeze()
    LitC_data_1D = inf.variables['Total_Litterfall_C'][:].squeeze() 
    Rs_data_1D   = inf.variables['Soil_Resp_Fact_'+SoilRespFunctions[iSRF]][:].squeeze() 

    if ('limatology' in SoilCarbon_Sources[iSource]):
        SC_data_1D = np.mean(SC_data_1D,axis=1)
    
    if not ('tatic' in SoilCarbon_Sources[iSource]):
        LitC_data_1D = np.mean(LitC_data_1D,axis=0)

    Rs_data_1D=np.mean(Rs_data_1D,axis=1)*3600.*24.*365.

elif (SoilCarbon_Sources[iSource][:5]=='JULES'):
    SC_data_1D = np.mean(inf.variables['cs'][:].squeeze(),axis=0)
    LitC_data_1D = np.sum( np.mean(  inf.variables['lit_c'][:].squeeze() \
                                   * inf.variables['frac'][:,:5,:,:].squeeze() \
                                   ,  axis=0 )  , axis=0 )
    Rs_data_1D = np.mean(inf.variables['resp_s'][:].squeeze(), axis=0) \
                * 3600.*24.*365./SC_data_1D


inf.close()

SC_data_2D       = SC_data_1D[:,grindex]*grimask
total_SC_data_2D = np.sum(SC_data_2D,axis=0)
plot_total_SC    = total_SC_data_2D
                                  
LitC_data_2D     = LitC_data_1D[grindex]*grimask
plot_LitC        = LitC_data_2D
                              
Rs_data_2D       = Rs_data_1D[:,grindex]*grimask
total_Rs_data_2D = np.sum(Rs_data_2D,axis=0)
plot_total_Rs    = total_Rs_data_2D
                                  




if (CSlims==None):   CSlims  =[0,np.ceil(np.mean(plot_total_SC)*2.0)]
if (LitClims==None): LitClims=[0,np.ceil(np.mean(plot_LitC)*2.0)]
if (Rslims==None):   Rslims  =[0,np.ceil(np.mean(plot_total_Rs)*2.0)]


print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'
print LitC_data_1D.shape
print np.sum(SC_data_1D,axis=0).shape
print np.sum(Rs_data_1D,axis=0).shape

CS_LitC_dex=np.where( (LitC_data_1D!=0)&(np.sum(SC_data_1D,axis=0)!=0) )
CS_LitC_CORRELATION = stats.pearsonr(np.sum(SC_data_1D,axis=0)[CS_LitC_dex],\
                                            LitC_data_1D[CS_LitC_dex]    )
print 'total soil carbon/litter correlation = '\
        + str(CS_LitC_CORRELATION)

CS_Rs_dex=np.where( (np.sum(Rs_data_1D,axis=0)!=0)&(np.sum(SC_data_1D,axis=0)!=0) )
CS_Rs_CORRELATION = stats.pearsonr(np.sum(SC_data_1D,axis=0)[CS_Rs_dex],\
                                   np.sum(Rs_data_1D,axis=0)[CS_Rs_dex])
print 'total soil carbon/respiration correlation = '\
        + str(CS_Rs_CORRELATION)

#################################################################################
print 'Plotting Litterfall C '+region+' maps - '+SoilCarbon_Sources[iSource]
FIG = plt.figure(figsize=(12,12))

# Plot map of Soil Carbon
#AXIS = plt.subplot2grid( (3,3),(0,0),colspan=2)
AXIS = FIG.add_subplot(3,2,1)
PTs.plot_map(plot_total_SC, \
             plot_lons,plot_lats, \
             AXIS=AXIS, \
             DATA_RANGE=CSlims,    \
             LAT_RANGE=lat_limits, LON_RANGE=lon_limits, \
             LONDEL=gridspc[0], LATDEL=gridspc[1], \
             MPL_CBAR='YlOrBr', NLEVELS=nlevels, \
             CBAR_SIZE='7%',CBAR_PAD=0.2,\
             PLOT_TITLE='Total Soil Carbon', \
             CBAR_LABEL='$kg$C $m^{-2}$', \
             RESOLUTION=resolution)

# Plot map of Litterfall
#AXIS = plt.subplot2grid( (3,3),(1,0),colspan=2)
AXIS=FIG.add_subplot(3,2,3)
PTs.plot_map(plot_LitC, \
             plot_lons,plot_lats, \
             AXIS=AXIS, \
             DATA_RANGE=LitClims,    \
             LAT_RANGE=lat_limits, LON_RANGE=lon_limits, \
             LONDEL=gridspc[0], LATDEL=gridspc[1], \
             MPL_CBAR='YlOrBr', NLEVELS=nlevels, \
             CBAR_SIZE='7%',CBAR_PAD=0.2,\
             PLOT_TITLE='Total Litterfall', \
             CBAR_LABEL='$kg$C $m^{-2}$ $yr^{-1}$', \
             RESOLUTION=resolution)

#Plot Scatter of Litter against CS
#AXIS = plt.subplot2grid( (3,3),(1,2),colspan=1)
AXIS=FIG.add_subplot(3,2,4)
AXIS.plot(np.sum(SC_data_1D,axis=0)[CS_LitC_dex],\
          LitC_data_1D[CS_LitC_dex],\
          ls='',marker='.')
AXIS.set_ylabel('Total Litterfall $kg$C $m^{-2}$ $yr^{-1}$')
AXIS.set_ylim([0,0.4])#  [LitClim*1 for LitClim in LitClims])
AXIS.set_xlabel('Total Soil Carbon $kg$C $m^{-2}$')
AXIS.set_xlim([CSlim*4 for CSlim in CSlims])


# Plot map of Soil respiration function
#AXIS = plt.subplot2grid( (3,3),(2,0),colspan=2)
AXIS=FIG.add_subplot(3,2,5)
PTs.plot_map(plot_total_Rs, \
             plot_lons,plot_lats, \
             AXIS=AXIS, \
             DATA_RANGE=Rslims, \
             LAT_RANGE=lat_limits, LON_RANGE=lon_limits, \
             LONDEL=gridspc[0], LATDEL=gridspc[1], \
             MPL_CBAR='YlOrBr', NLEVELS=nlevels, \
             CBAR_SIZE='7%',CBAR_PAD=0.2,\
             PLOT_TITLE='Total Soil Respiration Function', \
             CBAR_LABEL='$yr^{-1}$', \
             RESOLUTION=resolution)

#Plot Scatter of Soil Resp against CS
#AXIS = plt.subplot2grid( (3,3),(2,2),colspan=1)
AXIS=FIG.add_subplot(3,2,6)
AXIS.plot(np.sum(SC_data_1D,axis=0)[CS_Rs_dex],\
          np.sum(Rs_data_1D,axis=0)[CS_Rs_dex],\
          ls='',marker='.')
AXIS.set_ylabel('Total Soil Respiration Function $yr^{-1}$')
AXIS.set_ylim([Rslim*4 for Rslim in Rslims])
AXIS.set_xlabel('Total Soil Carbon $kg$C $m^{-2}$')
AXIS.set_xlim([CSlim*4 for CSlim in CSlims])


#Create axis to output text to
#AXIS=FIG.add_subplot(3,2,2)
FIG.text(0.55,0.8, "CS to Litter correlation = "+str(round(CS_LitC_CORRELATION[0],2)),fontsize=20)
FIG.text(0.55,0.74, "CS to Soil Resp correlation = "+str(round(CS_Rs_CORRELATION[0],2)),fontsize=20)

FIG.subplots_adjust(hspace=0.3)
FIG.suptitle('Soil Carbon Breakdown - \n'+SoilCarbon_Sources[iSource],fontsize=24)
print 'Writing image to: '+plots_DIR+region+'_LitC_mapNscat_'+plottag+'.png'
FIG.savefig(plots_DIR+region+'_CS_breakdown_mapNscat_'+plottag+'.png', bbox_inches='tight')
FIG.clear()






