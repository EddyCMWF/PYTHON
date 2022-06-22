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

CSlims = (0,12)

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
pool_maxvals = [ 0.08, 6.0, 0.3, 12 ]

pool_CLEVELS = [ [0,0.01,0.02,0.04,0.6,0.08,0.1,0.12], \
                 [0,0.5,1,2,4,6,8,10,12], \
                 [0,0.01,0.02,0.04,0.08,0.12,0.16,0.22,0.30], \
                 [0,1.,2.,4.,8.,12,16.,20,24]  ]

select_trans   = [ 0, 1, 2, 3 ]



########################################################################
print 'Setting Filenames and Directories'
# Input files:
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/'
SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'

plots_DIR          = '/users/eow/edwcom/SC_simulator/output/plots/WFDEI/'+region+'/CS/'
os.system('mkdir -p '+plots_DIR)

trans_file   = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5_orig.nc'

SoilCarbon_Sources = [ 'JULES (100 year spin)',                           \
                       'JULES-Triffid Big Spin (1000 year)',              \
                       'JULES-Triffid Quick Spin (100 year)' ,            \
                       'SCsim Static LAI',                                \
                       'SCsim pheno LAI (100 year spin)',                 \
                       'SCsim pheno LAI (10 year spin)',                  \
                       'SCsim pheno LAI with TRIF params (10 year spin)', \
                       'SCsim BigSpin LAI',                               \
                       'SCsim Big Spin LAI (all timesteps)',              \
                       'SCsim Constant LAI',                              \
                       'SCsim pheno LAI (100 year spin) Climatology',     \
                       'HWSD-NCSCD Soil Carbon Map'                       ]

SoilCarbon_files  = [ JULES_output_dir+'WFDEI_pheno/JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno.monthly_mean.nc',                 \
                      JULES_output_dir+'WFDEI_BigSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc',     \
                      JULES_output_dir+'WFDEI_QuickSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_QuickSpin.monthly_mean.nc', \
                      SCsimulator_DIR+'SCsim_WL_JULES-WFDEI-Zinke-hydro1k_Static_LAI_FastSoilCarbon.nc',                  \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_phenoLAI_MarthewsTI_FastSoilCarbon.nc',                           \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_QSphenoLAI_MarthewsTI_FastSoilCarbon.nc',                         \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_QSphenoLAI_trifparam_MarthewsTI_FastSoilCarbon.nc',               \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_BigSpinLAI_MarthewsTI_FastSoilCarbon.nc',                         \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_BigSpinLAI_alltimesteps_MarthewsTI_FastSoilCarbon.nc',            \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_ConstLAI_alltimesteps_MarthewsTI_FastSoilCarbon.nc',              \
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_phenoLAI_MarthewsTI_climatology_FastSoilCarbon.nc',               \
                      WFDEI_dir+'qrparm.soil_merge_HWSD_NCSCD_cont_cosbyWFDEI.nc'                                          ]

SoilCarbon_plottags= [ 'JULES_100yrspin',                   \
                       'JULES-Trif_1000yrspin',             \
                       'JULES-Trif_100yrspin',              \
                       'SCsim_StatLAI',                     \
                       'SCsim_phenLAI_100yrspin',           \
                       'SCsim_phenLAI_10yrspin',            \
                       'SCsim_phenLAI_TRIFparams_10yrspin', \
                       'SCsim_1000yrspinLAI',               \
                       'SCsim_1000yrspinLAI_allTS',         \
                       'SCsim_ConstLAI',                    \
                       'SCsim_phenLAI_100yrspin_climat',    \
                       'HWSDNCSCD_SC'                       ]


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
# Read in Transcom Regions
#print 'TRANSCOM file: ',trans_file
#inf=nc.Dataset(trans_file,'r')
#transcom_regions_2D_temp=inf.variables['transcom_regions'][:]
#inf.close()
#transcom_regions_2D=np.zeros_like(transcom_regions_2D_temp)
#transcom_regions_2D[:,:360]=transcom_regions_2D_temp[:,360:]
#transcom_regions_2D[:,360:]=transcom_regions_2D_temp[:,:360]
#del transcom_regions_2D_temp


# Get ICE mask from sthu in JULES BS file
print 'Getting ICEmask from: '+SoilCarbon_files[0]
ICEmask=np.ones_like( \
        np.ma.masked_equal( \
            nc.Dataset(SoilCarbon_files[0],'r').variables['sthu'][0,0,:].squeeze() \
          , 0.0 ) )

# Open and read in Soil Carbon source data
print 'Soil Carbon file: '+SoilCarbon_files[iSource]
inf=nc.Dataset(SoilCarbon_files[iSource],'r')

if (SoilCarbon_Sources[iSource][:5]=='SCsim'):
    SC_data_1D = inf.variables['C_1SC_'+SoilRespFunctions[iSRF]][:].squeeze() \
                * ICEmask
    if ('limatology' in SoilCarbon_Sources[iSource]):
        SC_data_1D = np.mean(SC_data_1D,axis=1)
elif (SoilCarbon_Sources[iSource][:5]=='JULES'):
    SC_data_1D = np.mean(inf.variables['cs'][:].squeeze()*ICEmask,axis=0)
    if len(SC_data_1D.shape)==2:
        SC_data_1D = np.sum(SC_data_1D,axis=0)/2.
elif (SoilCarbon_Sources[iSource][:4]=='HWSD'):
    SC_data_1D = inf.variables['field1397'][:].squeeze()*ICEmask#/7.

SC_data_1D.shape

SC_data_2D=SC_data_1D[grindex]*grimask

        
print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'
#############################
print 'Plotting Soil Carbon '+region+' Map - '+SoilCarbon_Sources[iSource]
FIG = plt.figure(figsize=(18,12))
AXIS=FIG.add_subplot(1,1,1)
PTs.plot_map(SC_data_2D[:,:], \
             lons_2D,lats_2D, \
             AXIS=AXIS, \
             LAT_RANGE=lat_limits, LON_RANGE=lon_limits, \
             DATA_RANGE=CSlims, \
             LONDEL=gridspc[0], LATDEL=gridspc[1], \
             MPL_CBAR='YlOrBr', NLEVELS=nlevels, \
             CBAR_SIZE='4%',CBAR_PAD=0.4,\
             PLOT_TITLE='Soil Carbon - '+SoilCarbon_Sources[iSource], \
             CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
             RESOLUTION=resolution, \
             FONTSIZES=[12,12,14,18])

#FIG.suptitle('Soil Carbon - '+SoilCarbon_Sources[iSource], fontsize=30 )
print 'Writing image to: '+plots_DIR+region+'_CS_singlepool_map_'+plottag+'.png'
FIG.savefig(plots_DIR+region+'_CS_singlepool_map_'+plottag+'.png', bbox_inches='tight')
FIG.clear()
#plt.show()


