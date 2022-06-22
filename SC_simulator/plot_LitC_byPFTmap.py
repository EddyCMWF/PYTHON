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
import sys

##################################################
#PLOTS=sys.argv[1]
#PLOTS[1] = JULES_BS, timeseries
#PLOTS[2] = JULES_BS, Climatology GLOBAL map
#PLOTS[3] = JULES_BS, mean global map
#PLOTS[4] = SCsim_staticLAI, mean global map
#PLOTS[5] = SCsim_phenoLAI, mean global map

##########################################
# Set parameters and variables
nYEARs=33
nMONTHs=12
nlandpoints=67209
nPOOLs=4
nPFTs=5

PFT_names      = [ 'BL','NL','C3','C4','Sh']
PFT_longnames  = [ 'Broadleaf','Needleleaf','C3 Grass','C4 Grass','Shrub' ]
PFT_colours    = [ 'lime','darkgreen', 'yellow','yellowgreen','orange' ]

LitC_maxvals   = [  0.5, 0.5, 0.5, 0.5, 0.5  ]

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



select_trans = range(12)        



########################################################################
print 'Setting Filenames and Directories'
# Input files:
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/'
SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'


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
                       'SCsim pheno LAI (100 year spin) Climatology'      ]

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
                      SCsimulator_DIR+'SCsim_J4.3_WFDEI_phenoLAI_MarthewsTI_climatology_FastSoilCarbon.nc'                ]

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
                       'SCsim_phenLAI_100yrspin_climat'     ]



print 'Available Soil Carbon Sources: '
for i in range(len(SoilCarbon_Sources)):
    print str(i)+': '+SoilCarbon_Sources[i]

iSource = input('Select a source: ')
plottag = SoilCarbon_plottags[iSource]

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

if (SoilCarbon_Sources[iSource][:5]=='SCsim'):
    SoilRespFunctions = ['Q10t','RothCit']
    print 'Available Soil Respiration Functions: '
    for i in range(len(SoilRespFunctions)):
        print str(i)+': '+SoilRespFunctions[i]

    iSRF=input('Select an Soil Respiration Funtion: ')

    plottag=plottag+SoilRespFunctions[iSRF]

plots_DIR          = '/users/eow/edwcom/SC_simulator/output/plots/WFDEI/GLOBAL/LitC/'

trans_file   = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5_orig.nc'
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
    SC_data_1D   = inf.variables['C_4pools_'+SoilRespFunctions[iSRF]][:].squeeze() \
                * ICEmask
    LitC_data_1D = inf.variables['Litterfall_C'][:].squeeze() \
                 * ICEmask

    if ('limatology' in SoilCarbon_Sources[iSource]):
        SC_data_1D = np.mean(SC_data_1D,axis=1)
    
    if not ('tatic' in SoilCarbon_Sources[iSource]):
        LitC_data_1D = np.mean(LitC_data_1D,axis=1)

elif (SoilCarbon_Sources[iSource][:5]=='JULES'):
    SC_data_1D = np.mean(inf.variables['cs'][:].squeeze()*ICEmask,axis=0)
    LitC_data_1D = np.mean(  inf.variables['lit_c'][:].squeeze() \
                           * inf.variables['frac'][:,:5,:,:].squeeze() \
                           * ICEmask,  axis=0 )

inf.close()

SC_data_2D=SC_data_1D[:,grindex]*grimask
LitC_data_2D=LitC_data_1D[:,grindex]*grimask

print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

#################################################################################
print 'Plotting Litterfall C Global Maps - '+SoilCarbon_Sources[iSource]
FIG = plt.figure(figsize=(18,12))
ipool=3    # Plot Hummus only
AXIS = FIG.add_subplot(3,2,1)
PTs.plot_map(SC_data_2D[ipool,:,:], \
             lons_2D,lats_2D, \
             AXIS=AXIS, \
             DATA_RANGE=[0,pool_maxvals[ipool]], \
             MPL_CBAR='YlOrBr', NLEVELS=10, \
             CBAR_SIZE='7%',CBAR_PAD=0.1,\
             PLOT_TITLE='Soil Carbon as '+pool_longnames[ipool], \
             CBAR_LABEL='Soil Carbon ($kg$C $m^{-2}$)', \
             RESOLUTION='c')

for ipft in range(nPFTs):
    AXIS=FIG.add_subplot(3,2,ipft+2)
    PTs.plot_map(LitC_data_2D[ipft,:,:], \
                 lons_2D,lats_2D, \
                 AXIS=AXIS, \
                 DATA_RANGE=[0,LitC_maxvals[ipft]], \
                 MPL_CBAR='YlOrBr', NLEVELS=10, \
                 CBAR_SIZE='7%',CBAR_PAD=0.1,\
                 PLOT_TITLE='Litterfall from '+PFT_longnames[ipft], \
                 CBAR_LABEL='Litterfall ($kg$C $m^{-2}$ $yr^{-1}$)', \
                 RESOLUTION='c')
    CORRELATION = stats.pearsonr(SC_data_2D[ipool,:,:].flatten(),\
                                 LitC_data_2D[ipft,:,:].flatten())
    AXIS.text(-175,-60,"Pearson's R = "+str(round(CORRELATION[0],2)))

FIG.subplots_adjust(wspace=-0.2,hspace=0.3)
FIG.suptitle('Litterfall Carbon by PFT - '+SoilCarbon_Sources[iSource],fontsize=24)
print 'Writing image to: '+plots_DIR+'Global_LitC_map_'+plottag+'.png'
FIG.savefig(plots_DIR+'Global_LitC_map_'+plottag+'.png', bbox_inches='tight')
FIG.clear()






