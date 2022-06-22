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
REG_CSlims   = [ (0,15.), None, (0,20) ]
REG_LitClims = [ (0,1.), None, (0,1.5) ]
REG_Rslims   = [ (0,15), None, (0,15) ]


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
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'
JULES_output_dir  = '/users/eow/edwcom/SC_simulator/JULES_output/'
SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'

plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/WFDEI/'+region+'/CS/'
plots_DIR = plots_DIR.replace('-fixed','')
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
print 'Grid file: ',WFDEI_gridfile
inf=nc.Dataset(WFDEI_gridfile,'r')
lats   =inf.variables['latitude'][:]
lons   =inf.variables['longitude'][:]
grindex=inf.variables['land_index'][:]-1   # minus 1 to convert from Fortran to C index
grifrac=inf.variables['land_fraction'][:]
grimask=np.ones_like(grindex)
inf.close()
lons_2D,lats_2D=np.meshgrid(lons,lats)
plot_lons = lons_2D[ lat_indices[0]:lat_indices[1], \
                     lon_indices[0]:lon_indices[1]  ]
plot_lats = lats_2D[ lat_indices[0]:lat_indices[1], \
                     lon_indices[0]:lon_indices[1]  ]

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
    LitC_data_1D = inf.variables['Total_Litterfall_C'][:].squeeze() \
                 * ICEmask
    Rs_data_1D   = inf.variables['Soil_Resp_Fact_'+SoilRespFunctions[iSRF]][:].squeeze() \
                 * ICEmask

    if ('limatology' in SoilCarbon_Sources[iSource]):
        SC_data_1D = np.mean(SC_data_1D,axis=1)
    
    if not ('tatic' in SoilCarbon_Sources[iSource]):
        LitC_data_1D = np.mean(LitC_data_1D,axis=0)

    Rs_data_1D=np.mean(Rs_data_1D,axis=1)*3600.*24.*365.

elif (SoilCarbon_Sources[iSource][:5]=='JULES'):
    SC_data_1D = np.mean(inf.variables['cs'][:].squeeze()*ICEmask,axis=0)
    LitC_data_1D = np.sum( np.mean(  inf.variables['lit_c'][:].squeeze() \
                                   * inf.variables['frac'][:,:5,:,:].squeeze() \
                                   * ICEmask,  axis=0 )  , axis=0 )
    Rs_data_1D = np.mean(inf.variables['resp_s'][:].squeeze()*ICEmask, axis=0) \
                * 3600.*24.*365./SC_data_1D


inf.close()

SC_data_2D       = SC_data_1D[:,grindex]*grimask
total_SC_data_2D = np.sum(SC_data_2D,axis=0)
plot_total_SC    = total_SC_data_2D[ lat_indices[0]:lat_indices[1], \
                                     lon_indices[0]:lon_indices[1]  ]
LitC_data_2D     = LitC_data_1D[grindex]*grimask
plot_LitC        = LitC_data_2D[ lat_indices[0]:lat_indices[1], \
                                 lon_indices[0]:lon_indices[1]  ]
Rs_data_2D       = Rs_data_1D[:,grindex]*grimask
total_Rs_data_2D = np.sum(Rs_data_2D,axis=0)
plot_total_Rs    = total_Rs_data_2D[ lat_indices[0]:lat_indices[1], \
                                     lon_indices[0]:lon_indices[1]  ]

if (CSlims==None):   CSlims  =[0,np.ceil(np.mean(plot_total_SC)*2.0)]
if (LitClims==None): LitClims=[0,np.ceil(np.mean(plot_LitC)*2.0)]
if (Rslims==None):   Rslims  =[0,np.ceil(np.mean(plot_total_Rs)*2.0)]


print '8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o 8=o'

#################################################################################
print 'Plotting Litterfall C '+region+' maps - '+SoilCarbon_Sources[iSource]
FIG = plt.figure(figsize=(18,12))

AXIS = FIG.add_subplot(3,1,1)
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

AXIS=FIG.add_subplot(3,1,2)
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
#CORRELATION = stats.pearsonr(total_SC_data_2D.flatten(),\
#                             LitC_data_2D.flatten())
#AXIS.text(textpos[0],textpos[1],"Pearson's R = "+str(round(CORRELATION[0],2)), fontsize=10)


AXIS=FIG.add_subplot(3,1,3)
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

#CORRELATION = stats.pearsonr(total_SC_data_2D.flatten(),\
#                             total_Rs_data_2D.flatten())
#AXIS.text(textpos[0],textpos[1],\
#          "Pearson's R = "+str(round(CORRELATION[0],2)),fontsize=10)
#CORRELATION = stats.pearsonr(LitC_data_2D.flatten(),\
#                             total_Rs_data_2D.flatten())
#AXIS.text(textpos[0],textpos[1]+textpos[2],\
#          "Pearson's R = "+str(round(CORRELATION[0],2)),fontsize=10)


FIG.subplots_adjust(hspace=0.3)
FIG.suptitle('Soil Carbon Breakdown - \n'+SoilCarbon_Sources[iSource],fontsize=24)
print 'Writing image to: '+plots_DIR+region+'_LitC_map_'+plottag+'.png'
FIG.savefig(plots_DIR+region+'_CS_breakdown_map_'+plottag+'.png', bbox_inches='tight')
FIG.clear()






