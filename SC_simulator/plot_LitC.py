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
import sys

##################################################
PLOTS=sys.argv[1]
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

LitC_maxvals   = [  1.0, 1.0, 1.0, 1.0, 1.0  ]

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
JULES_file        = JULES_output_dir+'WFDEI_pheno/JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno.monthly_mean.nc'
JULES_BSfile      = JULES_output_dir+'WFDEI_BigSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'
JULES_QSfile      = JULES_output_dir+'WFDEI_QuickSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_QuickSpin.monthly_mean.nc'
#

SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'
SCsim_staticLAI_file = 'SCsim_WL_JULES-WFDEI-Zinke-hydro1k_Static_LAI_FastSoilCarbon.nc'
SCsim_phenoLAI_file  = 'SCsim_J4.3_WFDEI_phenoLAI_MarthewsTI_FastSoilCarbon.nc'
SCsim_BigSpinLAI_file= 'SCsim_J4.3_WFDEI_BigSpinLAI_MarthewsTI_FastSoilCarbon.nc'

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

print 'TRANSCOM file: ',trans_file
inf=nc.Dataset(trans_file,'r')
transcom_regions_2D_temp=inf.variables['transcom_regions'][:]
inf.close()
transcom_regions_2D=np.zeros_like(transcom_regions_2D_temp)
transcom_regions_2D[:,:360]=transcom_regions_2D_temp[:,360:]
transcom_regions_2D[:,360:]=transcom_regions_2D_temp[:,:360]
del transcom_regions_2D_temp

# Open and read in simulated CS estimates
# 
# Static LAI simulation
#inf=nc.Dataset(SCsimulator_DIR+SCsim_staticLAI_file,'r')
#Loc_LitC_statLAI_1D  = inf.variables['Local_Litterfall_C'][:]
#Tot_LitC_statLAI_1D  = inf.variables['Total_Litterfall_C'][:]
#SoilResp_RsQ10_1D    = inf.variables['Soil_Resp_Fact_Q10t'][:]
#inf.close()

# Pheno LAI simulation
#inf=nc.Dataset(SCsimulator_DIR+SCsim_phenoLAI_file,'r')
#inf.close()

# BigSpin LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_BigSpinLAI_file,'r')
SC_BigSpinLAI_PFT_frac_1D      = inf.variables['PFT_frac'][:]
SC_BigSpinLAI_Loc_LitC_1D      = inf.variables['Local_Litterfall_C'][:] * \
                                    SC_BigSpinLAI_PFT_frac_1D
SC_BigSpinLAI_Tot_LitC_1D      = inf.variables['Total_Litterfall_C'][:]
SC_BigSpinLAI_SoilResp_RsQ10_1D= inf.variables['Soil_Resp_Fact_Q10t'][:]
SC_BigSpinLAI_C_4pools_Q10t_1D = inf.variables['C_4pools_Q10t'][:]
inf.close()

SC_BigSpinLAI_Loc_LitC_2D       = SC_BigSpinLAI_Loc_LitC_1D[:,:,grindex]*grimask
SC_BigSpinLAI_Tot_LitC_2D       = SC_BigSpinLAI_Tot_LitC_1D[:,grindex]*grimask
SC_BigSpinLAI_SoilResp_RsQ10_2D = SC_BigSpinLAI_SoilResp_RsQ10_1D[:,:,grindex]*grimask
SC_BigSpinLAI_C_4pools_Q10t_2D  = SC_BigSpinLAI_C_4pools_Q10t_1D[:,grindex]*grimask

SC_BigSpinLAI_Loc_LitC_2D_mean  = np.mean(SC_BigSpinLAI_Loc_LitC_2D,axis=1)
SC_BigSpinLAI_Tot_LitC_2D_mean  = np.mean(SC_BigSpinLAI_Tot_LitC_2D,axis=1)
SC_BigSpinLAI_SoilResp_RsQ10_2D_mean = np.mean(SC_BigSpinLAI_SoilResp_RsQ10_2D,axis=1)


# Open and read in the JULES cs data (and sthu for ICE/LAKE mask)
inf=nc.Dataset(JULES_BSfile)
JULES_BS_SoilCarbon=inf.variables['cs'][:].squeeze()
JULES_BS_LitC      =inf.variables['lit_c'][:].squeeze()
JULES_BS_SoilResp  =inf.variables['resp_s'][:].squeeze()
JULES_BS_frac      =inf.variables['frac'][:,:5,:,:].squeeze()
JULES_sthu         =inf.variables['sthu'][:].squeeze()
JULES_lats         =inf.variables['latitude'][:].squeeze()
JULES_lons         =inf.variables['longitude'][:].squeeze()
JULES_time         =nctime.num2date(inf.variables['time'][:], \
                           units=inf.variables['time'].units, \
                           calendar='standard'                )
inf.close()
JULES_sthu=np.ma.masked_equal(JULES_sthu,0.0)
ICEmask=np.ones_like(JULES_sthu[0,0,:])
JULES_BS_SoilCarbon=JULES_BS_SoilCarbon*ICEmask
JULES_BS_LitC      =JULES_BS_LitC*ICEmask*JULES_BS_frac
JULES_BS_SoilResp  =JULES_BS_SoilResp*ICEmask


####################################################################################
# Calculate JULES_BS means and Climatologies
JULES_BS_SoilCarbon_mean   = np.mean(JULES_BS_SoilCarbon,axis=0)
JULES_BS_LitC_mean         = np.mean(JULES_BS_LitC,axis=0)
JULES_BS_SoilResp_mean     = np.mean(JULES_BS_SoilResp,axis=0)

JULES_BS_SoilCarbon_climat = \
    np.mean(JULES_BS_SoilCarbon.reshape([nYEARs,nMONTHs,nPOOLs,nlandpoints]),axis=0)
JULES_BS_LitC_climat = \
    np.mean(JULES_BS_LitC.reshape([nYEARs,nMONTHs,nPFTs,nlandpoints]),axis=0)
JULES_BS_SoilResp_climat = \
    np.mean(JULES_BS_SoilResp.reshape([nYEARs,nMONTHs,nPOOLs,nlandpoints]),axis=0)

JULES_BS_SoilCarbon_mean_2D = JULES_BS_SoilCarbon_mean[:,grindex]*grimask
JULES_BS_LitC_mean_2D       = JULES_BS_LitC_mean[:,grindex]*grimask
JULES_BS_SoilResp_mean_2D   = JULES_BS_SoilResp_mean[:,grindex]*grimask

JULES_BS_SoilCarbon_climat_2D = JULES_BS_SoilCarbon_climat[:,:,grindex]*grimask
JULES_BS_LitC_climat_2D       = JULES_BS_LitC_climat[:,:,grindex]*grimask
JULES_BS_SoilResp_climat_2D   = JULES_BS_SoilResp_climat[:,:,grindex]*grimask

########################################################################

# create TRANSCOM indices
flindex = ( ((JULES_lats+89.75)*2.).astype('int'), \
            ((JULES_lons+179.75)*2.).astype('int') )

transcom_regions = np.zeros_like(JULES_lats)
transcom_regions = transcom_regions_2D[flindex]

TRANSCOM_names = ['Global',            \
                  'N. Amer Boreal',    \
                  'N. Amer Temperate', \
                  'S. Amer Tropical',  \
                  'S. Amer Temperate', \
                  'N. Africa',         \
                  'S. Africa',         \
                  'Asia Boreal',       \
                  'Asia Temperate',    \
                  'Asia Tropical',     \
                  'Australia',         \
                  'Europe'             ]
TRANSCOM_indices = [ np.where((transcom_regions>=1)&(transcom_regions<=11))[0], \
                     np.where(transcom_regions==1)[0], \
                     np.where(transcom_regions==2)[0], \
                     np.where(transcom_regions==3)[0], \
                     np.where(transcom_regions==4)[0], \
                     np.where(transcom_regions==5)[0], \
                     np.where(transcom_regions==6)[0], \
                     np.where(transcom_regions==7)[0], \
                     np.where(transcom_regions==8)[0], \
                     np.where(transcom_regions==9)[0], \
                     np.where(transcom_regions==10)[0],\
                     np.where(transcom_regions==11)[0] ]

###################################################################################


# Plot regional Climatology time series of JULES_BS
print 'Plotting Regional Climatology Time-series'
if (PLOTS[0]=='Y'):
    #JULES_BS_SoilCarbon_regional_TS=[]
    FIG_J=plt.figure(figsize=(18,18))
    FIG_S=plt.figure(figsize=(18,18))
    
    for iREGION in select_trans:
        AXIS_J=FIG_J.add_subplot(4,3,iREGION+1)
        AXIS_S=FIG_S.add_subplot(4,3,iREGION+1)

        REGIONname=TRANSCOM_names[iREGION]
        index=TRANSCOM_indices[iREGION]
        print REGIONname, index[0].shape
        JULES_temp_TS   = np.mean(JULES_BS_LitC_climat[:,:,index],axis=2)
        SC_temp_TS      = np.mean(SC_BigSpinLAI_Loc_LitC_1D[:,:,index],axis=2)
        
        for iPFT in range(nPFTs):
            AXIS_J.plot(np.arange(12)+1,JULES_temp_TS[:,iPFT],\
                        label=PFT_longnames[iPFT],lw=2)
            AXIS_J.set_ylabel('Local Litterfall Carbon ($kg$C $m^-2$)',fontsize=16)
            AXIS_J.set_title(REGIONname,fontsize=20)
            
            #print SC_temp_TS[iPFT,:]
            AXIS_S.plot(np.arange(12)+1,SC_temp_TS[iPFT,:],\
                        label=PFT_longnames[iPFT],lw=2)
            AXIS_S.set_ylabel('Local Litterfall Carbon ($kg$C $m^-2$)',fontsize=16)
            AXIS_S.set_title(REGIONname,fontsize=20)
    
    FONT = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 20 }
    
    FIG_J.suptitle('JULES Big Spin Litterfall',fontsize=30)
    HANDLES,LABELS=AXIS_J.get_legend_handles_labels()
    FIG_J.legend(HANDLES,LABELS,loc='lower center', \
                    ncol=len(select_trans), prop=FONT)
    FIG_J.savefig(plots_DIR+'JULES_WFDEI_BigSpin_Reg_LitC_TS.png')
    FIG_J.clear()

    FIG_S.suptitle('Soil Carbon Simulator Litterfall',fontsize=30)
    HANDLES,LABELS=AXIS_S.get_legend_handles_labels()
    FIG_S.legend(HANDLES,LABELS,loc='lower center', \
                    ncol=len(select_trans), prop=FONT)
    FIG_S.savefig(plots_DIR+'SCsim_BigSpin_Reg_LitC_TS.png')
    FIG_S.clear()



###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
#CLEVELS=pool_CLEVELS[ipool], \
if (PLOTS[1]=='Y'):
    for imonth in range(nMONTHs):
        print 'Plotting LitC Climatology Map for Month: '+str(imonth)
        #FIG,AXIS_set=plt.subplots(nrows=2,ncols=2,figsize=(18,12))
        FIG = plt.figure(figsize=(18,12))
        AXIS = FIG.add_subplot(3,2,i)
        PTs.plot_map(JULES_BS_SoilCarbon_climat_2D[imonth,ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')


        for ipft in range(nPFTs):
            AXIS=FIG.add_subplot(3,2,ipft+2)
            PTs.plot_map(JULES_BS_SoilCarbon_climat_2D[imonth,ipft,:,:], \
                         lons_2D,lats_2D, \
                         AXIS=AXIS, \
                         DATA_RANGE=[0,pool_maxvals[ipool]], \
                         MPL_CBAR='YlOrBr', NLEVELS=10, \
                         CBAR_SIZE='8%',CBAR_PAD=0.2,\
                         PLOT_TITLE=pool_longnames[ipool], \
                         CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                         RESOLUTION='c')

        FIG.suptitle('JULES Soil Carbon from 1000 year TRIFFID spin-up run, '+str(imonth+1) )
        FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_CLIMATOLOGY_'+str(imonth))
        FIG.clear()    #


#############################
if (PLOTS[2]=='Y'):
    print 'Plotting JULES_BS Mean Litterfall Global Maps'
    plot_tag='JULES_BigSpin'
    FIG = plt.figure(figsize=(18,12))
    ipool=3    # Plot Hummus only
    AXIS = FIG.add_subplot(3,2,1)
    PTs.plot_map(JULES_BS_SoilCarbon_mean_2D[ipool,:,:], \
                 lons_2D,lats_2D, \
                 AXIS=AXIS, \
                 DATA_RANGE=[0,pool_maxvals[ipool]], \
                 MPL_CBAR='YlOrBr', NLEVELS=10, \
                 CBAR_SIZE='8%',CBAR_PAD=0.2,\
                 PLOT_TITLE='Soil Carbon as '+pool_longnames[ipool], \
                 CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                 RESOLUTION='c')

    for ipft in range(nPFTs):
        AXIS=FIG.add_subplot(3,2,ipft+2)
        PTs.plot_map(JULES_BS_LitC_mean_2D[ipft,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE='Litterfall from '+PFT_longnames[ipft], \
                     CBAR_LABEL='Litterfall ($kg$C $m^-2$)', \
                     RESOLUTION='c')
    
    FIG.suptitle('JULES 1000 year TRIFFID spin-up run',fontsize=24)
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[3]=='Y'):
    print 'Plotting SCsim_BS Litterfall Global Maps'
    plot_tag='SCsim_BigSpin'
    FIG = plt.figure(figsize=(18,12))
    ipool=3    # Plot Hummus only
    AXIS = FIG.add_subplot(3,2,1)
    PTs.plot_map(SC_BigSpinLAI_C_4pools_Q10t_2D[ipool,:,:], \
                 lons_2D,lats_2D, \
                 AXIS=AXIS, \
                 DATA_RANGE=[0,pool_maxvals[ipool]], \
                 MPL_CBAR='YlOrBr', NLEVELS=10, \
                 CBAR_SIZE='8%',CBAR_PAD=0.2,\
                 PLOT_TITLE='Soil Carbon as '+pool_longnames[ipool], \
                 CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                 RESOLUTION='c')

    for ipft in range(nPFTs):
        AXIS=FIG.add_subplot(3,2,ipft+2)
        PTs.plot_map(SC_BigSpinLAI_Loc_LitC_2D_mean[ipft,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=PFT_longnames[ipft], \
                     CBAR_LABEL='Litterfall ($kg$C $m^-2$)', \
                     RESOLUTION='c')
    
    FIG.suptitle('Simulated Litterfall using JULES LAI and temperature data',fontsize=24)
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')


quit()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[4]=='Y'):
    print 'Plotting SCsim_pheno Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_PhenoLAI_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_phenoLAI_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')
    
    FIG.suptitle('Simulated Soil Carbon using JULES 10 year spin-up data')
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()



###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[5]=='Y'):
    print 'Plotting SCsim_BigSpin Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_BigSpinLAI_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_BigSpinLAI_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')

    FIG.suptitle('Simulated Soil Carbon using JULES-Triffid 1000 year spin-up data')
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()





