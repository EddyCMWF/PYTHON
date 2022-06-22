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
JULES_file        = JULES_output_dir+'WFDEI_pheno/JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno.monthly_mean.nc'
JULES_BSfile      = JULES_output_dir+'WFDEI_BigSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'
JULES_QSfile      = JULES_output_dir+'WFDEI_QuickSpin/JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_QuickSpin.monthly_mean.nc'
#

SCsimulator_DIR      = '/users/eow/edwcom/SC_simulator/output/'
SCsim_staticLAI_file = 'SCsim_WL_JULES-WFDEI-Zinke-hydro1k_Static_LAI_FastSoilCarbon.nc'
SCsim_phenoLAI_file  = 'SCsim_J4.3_WFDEI_phenoLAI_MarthewsTI_FastSoilCarbon.nc'
SCsim_QSphenoLAI_file= 'SCsim_J4.3_WFDEI_QSphenoLAI_MarthewsTI_FastSoilCarbon.nc'
SCsim_QSphenoLAI_trifparam_file= 'SCsim_J4.3_WFDEI_QSphenoLAI_trifparam_MarthewsTI_FastSoilCarbon.nc'
SCsim_BigSpinLAI_file= 'SCsim_J4.3_WFDEI_BigSpinLAI_MarthewsTI_FastSoilCarbon.nc'
SCsim_BSLAI_allTS_file='SCsim_J4.3_WFDEI_BigSpinLAI_alltimesteps_MarthewsTI_FastSoilCarbon.nc'
SCsim_constLAI_file   ='SCsim_J4.3_WFDEI_ConstLAI_alltimesteps_MarthewsTI_FastSoilCarbon.nc'

plots_DIR          = '/users/eow/edwcom/SC_simulator/output/plots/WFDEI/GLOBAL/CS/'

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
inf=nc.Dataset(SCsimulator_DIR+SCsim_staticLAI_file,'r')
SC_statLAI_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()

SC_statLAI_C_4pools_Q10t_2D  =SC_statLAI_C_4pools_Q10t_1D[:,grindex]*grimask
# 

# Pheno LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_phenoLAI_file,'r')
SC_phenoLAI_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()

SC_phenoLAI_C_4pools_Q10t_2D  =SC_phenoLAI_C_4pools_Q10t_1D[:,grindex]*grimask

#
# Pheno LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_QSphenoLAI_file,'r')
SC_QSphenoLAI_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()
SC_QSphenoLAI_C_4pools_Q10t_2D  =SC_QSphenoLAI_C_4pools_Q10t_1D[:,grindex]*grimask

#
# Pheno LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_QSphenoLAI_trifparam_file,'r')
SC_QSphenoLAI_trifparam_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()
SC_QSphenoLAI_trifparam_C_4pools_Q10t_2D  =SC_QSphenoLAI_trifparam_C_4pools_Q10t_1D[:,grindex]*grimask

#
# BigSpin LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_BigSpinLAI_file,'r')
SC_BigSpinLAI_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()
SC_BigSpinLAI_C_4pools_Q10t_2D  =SC_BigSpinLAI_C_4pools_Q10t_1D[:,grindex]*grimask


# BigSpin LAI allTS simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_BSLAI_allTS_file,'r')
SC_BSLAI_allTS_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()
SC_BSLAI_allTS_C_4pools_Q10t_2D  = SC_BSLAI_allTS_C_4pools_Q10t_1D[:,grindex]*grimask

# Constant LAI simulation
inf=nc.Dataset(SCsimulator_DIR+SCsim_constLAI_file,'r')
SC_ConstLAI_C_4pools_Q10t_1D  = inf.variables['C_4pools_Q10t'][:]
inf.close()
SC_ConstLAI_C_4pools_Q10t_2D  = SC_ConstLAI_C_4pools_Q10t_1D[:,grindex]*grimask


# Open and read in the JULES cs data (and sthu for ICE/LAKE mask)
inf=nc.Dataset(JULES_BSfile)
JULES_BS_SoilCarbon=inf.variables['cs'][:].squeeze()
JULES_sthu=inf.variables['sthu'][:].squeeze()
JULES_lats=inf.variables['latitude'][:].squeeze()
JULES_lons=inf.variables['longitude'][:].squeeze()
JULES_time=nctime.num2date(inf.variables['time'][:],          \
                           units=inf.variables['time'].units, \
                           calendar='standard'                )
inf.close()
JULES_sthu=np.ma.masked_equal(JULES_sthu,0.0)
ICEmask=np.ones_like(JULES_sthu)
JULES_BS_SoilCarbon=JULES_BS_SoilCarbon*ICEmask

########################################################################

# create TRANSCOM indices
flindex = ( ((JULES_lats+89.75)*2.).astype('int'), \
            ((JULES_lons+179.75)*2.).astype('int') )

transcom_regions = np.zeros_like(JULES_lats)
transcom_regions = transcom_regions_2D[flindex]

TRANSCOM_names = ['N. Amer Boreal',    \
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
TRANSCOM_indices = [ np.where(transcom_regions==1)[0], \
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


####################################################################################
# Calculate JULES_BS means and Climatologies
JULES_BS_SoilCarbon_mean   = np.mean(JULES_BS_SoilCarbon,axis=0)
JULES_BS_SoilCarbon_climat = \
    np.mean(JULES_BS_SoilCarbon.reshape([nYEARs,nMONTHs,nPOOLs,nlandpoints]),axis=0)

JULES_BS_SoilCarbon_mean_2D   = JULES_BS_SoilCarbon_mean[:,grindex]*grimask
JULES_BS_SoilCarbon_climat_2D = JULES_BS_SoilCarbon_climat[:,:,grindex]*grimask


JULES_BS_SoilCarbon_Deseasoned = \
       ( JULES_BS_SoilCarbon.reshape([nYEARs,nMONTHs,nPOOLs,nlandpoints])  \
       - JULES_BS_SoilCarbon_climat).reshape([nYEARs*nMONTHs,nPOOLs,nlandpoints]) 


###################################################################################
# Plot regional time series of JULES_BS
plot_tag='JULES_WFDEI_BigSpin'
if (PLOTS[0]=='Y'):
    #JULES_BS_SoilCarbon_regional_TS=[]
    fig_Tots=plt.figure(figsize=(18,12))
    fig_PAns=plt.figure(figsize=(18,12))
    fig_Ans=plt.figure(figsize=(18,12))
    fig_DeAns=plt.figure(figsize=(18,12))
    fig_DePAns=plt.figure(figsize=(18,12))
    
    for iREGION in select_trans:
        name=TRANSCOM_names[iREGION]
        index=TRANSCOM_indices[iREGION]
        print name, index[0].shape
        temp_TS   = np.mean(JULES_BS_SoilCarbon[:,:,index],axis=2)
        temp_An   = temp_TS-np.mean(temp_TS,axis=0)
        temp_PAn  = ((temp_TS*100.)/np.mean(temp_TS,axis=0))-100.
        temp_DeTS = np.mean(JULES_BS_SoilCarbon_Deseasoned[:,:,index],axis=2)
        temp_DeAn = temp_DeTS-np.mean(temp_DeTS,axis=0)
        temp_DePAn= ((temp_DeTS*100.)/np.mean(temp_TS,axis=0))
        #JULES_BS_SoilCarbon_regional_TS.append(temp_TS)
        for iPOOL in range(nPOOLs):
            
            # Plot Totals
            ax_tot=fig_Tots.add_subplot(2,2,iPOOL+1)
            ax_tot.plot(JULES_time,temp_TS[:,iPOOL],label=name,lw=2)
            ax_tot.set_ylabel('Soil Carbon ($kg$C $m^-2$)',fontsize=16)
            ax_tot.set_title(pool_names[iPOOL],fontsize=20)
            
            # Plot  anomaly
            ax_an=fig_Ans.add_subplot(2,2,iPOOL+1)
            ax_an.plot(JULES_time,temp_An[:,iPOOL],label=name,lw=2)
            ax_an.set_ylabel('Soil Carbon Anomaly ($kg$C $m^-2$)',fontsize=16)
            ax_an.set_title(pool_names[iPOOL],fontsize=20)
            
            # Plot percentage anomaly
            ax_pan=fig_PAns.add_subplot(2,2,iPOOL+1)
            ax_pan.plot(JULES_time,temp_PAn[:,iPOOL],label=name,lw=2)
            ax_pan.set_ylabel('Soil Carbon Percentage Anomaly (\%)',fontsize=16)
            ax_pan.set_title(pool_names[iPOOL],fontsize=20)
            
            # Plot deseasoned anomaly
            ax_dean=fig_DeAns.add_subplot(2,2,iPOOL+1)
            ax_dean.plot(JULES_time,temp_DeAn[:,iPOOL],label=name,lw=2)
            ax_dean.set_ylabel('Deseasoned Soil Carbon Anom. ($kg$C $m^-2$)',fontsize=16)
            ax_dean.set_title(pool_names[iPOOL],fontsize=20)
            
            # Plot deseasoned percentage anomaly
            ax_depan=fig_DePAns.add_subplot(2,2,iPOOL+1)
            ax_depan.plot(JULES_time,temp_DePAn[:,iPOOL],label=name,lw=2)
            ax_depan.set_ylabel('Deseasoned Soil Carbon Perc. Anom. (\%)',fontsize=16)
            ax_depan.set_title(pool_names[iPOOL],fontsize=20)

        
    FONT = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 20 }
    
    HANDLES,LABELS=ax_tot.get_legend_handles_labels()
    fig_Tots.legend(HANDLES,LABELS,loc='lower center', \
                    ncol=len(select_trans), prop=FONT)
    fig_Tots.savefig(plots_DIR+plot_tag+'_Reg_cs_TS.png')
    fig_Tots.clear()

    HANDLES,LABELS=ax_an.get_legend_handles_labels()
    fig_Ans.legend(HANDLES,LABELS,loc='lower center', \
                   ncol=len(select_trans), prop=FONT)
    fig_Ans.savefig(plots_DIR+plot_tag+'_Reg_cs_TS_anom.png')
    fig_Ans.clear()

    HANDLES,LABELS=ax_pan.get_legend_handles_labels()
    fig_PAns.legend(HANDLES,LABELS,loc='lower center', \
                    ncol=len(select_trans), prop=FONT)
    fig_PAns.savefig(plots_DIR+plot_tag+'_Reg_cs_TS_perc_anom.png')
    fig_PAns.clear()
    
    HANDLES,LABELS=ax_dean.get_legend_handles_labels()
    fig_DeAns.legend(HANDLES,LABELS,loc='lower center', \
                     ncol=len(select_trans), prop=FONT)
    fig_DeAns.savefig(plots_DIR+plot_tag+'_Reg_cs_TS_desea_anom.png')
    fig_DeAns.clear()

    HANDLES,LABELS=ax_depan.get_legend_handles_labels()
    fig_DePAns.legend(HANDLES,LABELS,loc='lower center', \
                      ncol=len(select_trans), prop=FONT)
    fig_DePAns.savefig(plots_DIR+plot_tag+'_Reg_cs_TS_desea_perc_anom.png')
    fig_DePAns.clear()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
#CLEVELS=pool_CLEVELS[ipool], \
if (PLOTS[1]=='Y'):
    for imonth in range(nMONTHs):
        print 'Plotting Climatology Map for Month: '+str(imonth)
        #FIG,AXIS_set=plt.subplots(nrows=2,ncols=2,figsize=(18,12))
        FIG = plt.figure(figsize=(18,12))
        for ipool in range(nPOOLs):
            AXIS=FIG.add_subplot(2,2,ipool+1)
            PTs.plot_map(JULES_BS_SoilCarbon_climat_2D[imonth,ipool,:,:], \
                         lons_2D,lats_2D, \
                         AXIS=AXIS, \
                         DATA_RANGE=[0,pool_maxvals[ipool]], \
                         MPL_CBAR='YlOrBr', NLEVELS=10, \
                         CBAR_SIZE='8%',CBAR_PAD=0.2,\
                         PLOT_TITLE=pool_longnames[ipool], \
                         CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                         RESOLUTION='c')

        FIG.suptitle('JULES Soil Carbon from 1000 year TRIFFID spin-up run, '+str(imonth+1), fontsize=30 )
        FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_CLIMATOLOGY_'+str(imonth))
        FIG.clear()    #


#############################
if (PLOTS[2]=='Y'):
    print 'Plotting JULES_BS Mean Soil Carbon Global Map'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(JULES_BS_SoilCarbon_mean_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')
    
    FIG.suptitle('JULES Soil Carbon from 1000 year TRIFFID spin-up run', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[3]=='Y'):
    print 'Plotting SCsim_static Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_StaticLAI_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_statLAI_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')

    FIG.suptitle('Simulated Soil Carbon using JULES 100 year spin-up data, static LAI map', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()


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
    
    FIG.suptitle('Simulated Soil Carbon using JULES 100 year pheno spin-up data', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[5]=='Y'):
    print 'Plotting SCsim_QSpheno Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_QSphenoLAI_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_QSphenoLAI_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')
    
    FIG.suptitle('Simulated Soil Carbon using JULES 10 year pheno spin-up data', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()



###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[6]=='Y'):
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

    FIG.suptitle('Simulated Soil Carbon using JULES-Triffid 1000 year spin-up data', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()



###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[7]=='Y'):
    print 'Plotting SCsim_BigSpin Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_BSLAI_allTS_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_BSLAI_allTS_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')

    FIG.suptitle('Simulated Soil Carbon using JULES-Triffid 1000 year spin-up data (all TS)', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP_mean')
    FIG.clear()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[8]=='Y'):
    print 'Plotting SCsim_BigSpin Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_ConstLAI_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_ConstLAI_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')

    FIG.suptitle('Simulated Soil Carbon using Constant LAI values [9,5,4,4,3]', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP')
    FIG.clear()


###################################################################################
# Plot GLOBAL Maps of JULES_BS (Climatology and Mean)
if (PLOTS[9]=='Y'):
    print 'Plotting SCsim_QSphenoLAI_trifparams Mean Soil Carbon Global Map'
    plot_tag = 'SCsim_QSpheno_trifparam_RsQ10'
    FIG = plt.figure(figsize=(18,12))
    for ipool in range(nPOOLs):
        AXIS=FIG.add_subplot(2,2,ipool+1)
        PTs.plot_map(SC_QSphenoLAI_trifparam_C_4pools_Q10t_2D[ipool,:,:], \
                     lons_2D,lats_2D, \
                     AXIS=AXIS, \
                     DATA_RANGE=[0,pool_maxvals[ipool]], \
                     MPL_CBAR='YlOrBr', NLEVELS=10, \
                     CBAR_SIZE='8%',CBAR_PAD=0.2,\
                     PLOT_TITLE=pool_longnames[ipool], \
                     CBAR_LABEL='Soil Carbon ($kg$C $m^-2$)', \
                     RESOLUTION='c')

    FIG.suptitle('Simulated Soil Carbon using JULES pheno run, trif LAI params, 10 year spin', fontsize=30 )
    FIG.savefig(plots_DIR+plot_tag+'_GLOBAL_MAP')
    FIG.clear()


