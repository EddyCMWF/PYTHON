#!/usr/bin/env python

################################################################################
# 
# Program: SC_simulator.py
# 
# Python Script to simulate soil carbon estimate for a long-term equilibrium
#   condition, i.e. dCs/dt=0 .'. litterfall=microbial_soil_respiraiton       
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
################################################################################
#
# 
#  
################################################################################
#
#
import pylab as plt
import numpy as np
import netCDF4 as nc
import JULES_SC_functions as J_SC_F
import data_info_SC as di_SC
import plot_tools as PT
#
#####################################################################################
# 0. Set flags and diags
###########################
plot_SMdiags=False
plot_VEGdiags=False
plot_LITdiags=False
plot_LAIdiags=False
plot_SOILdiags=False
plot_NFLITdiags=False
plot_Cq10diags=False
plot_Crothdiags=False
nPOOLs=4
fill_value=-999
#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:

JULES_sources, JULES_sources_info = di_SC.jules_sources_info()

iJULES = di_SC.select_source(JULES_sources,Message='Select JULES source:')
JULES_source=JULES_sources[iJULES]
JULES_source_info=JULES_sources_info[JULES_source]
print 'JULES file = ',JULES_source_info['file']

LAI_source_info = di_SC.LAI_source_info(JULES_source_info)

PFT_source_info = di_SC.PFT_source_info(JULES_source_info)

SOIL_source_info = di_SC.SOIL_source_info(JULES_source_info)

#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
print 'Current outdir: '+OUT_DIR
temp = raw_input('Enter alternative OUT_DIR or hit return to continue: \n')

if temp!='':
    OUT_DIR=temp.copy()

out_tag = 'SCsim_'+JULES_source_info['tag'] + \
          '_'+LAI_source_info['tag']        + \
          '_'+PFT_source_info['tag']

#quit()

#

#####################################################################################
# 2. Define parameters:
##################################
#
print 'Defining parameters'
# 2.0 Define Variables:
# 2.0.1 Soil Specific Respiration Rate for RothC pools and single soil pool
kappa_s_RothC  = np.array([ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10 ])
kappa_s_1SC    = [0.5e-8]
# 2.0.2 Soil_Resp_Fractor beta_r in documetation, formulation from JULES src
clay_content=0.23   # as in JULES
beta_r  = 1.0 / (4.0895 + 2.672*  np.exp(-0.0786 * 100.0 * clay_content))
#
dz_soil = 0.1   # top soil layer is 10cm  
#

#####################################################################################
# 3. Read in data:
##################################
print 'Reading LST, T_soil and SM_sthu from JULES output data'
print JULES_source_info['file']
inf        = nc.Dataset(JULES_source_info['file'],'r')
# Get lat and lon data
lat        = inf.variables['latitude'][:].squeeze()
lon        = inf.variables['longitude'][:].squeeze()
# get LST, T_soil and SM_sthu
LST        = inf.variables['tstar'][:,:5,:,:].squeeze()
T_soil     = inf.variables['t_soil'][:,0,:,:].squeeze()
SM_sthu    = inf.variables['sthu'][:,0,:,:].squeeze()
inf.close()

# Get sm_wilt as a fraction of saturation
print 'Reading Soil saturation and wilting point data from file'
print SOIL_source_info['file']
inf=nc.Dataset(SOIL_source_info['file'],'r')
SMwilt_sthu= inf.variables[SOIL_source_info['wilt_nc_name']][:].squeeze() / \
             inf.variables[SOIL_source_info['sat_nc_name']][:].squeeze()
# If Jules soil, take top soil layer
if (SOIL_source_info['source']=='JULES_SOIL'):
    SMwilt_sthu = SMwilt_sthu[:,0,:]

inf.close()

# Get LAI data from file
print 'Reading LAI data from file'
print LAI_source_info['file']
LAI        = nc.Dataset(LAI_source_info['file'],'r').\
             variables[LAI_source_info['nc_name']][:]
LAI=LAI.squeeze()
#print 'LAI shape = ',LAI.shape

# Get PFT data from file
print 'Reading PFT data from file'
print PFT_source_info['file']
PFTfrac    = nc.Dataset(PFT_source_info['file'],'r').\
             variables[PFT_source_info['nc_name']][:,:5,:] 
PFTfrac=PFTfrac.squeeze()
#print 'PFT shape = ',PFTfrac.shape

index=np.where(SM_sthu==0)


# Ice points for masking
SM_sthu  = np.ma.masked_equal(SM_sthu,0)
SM_sthu.fill_value=fill_value

ICEmask  = np.ones_like(SM_sthu[0,:])


SMwilt_sthu = SMwilt_sthu*ICEmask
SMwilt_sthu.fill_value=-999.

# Mask all other arrays to match SM
PFTfrac= PFTfrac*ICEmask
PFTfrac.fill_value=fill_value
LAI    = LAI*ICEmask
LAI.fill_value=fill_value
LST    = LST*ICEmask
LST.fill_value=fill_value
T_soil = T_soil*ICEmask
T_soil.fill_value=fill_value


# Calculate dimension lengths
nMONTHs     = JULES_source_info['tsteps']   
nYEARs      = nMONTHs/12.
nPFTs       = LAI.shape[1]
nlandpoints = LAI.shape[2]

######################################################
# 2.9 Read in grid file for plotting and diagnostics
####################################################
if (JULES_source_info['grid']=='WFDEI'):
    gridfile = '/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
    inf=nc.Dataset(gridfile,'r')
    grindex=inf.variables['land_index'][:]-1
    grifrac=inf.variables['land_fraction'][:]
    grimask=np.ones_like(grindex)
    inf.close()
elif (JULES_source_info['grid']=='CHESS'):
    gridfile = '/users/eow/edwcom/CHESS/chess_jules_land_index.nc'
    inf=nc.Dataset(gridfile,'r')
    grindex=inf.variables['index_2D'][:]
    onedex =inf.variables['index_to1D'][:]
    grimask=np.ones_like(grindex)
    inf.close()



#####################################################################################
# 3. Choose tstep option and calculate required setups
######################################
# 3.1 Average
# PFT params
calc_types = ['All_timesteps','Monthly_climatology','Mean_all_timesteps']
iCALC = di_SC.select_source(calc_types,Message=' \nCalculation Types: ')
calc_type=calc_types[iCALC]

out_tag = out_tag+'_'+calc_type

if (calc_type=='All_timesteps'):
    # Temp and SM parameters copy straight to work params
    LST_work = LST
    T_soil_work = T_soil
    SM_sthu_work = SM_sthu
    
    # LAI, if same number of tsteps copy straight over
    if (LAI_source_info['tsteps']==JULES_source_info['tsteps']):
        LAI_work = LAI
    # Else if single time step replicate on to all timesteps
    elif (LAI_source_info['tsteps']==1):
        LAI_work = np.ma.array([LAI.filled() for i in range(nMONTHs)])
        LAI_work=LAI_work*ICEmask
        LAI_work.fill_value=fill_value
    # Else if a 12 month climatology replicate for all years
    elif (LAI_source_info['tsteps']==12):
        LAI_work = np.ma.array([LAI.filled() for i in range(nYEARs)]).\
                        reshape(nMONTHs,nPFTs,nlandpoints)
        LAI_work = LAI_work*ICEmask
        LAI_work.fill_value=fill_value

    # PFT, if same number of tsteps copy straight over
    if (PFT_source_info['tsteps']==PFT_source_info['tsteps']):
        PFTfrac_work = PFTfrac
    # Else if single time step replicate on to all timesteps
    elif (PFT_source_info['tsteps']==1):
        PFTfrac_work = np.ma.array([PFTfrac.filled() for i in range(nMONTHs)])
        PFTfrac_work=PFTfrac_work*ICEmask
        PFTfrac_work.fill_value=fill_value 
    
    # SOIL, if same number of tsteps copy straight over
    if (SOIL_source_info['tsteps']==SOIL_source_info['tsteps']):
        SMwilt_work = SMwilt_sthu
    # Else if single time step replicate on to all timesteps
    elif (SOIL_source_info['tsteps']==1):
        SMwilt_work = np.ma.array([SMwilt_sthu.filled() for i in range(nMONTHs)])
        SMwilt_work=SMwilt_work*ICEmask
        SMwilt_work.fill_value=fill_value 
    

    # Move PFT dim to first position
    LST_work=LST_work.transpose(1,0,2)
    LAI_work=LAI_work.transpose(1,0,2)
    PFTfrac_work=PFTfrac_work.transpose(1,0,2)

    # Create Dummy array of ones for input into soil resp functions
    Carbon_Ones     = np.ones([nPOOLs,nMONTHs,nlandpoints])
    Carbon_Ones_1SP = np.ones([1,nMONTHs,nlandpoints]) 
    
elif (calc_type=='Monthly_climatology'):
    # Temp and SM regrid and then mean over year dimension
    LST_work     = np.mean(LST.reshape(nYEARs,12,nPFTs,nlandpoints),axis=0)
    T_soil_work  = np.mean(T_soil.reshape(nYEARs,12,nlandpoints),axis=0)
    SM_sthu_work = np.mean(SM_sthu.reshape(nYEARs,12,nlandpoints),axis=0)

    # LAI, if same number of tsteps regrid and then mean over year dimension
    if (LAI_source_info['tsteps']==JULES_source_info['tsteps']):
        LAI_work = np.mean(LAI.reshape(nYEARs,12,nPFTs,nlandpoints),axis=0)
    # Else if single time step replicate on to all timesteps
    elif (LAI_source_info['tsteps']==1):
        LAI_work = np.ma.array([LAI.filled() for i in range(12)])
        LAI_work=LAI_work*ICEmask
        LAI_work.fill_value=fill_value
    # Else if already a climatology just copy straight over
    elif (LAI_source_info['tsteps']==12):
        LAI_work = LAI
    
    # PFT, if same number of tsteps regrid and then mean over year dimension
    if (PFT_source_info['tsteps']==JULES_source_info['tsteps']):
        PFTfrac_work = np.mean(PFTfrac.reshape(nYEARs,12,nPFTs,nlandpoints),axis=0)
    # Else if single time step replicate on to all timesteps
    elif (PFT_source_info['tsteps']==1):
        PFTfrac_work = np.ma.array([PFTfrac.filled() for i in range(12)]) 
        PFTfrac_work=PFTfrac_work*ICEmask
        PFTfrac_work.fill_value=fill_value 

    # SOIL, if same number of tsteps regrid and then mean over year dimension
    if (SOIL_source_info['tsteps']==JULES_source_info['tsteps']):
        SMwilt_work = np.mean(SMwilt_sthu.reshape(nYEARs,12,nSOILs,nlandpoints),axis=0)
    # Else if single time step replicate on to all timesteps
    elif (SOIL_source_info['tsteps']==1):
        SMwilt_work = np.ma.array([SMwilt_sthu.filled() for i in range(12)]) 
        SMwilt_work=SMwilt_work*ICEmask
        SMwilt_work.fill_value=fill_value 

    # Move PFT dim to first position
    LST_work=LST_work.transpose(1,0,2)
    LAI_work=LAI_work.transpose(1,0,2)
    PFTfrac_work=PFTfrac_work.transpose(1,0,2)

    # Create Dummy array of ones for input into soil resp functions
    Carbon_Ones=np.ones([nPOOLs,12,nlandpoints])
    Carbon_Ones_1SP = np.ones([1,12,nlandpoints]) 
    
elif (calc_type=='Mean_all_timesteps'):
    # Temp and SM regrid and then mean over year dimension
    LST_work     = np.mean(LST,axis=0)
    T_soil_work  = np.mean(T_soil,axis=0)
    SM_sthu_work = np.mean(SM_sthu,axis=0)
    SMwilt_work  = np.mean(SMwilt_sthu,axis=0)

    # LAI, if same number of tsteps mean over  first dimension
    if (LAI_source_info['tsteps']==JULES_source_info['tsteps']):
        LAI_work = np.mean(LAI,axis=0)
    # Else if single time step replicate on to all timesteps
    elif (LAI_source_info['tsteps']==1):
        LAI_work = LAI
    # Else if 12 month climatology mean over first dimension
    elif (LAI_source_info['tsteps']==12):
        LAI_work = np.mean(LAI,axis=0)
    
    # PFT, if same number of tsteps copy straight over
    if (PFT_source_info['tsteps']==JULES_source_info['tsteps']):
        PFTfrac_work = np.mean(PFTfrac,axis=0)
    # Else if single time step replicate on to all timesteps
    elif (PFT_source_info['tsteps']==1):
        PFTfrac_work = PFTfrac

    # SOIL, if same number of tsteps copy straight over
    if (SOIL_source_info['tsteps']==JULES_source_info['tsteps']):
        SMwilt_work = np.mean(SMwilt_sthu,axis=0)
    # Else if single time step replicate on to all timesteps
    elif (SOIL_source_info['tsteps']==1):
        SMwilt_work = SMwilt_sthu

    # Create Dummy array of ones for input into soil resp functions
    Carbon_Ones=np.ones([nPOOLs,nlandpoints])
    Carbon_Ones_1SP = np.ones([1,nlandpoints]) 
    
# 3.3 calculate veg fraction (sum of all PFTs)
VEGfrac_work = np.sum(PFTfrac_work,axis=0)

#####################################################################################
# 4. Litterfall component
##########################
# 
print 'Litterfall Calculations'
# 4.1 Calculate local litterfall
LocLit_C_ECP, Cveg_ECP = J_SC_F.JULES_LOCAL_LITTERFALL( LAI_work, LST_work,  \
                                                        get_Cv=True,with_phenol=False)
LocLit_C_ECP = LocLit_C_ECP*ICEmask
LocLit_C_ECP.fill_value=fill_value
Cveg_ECP     = Cveg_ECP*ICEmask
Cveg_ECP.fill_value=fill_value
#
# 4.2 Calculate the total litterfall per PFT
Lit_C_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_work, LocLit_C_ECP, Cveg_ECP, \
                                           Per_PFT=True )
Lit_C_ECP=Lit_C_ECP*ICEmask
Lit_C_ECP.fill_value=fill_value
#
# 4.3 Calculate the total litterfall for all PFTs
Lit_C_total_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_work, LocLit_C_ECP, Cveg_ECP )
Lit_C_total_ECP = Lit_C_total_ECP*ICEmask
Lit_C_total_ECP.fill_value = fill_value
#
# 4.4 Split the total litter fall into DPM and RPM components
Lit_SCpools_ECP = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_ECP)
Lit_SCpools_ECP = Lit_SCpools_ECP*ICEmask
Lit_SCpools_ECP.fill_value=fill_value


mnth=6
#plot_LITdiags=False
if (plot_LITdiags==True):
    print Lit_C_ECP.shape
    temp_data=Lit_C_ECP*PFTfrac_work
    for i in range(5):
        plt.subplot(3,2,i+1)
        temp_plotdata=np.max(temp_data[i,:,:],axis=0)[grindex]*grimask
        plt.imshow(temp_plotdata,origin='bottom',vmin=0,vmax=1,cmap='YlOrBr')
        plt.colorbar()
        
    plt.subplot(3,2,6)
    temp_plotdata=np.max(Lit_C_total_ECP,axis=0)[grindex]*grimask
    plt.imshow(Lit_C_total_ECP[mnth,grindex]*grimask,origin='bottom',vmin=0,vmax=1,cmap='YlOrBr')
    plt.colorbar()
    plt.show()

#quit()
#####################################################################################
# 5. Soil Respiration Factor component
#######################################
print 'Soil Respiration Calculations'
# As solving for Soil Carbon we only want the factor of the Soil 
#  respiration function at this point. Simple hack is to call my soil respiration
#  function with Soil Carbon = 1
#  Soil_Resp_Fact= kappa * F(Ts) * F(SM) * F(Vf)
# 5.1 Create Dummy soil carbon array

#
# The input data all have 4 soil layers where the functions only require the surface layer,
# hence we zero index the first dimension of the soil parameters
# 5.2a Soil_Resp_factor for Q10 temperature function
print 'Soil_Resp_Fact_Q10t'

Soil_Resp_Fact_Q10t = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_work,     \
                                                     SM_sthu_work, SMwilt_work,  \
                                                     VEGfrac_work, Tfunc='Q10',    \
                                                     OUTPUT_opt='pools'  )
# multiply by ICEmask
Soil_Resp_Fact_Q10t = Soil_Resp_Fact_Q10t * ICEmask
Soil_Resp_Fact_Q10t.fill_value=fill_value
#
# 5.2b Soil_Resp_factor for RothC temperature function
print 'Soil_Resp_Fact_RothCt'
Soil_Resp_Fact_RothCt = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_work,    \
                                                       SM_sthu_work, SMwilt_work, \
                                                       VEGfrac_work, Tfunc='RothC', \
                                                       OUTPUT_opt='pools')
# multiply by ICEmask
Soil_Resp_Fact_RothCt = Soil_Resp_Fact_RothCt * ICEmask
Soil_Resp_Fact_RothCt.fill_value=fill_value
#

# 5.3 Calculate the Soil_Resp_Fact for a single carbon pool:
print 'Soil_Resp_Fact_singlepool_Q10t'
Soil_Resp_Fact_singlepool_Q10t = \
    J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones_1SP,T_soil_work,\
                                   SM_sthu_work,SMwilt_work, \
                                   VEGfrac_work, Tfunc='Q10',OUTPUT_opt='pools',  \
                                   kappa=kappa_s_1SC )
Soil_Resp_Fact_singlepool_Q10t = Soil_Resp_Fact_singlepool_Q10t.squeeze() * ICEmask
Soil_Resp_Fact_singlepool_Q10t.fill_value=fill_value

# 5.3b RothC temperature Function:
print 'Soil_Resp_Fact_singlepool_RothCt'
Soil_Resp_Fact_singlepool_RothCt = \
    J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones_1SP, T_soil_work,\
                                   SM_sthu_work,SMwilt_work, \
                                   VEGfrac_work, Tfunc='RothC',OUTPUT_opt='pools', \
                                    kappa=kappa_s_1SC )
Soil_Resp_Fact_singlepool_RothCt = Soil_Resp_Fact_singlepool_RothCt.squeeze() * ICEmask
Soil_Resp_Fact_singlepool_RothCt.fill_value=fill_value
#


#plot_SOILdiags=False
if (plot_SOILdiags==True):
    Lit_SCpools_ECP.shape
    temp_data=Lit_SCpools_ECP
    for i in range(2):
        plt.subplot(3,2,i+1)
        temp_plotdata=np.max(temp_data[i,:,:],axis=0)[grindex]*grimask
        plt.imshow(temp_plotdata,origin='bottom',vmin=0,vmax=1,cmap='YlOrBr')
        plt.colorbar()
    
    print Soil_Resp_Fact_Q10t.shape
    temp_data=Soil_Resp_Fact_Q10t*3600.*24.*365
    for i in range(4):
        plt.subplot(3,2,i+3)
        temp_plotdata=np.max(temp_data[i,:,:],axis=0)[grindex]*grimask
        plt.imshow(temp_plotdata,origin='bottom',vmin=0,vmax=5,cmap='YlOrBr')
        plt.colorbar()
        
    plt.show()


######################################################################################
# 6. Calculate the soil carbon for each pool using the albmar/ECP equations
#########################################################################
print 'RothC Soil Carbon Pools Calculations'
#
# 6.1 Decomposable Plant Material
#  C_dpm = (f_dpm*Lit_c)/ [Soil_resp_factor]_dpm
#      Calculate Cpool based on Q10 temp function and RothC temp function
#      Correcting for unit difference, i.e. seconds in year
#print Lit_SCpools_ECP.shape
#print Soil_Resp_Fact_Q10t.shape
C_dpm_Q10t   = (np.sum(Lit_SCpools_ECP[0,:,:], axis=0)       /  \
                np.sum(Soil_Resp_Fact_Q10t[0,:,:], axis=0) ) /  \
                (3600.*24.*360.)


C_dpm_RothCt = (np.sum(Lit_SCpools_ECP[0,:,:], axis=0)       /  \
                np.sum(Soil_Resp_Fact_RothCt[0,:,:], axis=0) ) /  \
                (3600.*24.*360.)

#
# 6.2 Resistant Plant Material:
#  C_rpm = (f_dpm*Lit_c)/ [Soil_resp_factor]_dpm
#      Calculate Cpool based on Q10 temp function and RothC temp function
C_rpm_Q10t   = (np.sum(Lit_SCpools_ECP[1,:,:], axis=0)       /  \
                       np.sum(Soil_Resp_Fact_Q10t[1,:,:], axis=0) ) /  \
                       (3600.*24.*360.)

C_rpm_RothCt = (np.sum(Lit_SCpools_ECP[1,:,:], axis=0)       /  \
                       np.sum(Soil_Resp_Fact_RothCt[1,:,:], axis=0) ) /  \
                       (3600.*24.*360.)

#
# See Eddy Equations:
# 6.3 Biomass:
#  C_bio =  ( 0.46 / ( (1/beta_r)-1 )*kappa_bio ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
# 6.4 Hummus:
#  C_hum =  ( 0.54 / ( (1/beta_r)-1 )*kappa_hum ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
#      Calculate Cpool based on Q10 temp function and RothC temp function
dpm_rpm_term_Q10t   = (kappa_s_RothC[0]*C_dpm_Q10t)+(kappa_s_RothC[1]*C_rpm_Q10t)
dpm_rpm_term_RothCt = (kappa_s_RothC[0]*C_dpm_RothCt)+(kappa_s_RothC[1]*C_rpm_RothCt)
bio_factor =  0.46 / ( ( (1/beta_r)-1 )*kappa_s_RothC[2] )
hum_factor =  0.54 / (( (1/beta_r)-1 )*kappa_s_RothC[3] )
# 6.3:
C_bio_Q10t   = bio_factor*dpm_rpm_term_Q10t
C_bio_RothCt = bio_factor*dpm_rpm_term_RothCt
# 6.4:
C_hum_Q10t   = hum_factor*dpm_rpm_term_Q10t
C_hum_RothCt = hum_factor*dpm_rpm_term_RothCt
#
# 6.5 append to single array:
C_4pools_Q10t   = np.array( [C_dpm_Q10t,C_rpm_Q10t,C_bio_Q10t,C_hum_Q10t] )*ICEmask
C_4pools_RothCt = np.array( [C_dpm_RothCt,C_rpm_RothCt,C_bio_RothCt,C_hum_RothCt] )*ICEmask
#

C_4pools_Q10temp=np.ma.masked_equal(C_4pools_Q10t,0)
C_4pools_RothCtemp=np.ma.masked_equal(C_4pools_Q10t,0)

vmaxi = [1,20,5,50]
#plot_Cq10diags=False
if (plot_Cq10diags==True):
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.imshow(C_4pools_Q10t[i,grindex-1]*grimask,\
                   origin='bottom',cmap='YlOrBr')
        plt.colorbar()
    plt.show()
    

if (plot_Crothdiags==True):
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.imshow((C_4pools_RothCt[i,grindex-1]*grimask), \
                 origin='bottom',cmap='YlOrBr')
        plt.colorbar()
    plt.show()

######################################################################################
# 7. Calculate the soil carbon for a single pool using the albmar/ECP equations
#########################################################################
print 'Single Soil Carbon Pool Calculations'
#
# 7.1 
# Carbon for single pool = total litterfall / soil respiration factor
C_1SC_Q10t   =  (np.sum(Lit_C_total_ECP, axis=0)       /  \
                 np.sum(Soil_Resp_Fact_singlepool_Q10t, axis=0) ) /  \
                 (3600.*24.*360.)
C_1SC_Q10temp = np.ma.masked_equal(C_1SC_Q10t,0)

C_1SC_RothCt =  (np.sum(Lit_C_total_ECP, axis=0)       /  \
                 np.sum(Soil_Resp_Fact_singlepool_RothCt, axis=0) ) /  \
                 (3600.*24.*360.)
C_1SC_RothCtemp = np.ma.masked_equal(C_1SC_RothCt,0)

#plot_Cq10diags=True
if (plot_Cq10diags==True):
    plt.subplot(2,1,1)
    plt.imshow(np.sum(C_4pools_Q10t,axis=0)[grindex-1]*grimask, \
               origin='bottom',cmap='YlOrBr')
    plt.colorbar()
    print 'C_1SC_Q10t.shape =', C_1SC_Q10t.shape
    plt.subplot(2,1,2)
    plt.imshow(C_1SC_Q10t[grindex-1]*grimask, \
               origin='bottom',cmap='YlOrBr')
    plt.colorbar()
    plt.show()
#quit()
######################################################################################
# 9. Output the data to netCDF files
#########################################################################
#
print 'Output data to netCDF'
outfile =  OUT_DIR+out_tag+'.nc'
print outfile
# 9.1 open output file for 4 pools data.
outf=nc.Dataset(outfile,'w')

# else output a single snap shot
outf.createDimension('x',nlandpoints)
outf.createDimension('month',nMONTHs)
outf.createDimension('pool',nPOOLs)
outf.createDimension('pft',nPFTs)
outf.createDimension('LitSCpool',2)
#
outvar=outf.createVariable('lats','float32',('x'))
outvar[:]=lat
outvar=outf.createVariable('lons','float32',('x'))
outvar[:]=lon

# LST
outvar=outf.createVariable('LST','float32',('pft','month','x'))
outvar.units='K'
outvar.note='LST for the 5 standard JULES PFTs'
outvar[:]=LST_work  

# Monthly T_soil
outvar=outf.createVariable('T_SOIL','float32',('month','x'))
outvar.units='K'
outvar.note='Monthly Soil Temperature (Top layer 10cm)'
outvar[:]=T_soil_work  

# Monthly Soil Moisture
outvar=outf.createVariable('SM','float32',('month','x'))
outvar.units='-'
outvar.note='Monthly Soil Moisture as fraction of saturation (Top layer 10cm)'
outvar[:]=SM_sthu_work  

#  Soil Q10t Resp Factor
outvar=outf.createVariable('Soil_Resp_Fact_Q10t','float32',('pool','month','x'))
outvar.units='kgC m-2 s-1'
outvar.note='Monthly Soil Respiration Factor using Q10 temperature function'
outvar[:]=Soil_Resp_Fact_Q10t

#  Soil RothCt Resp Factor
outvar=outf.createVariable('Soil_Resp_Fact_RothCt','float32',('pool','month','x'))
outvar.units='kgC m-2 s-1'
outvar.note='Monthly Soil Respiration Factor using RothC temperature function'
outvar[:]=Soil_Resp_Fact_RothCt

# Local Litterfall
outvar=outf.createVariable('Local_Litterfall_C','float32',('pft','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Local Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=LocLit_C_ECP

# Total Litterfall per PFT
outvar=outf.createVariable('Litterfall_C','float32',('pft','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=Lit_C_ECP
    
# Total Litterfall
outvar=outf.createVariable('Total_Litterfall_C','float32',('month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon summed over all 5 standard JULES PFTs'
outvar[:]=Lit_C_total_ECP

# Lit_SCpools
outvar=outf.createVariable('Lit_SCpools','float32',('LitSCpool','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon that goes into each Carbon Pool'
outvar[:]= Lit_SCpools_ECP
    
# Cveg per PFT
outvar=outf.createVariable('c_veg','float32',('pft','month','x'))
outvar.units='kgC m-2'
outvar.note='Total vegetation Carbon for each of the 5 standard JULES PFTs'
outvar[:]=Cveg_ECP
    
# PFT fraction
outvar=outf.createVariable('PFT_frac','float32',('pft','month','x'))
outvar.units='kgC m-2'
outvar.note='Fractional cover of the 5 standard JULES PFTs'
outvar[:]=PFTfrac_work   
    
# LAI
outvar=outf.createVariable('LAI','float32',('pft','month','x'))
outvar.units='m2 / m2'
outvar.note='LAI of the 5 standard JULES PFTs'
outvar[:]=LAI_work  
    
# Carbon 4 pools Q10t
outvar=outf.createVariable('C_4pools_Q10t','float32',('pool','x'))
outvar.units='kgC m-2'
outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
outvar[:]=C_4pools_Q10t
    
# Carbon 4 pools RothC
outvar=outf.createVariable('C_4pools_RothCt','float32',('pool','x'))
outvar.units='kgC m-2'
outvar.note='Soil Carbon for the RothC pools calculated using the RothC temperature function'
outvar[:]=C_4pools_RothCt
    
# Carbon single pool Q10
outvar=outf.createVariable('C_1SC_Q10t','float32',('x'))
outvar.units='kgC m-2'
outvar.note='Soil Carbon for a single carbon pool calculated using the Q10 temperature function'
outvar[:]=C_1SC_Q10t

# Carbon single pool RothC
outvar=outf.createVariable('C_1SC_RothCt','float32',('x'))
outvar.units='kgC m-2'
outvar.note='Soil Carbon for a single carbon pool calculated using the RothC temperature function'
outvar[:]=C_1SC_RothCt
#

outf.author    = 'Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.JULES_RUN = JULES_source 
outf.LAI_source= LAI_source_info['source']
outf.PFT_source= PFT_source_info['source']
outf.calc_type = calc_type
outf.close()



