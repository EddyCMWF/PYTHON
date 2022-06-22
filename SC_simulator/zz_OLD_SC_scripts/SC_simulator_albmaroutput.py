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
#

#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_frac_file   = WFDEI_dir+'frac_igbp_watch_0p5deg_capUM6.6_E2OBS.nc'
WFDEI_soil_file   = WFDEI_dir+'soil_igbp_bc_watch_0p5deg_capUM6.6_2D_EMMA.nc'

WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'
JULES_output_dir  = WFDEI_dir+'JULES_output/phen/'
JULES_file_tag    = 'e2o_nerc_phen_glob30_PARAM.nc'
#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
out_tag = 'SCsim_albmar_JULES-WFDEI-1979-2012'
#

#####################################################################################
# 2. Define or Read in the various data:
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
#
# List of Parameters required from JULES output:
#JULES_params= [ 'LAI', 'AvgSurfT', 'TotMoist', 'TotSoilSat' ]

# 2.4 Read in T_soil, SM, SC, LAI, CanHt and PFT from JULES output
print 'Reading in JULES output data'
#JULES_data_dict = {}
LAI   = nc.Dataset(JULES_output_dir+JULES_file_tag.replace('PARAM','mon_LAI_1979-2012'),\
                   'r').variables['LAI'][:]
LST   = nc.Dataset(JULES_output_dir+JULES_file_tag.replace('PARAM','mon_AvgSurfT_1979-2012'),\
                   'r').variables['AvgSurfT'][:]
SM    = nc.Dataset(JULES_output_dir+JULES_file_tag.replace('PARAM','mon_TotMoist_1979-2012'),\
                   'r').variables['TotMoist'][:]
SMsat = nc.Dataset(JULES_output_dir+JULES_file_tag.replace('PARAM','fix_TotSoilSat'),\
                   'r').variables['TotSoilSat'][:]

# Soil moisture as a fraction of saturation:
SM_STHF = SM/SMsat

# 2.5 read in pft fractions and sm_wilt from 2D ancil files
PFTfrac = nc.Dataset(WFDEI_frac_file,'r').variables['frac'][:5,:,:]

SMwilt_V = np.zeros_like(SMsat)
SMsat_V = np.zeros_like(SMsat)
# The following indices were calculated offline, 
# places subsection into correct location
SMwilt_V[68:348] = nc.Dataset(WFDEI_soil_file,'r').variables['vwilt'][:]
SMsat_V[68:348] = nc.Dataset(WFDEI_soil_file,'r').variables['vsat'][:]

SMwilt_sthf = SMwilt_V/SMsat_V 

#####################################################################################
# 3. Calculate/Create monthly climatologies
######################################
LAI_climat     = np.mean(LAI.reshape(34,12,360,720),axis=0)
LST_climat     = np.mean(LST.reshape(34,12,360,720),axis=0)
SM_STHF_climat = np.mean(SM_STHF.reshape(34,12,360,720),axis=0)

PFTfrac_climat = 

#####################################################################################
# 4. Litterfall component
##########################
# 
print 'Litterfall Calculations'
# 4.1 Calculate local litterfall
LocLit_C_ECP, Cveg_ECP = J_SC_F.JULES_LOCAL_LITTERFALL( LAI, LST, \
                                                        LAI_b=LAI_BAL_JULES,    \
                                                        get_Cv=True,with_phenol=False     )
#
# 4.2 Calculate the total litterfall per PFT
Lit_C_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_JULES, LocLit_C_ECP, Cveg_ECP, \
                                           Per_PFT=True )
#
# 4.3 Calculate the total litterfall for all PFTs
Lit_C_total_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_JULES, LocLit_C_ECP, Cveg_ECP )
#
# 4.4 Split the total litter fall into DPM and RPM components
Lit_SCpools_ECP = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_ECP)

#####################################################################################
# 5. Soil Respiration Factor component
#######################################
print 'Soil Respiration Calculations'
# As solving for Soil Carbon we only want the factor of the Soil 
#  respiration function at this point. Simple hack is to call my soil respiration
#  function with Soil Carbon = 1
#  Soil_Resp_Fact= kappa * F(Ts) * F(SM) * F(Vf)
# 5.1 Create Dummy soil carbon array
Carbon_Ones=np.ones_like(T_soil_JULES)
#
# The input data all have 4 soil layers where the functions only require the surface layer,
# hence we zero index the first dimension of the soil parameters
# 5.2a Soil_Resp_factor for Q10 temperature function
print 'Soil_Resp_Fact_Q10t'
Soil_Resp_Fact_Q10t = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_JULES[0,:],           \
                                                     sthu_JULES[0,:], sm_wilt_frac_JULES[0,:], \
                                                     FracVeg, Tfunc='Q10', OUTPUT_opt='pools'  )
#
# 5.2b Soil_Resp_factor for RothC temperature function
print 'Soil_Resp_Fact_RothCt'
Soil_Resp_Fact_RothCt = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_JULES[0,:],           \
                                                       sthu_JULES[0,:], sm_wilt_frac_JULES[0,:], \
                                                       FracVeg, Tfunc='RothC', OUTPUT_opt='pools')
#
#
# 5.3 Calculate the Soil_Resp_Fact for a single carbon pool:
# Create temporary dimension list to maintain a single SC pool as first dim.
temp_CarbonShape=[1]
for dim in Carbon_Ones[0,:].shape:
    temp_CarbonShape.append(dim)
# 5.3a Q10 temperature Function:
print 'Soil_Resp_Fact_singlepool_Q10t'
Soil_Resp_Fact_singlepool_Q10t = \
    J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones[0,:].reshape(temp_CarbonShape), \
                                   T_soil_JULES[0,:],sthu_JULES[0,:], sm_wilt_frac_JULES[0,:], \
                                   FracVeg, Tfunc='Q10',OUTPUT_opt='pools',  \
                                   kappa=kappa_s_1SC )
# 5.3b RothC temperature Function:
print 'Soil_Resp_Fact_singlepool_RothCt'
Soil_Resp_Fact_singlepool_RothCt = \
    J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones[0,:].reshape(temp_CarbonShape), \
                                   T_soil_JULES[0,:],sthu_JULES[0,:], sm_wilt_frac_JULES[0,:],  \
                                    FracVeg, Tfunc='RothC',OUTPUT_opt='pools', \
                                    kappa=kappa_s_1SC )
#

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
C_dpm_Q10t   = (Lit_SCpools_ECP[0,:] / Soil_Resp_Fact_Q10t[0,:]) / (3600.*24.*360.)
C_dpm_RothCt = (Lit_SCpools_ECP[0,:] / Soil_Resp_Fact_RothCt[0,:]) / (3600.*24.*360.)

#
# 6.2 Resistant Plant Material:
#  C_rpm = (f_dpm*Lit_c)/ [Soil_resp_factor]_dpm
#      Calculate Cpool based on Q10 temp function and RothC temp function
C_rpm_Q10t   = (Lit_SCpools_ECP[1,:] / Soil_Resp_Fact_Q10t[1,:]) / (3600.*24.*360.)
C_rpm_RothCt = (Lit_SCpools_ECP[1,:] / Soil_Resp_Fact_RothCt[1,:]) / (3600.*24.*360.)
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
C_4pools_Q10t   = np.array( [C_dpm_Q10t,C_rpm_Q10t,C_bio_Q10t,C_hum_Q10t] )
C_4pools_RothCt = np.array( [C_dpm_RothCt,C_rpm_RothCt,C_bio_RothCt,C_hum_RothCt] )
#

######################################################################################
# 7. Calculate the soil carbon for a single pool using the albmar/ECP equations
#########################################################################
print 'Single Soil Carbon Pool Calculations'
#
# 7.1 
# Carbon for single pool = total litterfall / soil respiration factor
C_1SC_Q10t   = (Lit_C_total_ECP / Soil_Resp_Fact_singlepool_Q10t) / (3600.*24.*365.)
C_1SC_RothCt = (Lit_C_total_ECP / Soil_Resp_Fact_singlepool_RothCt) / (3600.*24.*365.)
#
#print C_1SC_Q10t.shape
#print C_1SC_RothCt.shape

######################################################################################
# 8. Put soil carbon onto original land points
##############################################
# If outputting a monthly data time dimension required
if 'monthly' in INTERVAL:
    # 8.1 Create arrays based on dimensions of original lat array
    # Soil Carbon 
    C_4pools_Q10t_origgrid=np.zeros([4,12,lats_orig.shape[0]])+fill_value
    C_4pools_RothCt_origgrid=np.zeros([4,12,lats_orig.shape[0]])+fill_value
    C_1SC_Q10t_origgrid=np.zeros([1,12,lats_orig.shape[0]])+fill_value
    C_1SC_RothCt_origgrid=np.zeros([1,12,lats_orig.shape[0]])+fill_value
    # Supplementary data
    Soil_Resp_Fact_Q10t_origgrid=np.zeros([4,12,lats_orig.shape[0]])+fill_value
    Soil_Resp_Fact_RothCt_origgrid=np.zeros([4,12,lats_orig.shape[0]])+fill_value
    LocLit_C_ECP_origgrid=np.zeros([5,12,lats_orig.shape[0]])+fill_value
    Lit_C_ECP_origgrid=np.zeros([5,12,lats_orig.shape[0]])+fill_value
    Lit_C_total_ECP_origgrid=np.zeros([12,lats_orig.shape[0]])+fill_value
    Lit_SCpools_origgrid=np.zeros([4,12,lats_orig.shape[0]])+fill_value
    Cveg_ECP_origgrid=np.zeros([5,12,lats_orig.shape[0]])+fill_value
    PFTfrac_origgrid=np.zeros([5,12,lats_orig.shape[0]])+fill_value
    LAI_origgrid=np.zeros([5,12,lats_orig.shape[0]])+fill_value
    # 
    # 
    # 8.2 Fill arrays with data and mask
    C_4pools_Q10t_origgrid[:,:,soil_points]=C_4pools_Q10t
    C_4pools_Q10t_origgrid=np.ma.masked_equal(C_4pools_Q10t_origgrid,fill_value)
    #
    C_4pools_RothCt_origgrid[:,:,soil_points]=C_4pools_RothCt
    C_4pools_RothCt_origgrid=np.ma.masked_equal(C_4pools_RothCt_origgrid,fill_value)
    #
    C_1SC_Q10t_origgrid[:,:,soil_points]=C_1SC_Q10t
    C_1SC_Q10t_origgrid=np.ma.masked_equal(C_1SC_Q10t_origgrid,fill_value)
    #
    C_1SC_RothCt_origgrid[:,:,soil_points]=C_1SC_RothCt
    C_1SC_RothCt_origgrid=np.ma.masked_equal(C_1SC_RothCt_origgrid,fill_value)
    #
    # Supplementary data
    Soil_Resp_Fact_Q10t_origgrid[:,:,soil_points]=Soil_Resp_Fact_Q10t
    Soil_Resp_Fact_Q10t_origgrid=np.ma.masked_equal(Soil_Resp_Fact_Q10t_origgrid,fill_value)
    #
    Soil_Resp_Fact_RothCt_origgrid[:,:,soil_points]=Soil_Resp_Fact_RothCt
    Soil_Resp_Fact_RothCt_origgrid=np.ma.masked_equal(Soil_Resp_Fact_RothCt_origgrid,fill_value)
    #
    LocLit_C_ECP_origgrid[:,:,soil_points]=LocLit_C_ECP
    LocLit_C_ECP_origgrid=np.ma.masked_equal(LocLit_C_ECP_origgrid,fill_value)
    #
    Lit_C_ECP_origgrid[:,:,soil_points]=Lit_C_ECP
    Lit_C_ECP_origgrid=np.ma.masked_equal(Lit_C_ECP_origgrid,fill_value)
    #
    Lit_SCpools_origgrid[:2,:,soil_points]=Lit_SCpools_ECP
    Lit_SCpools_origgrid=np.ma.masked_equal(Lit_SCpools_origgrid,fill_value)
    #
    Lit_C_total_ECP_origgrid[:,soil_points]=Lit_C_total_ECP
    Lit_C_total_ECP_origgrid=np.ma.masked_equal(Lit_C_total_ECP_origgrid,fill_value)
    #
    Cveg_ECP_origgrid[:,:,soil_points]=Cveg_ECP
    Cveg_ECP_origgrid=np.ma.masked_equal(Cveg_ECP_origgrid,fill_value)
    #
    PFTfrac_origgrid[:,:,soil_points]=PFTfrac_JULES
    PFTfrac_origgrid=np.ma.masked_equal(PFTfrac_origgrid,fill_value)
    #
    LAI_origgrid[:,:,soil_points]=LAI_JULES
    LAI_origgrid=np.ma.masked_equal(LAI_origgrid,fill_value)

else:
    # 8.1 Create arrays based on dimensions of original lat array
    C_4pools_Q10t_origgrid=np.zeros([4,lats_orig.shape[0]])+fill_value
    C_4pools_RothCt_origgrid=np.zeros([4,lats_orig.shape[0]])+fill_value
    Soil_Resp_Fact_Q10t_origgrid=np.zeros([4,lats_orig.shape[0]])+fill_value
    Soil_Resp_Fact_RothCt_origgrid=np.zeros([4,lats_orig.shape[0]])+fill_value
    LocLit_C_ECP_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
    Lit_C_ECP_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
    Lit_C_total_ECP_origgrid=np.zeros([lats_orig.shape[0]])+fill_value
    Lit_SCpools_origgrid=np.zeros([4,lats_orig.shape[0]])+fill_value
    Cveg_ECP_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
    PFTfrac_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
    LAI_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
    # 
    C_1SC_Q10t_origgrid=np.zeros(lats_orig.shape[0])+fill_value
    C_1SC_RothCt_origgrid=np.zeros(lats_orig.shape[0])+fill_value
    #
    # 8.2 Fill arrays with data and mask
    C_4pools_Q10t_origgrid[:,soil_points]=C_4pools_Q10t
    C_4pools_Q10t_origgrid=np.ma.masked_equal(C_4pools_Q10t_origgrid,fill_value)
    #
    C_4pools_RothCt_origgrid[:,soil_points]=C_4pools_RothCt
    C_4pools_RothCt_origgrid=np.ma.masked_equal(C_4pools_RothCt_origgrid,fill_value)
    #
    C_1SC_Q10t_origgrid[soil_points]=C_1SC_Q10t
    C_1SC_Q10t_origgrid=np.ma.masked_equal(C_1SC_Q10t_origgrid,fill_value)
    #
    C_1SC_RothCt_origgrid[soil_points]=C_1SC_RothCt
    C_1SC_RothCt_origgrid=np.ma.masked_equal(C_1SC_RothCt_origgrid,fill_value)
    #
    # Supplementary data
    Soil_Resp_Fact_Q10t_origgrid[:,soil_points]=Soil_Resp_Fact_Q10t
    Soil_Resp_Fact_Q10t_origgrid=np.ma.masked_equal(Soil_Resp_Fact_Q10t_origgrid,fill_value)
    #
    Soil_Resp_Fact_RothCt_origgrid[:,soil_points]=Soil_Resp_Fact_RothCt
    Soil_Resp_Fact_RothCt_origgrid=np.ma.masked_equal(Soil_Resp_Fact_RothCt_origgrid,fill_value)
    #
    LocLit_C_ECP_origgrid[:,soil_points]=LocLit_C_ECP
    LocLit_C_ECP_origgrid=np.ma.masked_equal(LocLit_C_ECP_origgrid,fill_value)
    #
    Lit_C_ECP_origgrid[:,soil_points]=Lit_C_ECP
    Lit_C_ECP_origgrid=np.ma.masked_equal(Lit_C_ECP_origgrid,fill_value)
    #
    Lit_C_total_ECP_origgrid[soil_points]=Lit_C_total_ECP
    Lit_C_total_ECP_origgrid=np.ma.masked_equal(Lit_C_total_ECP_origgrid,fill_value)
    #
    Lit_SCpools_origgrid[:2,soil_points]=Lit_SCpools_ECP
    Lit_SCpools_origgrid=np.ma.masked_equal(Lit_SCpools_origgrid,fill_value)
    #
    Cveg_ECP_origgrid[:,soil_points]=Cveg_ECP
    Cveg_ECP_origgrid=np.ma.masked_equal(Cveg_ECP_origgrid,fill_value)
    #
    PFTfrac_origgrid[:,soil_points]=PFTfrac_JULES
    PFTfrac_origgrid=np.ma.masked_equal(PFTfrac_origgrid,fill_value)
    #
    LAI_origgrid[:,soil_points]=LAI_JULES
    LAI_origgrid=np.ma.masked_equal(LAI_origgrid,fill_value)
    #
#

######################################################################################
# 9. Output the data to netCDF files
#########################################################################
#
print 'Output data to netCDF'
# 9.1 open output file for 4 pools data.
outf=nc.Dataset(OUT_DIR+out_tag+'_FastSoilCarbon.nc','w')
# 
# If outputting a monthly data time dimension required
if 'monthly' in INTERVAL:
    # 9.2 declare dimensions
    outf.createDimension('month',12)
    outf.createDimension('x',lats_orig.shape[0])
    outf.createDimension('pool',4)
    outf.createDimension('pft',5)
    #
    outvar=outf.createVariable('lats','float32',('x'))
    outvar[:]=lats_orig
    outvar=outf.createVariable('lons','float32',('x'))
    outvar[:]=lons_orig
    #
    #  Soil Q10t Resp Factor
    outvar=outf.createVariable('Soil_Resp_Fact_Q10t','float32',('pool','month','x'))
    outvar.units='kgC m-2 s-1'
    outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
    outvar[:]=Soil_Resp_Fact_Q10t_origgrid
    #
    #  Soil RothCt Resp Factor
    outvar=outf.createVariable('Soil_Resp_Fact_RothCt','float32',('pool','month','x'))
    outvar.units='kgC m-2 s-1'
    outvar.note='Soil Carbon for the RothC pools calculated using the RothC temperature function'
    outvar[:]=Soil_Resp_Fact_RothCt_origgrid
    #
    # Local Litterfall
    outvar=outf.createVariable('Local_Litterfall_C','float32',('pft','month','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Local Litterfall Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=LocLit_C_ECP_origgrid
    #
    # Total Litterfall per PFT
    outvar=outf.createVariable('Litterfall_C','float32',('pft','month','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=Lit_C_ECP_origgrid
    #
    # Total Litterfall
    outvar=outf.createVariable('Total_Litterfall_C','float32',('month','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon summed over all 5 standard JULES PFTs'
    outvar[:]=Lit_C_total_ECP_origgrid
    #
    # Lit_SCpools
    outvar=outf.createVariable('Lit_SCpools','float32',('pool','month','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon that goes into each Carbon Pool'
    outvar[:]=Lit_SCpools_origgrid
    #
    # Cveg per PFT
    outvar=outf.createVariable('c_veg','float32',('pft','month','x'))
    outvar.units='kgC m-2'
    outvar.note='Total vegetation Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=Cveg_ECP_origgrid
    #
    # PFT fraction
    outvar=outf.createVariable('PFT_frac','float32',('pft','month','x'))
    outvar.units='-'
    outvar.note='Fractional cover of the 5 standard JULES PFTs'
    outvar[:]=PFTfrac_origgrid
    #
    # LAI
    outvar=outf.createVariable('LAI','float32',('pft','month','x'))
    outvar.units='m2 / m2'
    outvar.note='LAI of the 5 standard JULES PFTs'
    outvar[:]=LAI_origgrid
    #
    outvar=outf.createVariable('C_4pools_Q10t','float32',('pool','month','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
    outvar[:]=C_4pools_Q10t_origgrid
    #
    outvar=outf.createVariable('C_4pools_RothCt','float32',('pool','month','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for the RothC pools calculated using the RothC temperature function'
    outvar[:]=C_4pools_RothCt_origgrid
    #
    outvar=outf.createVariable('C_1SC_Q10t','float32',('month','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for a single carbon pool calculated using the Q10 temperature function'
    outvar[:]=C_1SC_Q10t_origgrid
    #
    outvar=outf.createVariable('C_1SC_RothCt','float32',('month','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for a single carbon pool calculated using the RothC temperature function'
    outvar[:]=C_1SC_RothCt_origgrid

# else output a single snap shot
else:
    # 9.2b declare dimensions
    outf.createDimension('x',lats_orig.shape[0])
    outf.createDimension('pool',4)
    outf.createDimension('pft',5)
    #
    outvar=outf.createVariable('lats','float32',('x'))
    outvar[:]=lats_orig
    outvar=outf.createVariable('lons','float32',('x'))
    outvar[:]=lons_orig
    #
    #  Soil Q10t Resp Factor
    outvar=outf.createVariable('Soil_Resp_Fact_Q10t','float32',('pool','x'))
    outvar.units='kgC m-2 s-1'
    outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
    outvar[:]=Soil_Resp_Fact_Q10t_origgrid
    #
    #  Soil RothCt Resp Factor
    outvar=outf.createVariable('Soil_Resp_Fact_RothCt','float32',('pool','x'))
    outvar.units='kgC m-2 s-1'
    outvar.note='Soil Carbon for the RothC pools calculated using the RothC temperature function'
    outvar[:]=Soil_Resp_Fact_RothCt_origgrid
    #
    # Local Litterfall
    outvar=outf.createVariable('Local_Litterfall_C','float32',('pft','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Local Litterfall Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=LocLit_C_ECP_origgrid
    #
    # Total Litterfall per PFT
    outvar=outf.createVariable('Litterfall_C','float32',('pft','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=Lit_C_ECP_origgrid
    #
    # Total Litterfall
    outvar=outf.createVariable('Total_Litterfall_C','float32',('x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon summed over all 5 standard JULES PFTs'
    outvar[:]=Lit_C_total_ECP_origgrid
    #
    # Lit_SCpools
    outvar=outf.createVariable('Lit_SCpools','float32',('pool','x'))
    outvar.units='kgC m-2 year-1'
    outvar.note='Total Litterfall Carbon that goes into each Carbon Pool'
    outvar[:]=Lit_SCpools_origgrid
    #
    # Cveg per PFT
    outvar=outf.createVariable('c_veg','float32',('pft','x'))
    outvar.units='kgC m-2'
    outvar.note='Total vegetation Carbon for each of the 5 standard JULES PFTs'
    outvar[:]=Cveg_ECP_origgrid
    #
    # PFT fraction
    outvar=outf.createVariable('PFT_frac','float32',('pft','x'))
    outvar.units='kgC m-2'
    outvar.note='Fractional cover of the 5 standard JULES PFTs'
    outvar[:]=PFTfrac_origgrid
    #
    # LAI
    outvar=outf.createVariable('LAI','float32',('pft','x'))
    outvar.units='m2 / m2'
    outvar.note='LAI of the 5 standard JULES PFTs'
    outvar[:]=LAI_origgrid
    #
    outvar=outf.createVariable('C_4pools_Q10t','float32',('pool','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
    outvar[:]=C_4pools_Q10t_origgrid
    #
    outvar=outf.createVariable('C_4pools_RothCt','float32',('pool','x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for the RothC pools calculated using the RothC temperature function'
    outvar[:]=C_4pools_RothCt_origgrid
    #
    outvar=outf.createVariable('C_1SC_Q10t','float32',('x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for a single carbon pool calculated using the Q10 temperature function'
    outvar[:]=C_1SC_Q10t_origgrid
    #
    outvar=outf.createVariable('C_1SC_RothCt','float32',('x'))
    outvar.units='kgC m-2'
    outvar.note='Soil Carbon for a single carbon pool calculated using the RothC temperature function'
    outvar[:]=C_1SC_RothCt_origgrid
#

outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.close()

quit()


colors=['y','g','orange','saddlebrown']
labels=['DPM','RPM','BIO','HUM']

fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothCt[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothCt[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()



fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10t[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10t[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()



lat_binsize=5.
minlat=-60.
maxlat=90.
lat_bins=np.arange(minlat,maxlat,lat_binsize)

Lit_SCpools_Zonal=np.zeros([2,len(lat_bins)])
Soil_Resp_Fact_Q10_Zonal=np.zeros([4,len(lat_bins)])
Soil_Resp_Fact_RothC_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_Q10_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_RothC_Zonal=np.zeros([4,len(lat_bins)])


for bin in range(len(lat_bins)):
    tempdex=np.where( (lats_JULES>=lat_bins[bin]        ) & \
                    (lats_JULES< lat_bins[bin]+lat_binsize) )
    #
    Lit_SCpools_Zonal[:,bin]=np.mean(Lit_SCpools_ECP[:,tempdex],axis=1)
    Soil_Resp_Fact_Q10_Zonal[:,bin]=np.mean(Soil_Resp_Fact_Q10t[:,tempdex],axis=1)
    Soil_Resp_Fact_RothC_Zonal[:,bin]=np.mean(Soil_Resp_Fact_RothCt[:,tempdex],axis=1)
    C_4pools_Q10_Zonal[:,bin]=np.mean(C_4pools_Q10t[:,tempdex],axis=1)
    C_4pools_RothC_Zonal[:,bin]=np.mean(C_4pools_RothCt[:,tempdex],axis=1)




fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
#ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()



fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothC_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothC_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()

