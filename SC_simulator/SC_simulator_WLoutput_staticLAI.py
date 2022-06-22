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
# 0. Set flags and diags
###########################
nMONTHs=12
nPOOLs=4

#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:
WFDEI_dir         = '/users/eow/edwcom/WFD_EI/'
WFDEI_frac_file   = WFDEI_dir+'qrparm.veg.fracNew.nc'
frac_varname      = 'field1391'
WFDEI_soil_file   = WFDEI_dir+'qrparm.soil_HWSD_class3_van_genuchtenNew.nc'
smwilt_varname    = 'field329'
smsat_varname     = 'field332'
WFDEI_gridfile    = WFDEI_dir+'wfdei-land-mask.nc'

JULES_output_dir  = WFDEI_dir+'JULES_output/'
JULES_file        = JULES_output_dir+'JULES_v42_WFD-EI-GPCC_Zinke_global_DD_newtopo.monthly_veg.nc'
#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
out_tag = 'SCsim_WL_JULES-WFDEI-Zinke-hydro1k_Static_LAI'
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
dz_soil = 0.1   # top soil layer is 10cm  
#
# List of Parameters required from JULES output:

# 2.4 Read in T_soil, SM, SC, LAI, CanHt and PFT from JULES output
print 'Reading in JULES output data'
inf    = nc.Dataset(JULES_file,'r')
lat    = inf.variables['latitude'][:].squeeze()
lon    = inf.variables['longitude'][:].squeeze()
LAI    = inf.variables['lai'][:].squeeze()
LST    = inf.variables['tstar'][:,:5,:,:].squeeze()
SM     = inf.variables['smcl'][:,0,:,:].squeeze()
T_soil = inf.variables['t_soil'][:,0,:,:].squeeze()
inf.close()

# Convert SM from kg of water in top soil layer to
#      volumtric water content
SM = SM * (1/dz_soil) * 1e-3

# 2.5 read in pft fractions ancil file
inf=nc.Dataset(WFDEI_frac_file,'r')
PFTfrac = inf.variables[frac_varname][:5,:]
inf.close()

# places subsection into correct location
inf=nc.Dataset(WFDEI_soil_file,'r')
SMwilt  = inf.variables[smwilt_varname][:]
SMsat = inf.variables[smsat_varname][:]
inf.close()

# Ice points for masking
SMsat   = np.ma.masked_equal(SMsat,0)
ICEmask = np.ones_like(SMsat)
SMwilt  = SMwilt*ICEmask

SMwilt_sthf = SMwilt/SMsat
SM_sthf     = (SM/SMsat)

# Mask all other arrays to match SM
LAI    = LAI*ICEmask
LST    = LST*ICEmask
T_soil = T_soil*ICEmask


# Calculate dimension lengths
nMONTHs     = nMONTHs   # set at top
nYEARs      = int(LAI.shape[0]/nMONTHs)
nPFTs       = LAI.shape[1]
nlandpoints = LAI.shape[2]

######################################################
# 2.9 Read in grid file for plotting and diagnostics
####################################################
inf=nc.Dataset(WFDEI_gridfile,'r')
grindex=inf.variables['land_index'][:]
grifrac=inf.variables['land_fraction'][:]
grimask=np.ones_like(grindex)
inf.close()

#####################################################################################
# 3. Calculate/Create monthly climatologies
######################################
# 3.1 Average
# PFT params
LAI_climat     = np.mean(LAI.reshape(nYEARs,nMONTHs,nPFTs,nlandpoints),axis=0)
LAI_mean       = LAI_climat[0,:,:]
LST_climat     = np.mean(LST.reshape(nYEARs,nMONTHs,nPFTs,nlandpoints),axis=0)
LST_mean       = np.mean(LST,axis=0)
PFTfrac_climat = np.ma.array([PFTfrac for i in range(nMONTHs)])
# Soil Params
SM_sthf_climat = np.mean(SM_sthf.reshape(nYEARs,nMONTHs,nlandpoints),axis=0)
T_soil_climat  = np.mean(T_soil.reshape(nYEARs,nMONTHs,nlandpoints),axis=0)
SMwilt_climat  = np.ma.array([SMwilt_sthf for i in range(nMONTHs)])

# 3.2 transpose PFT dim to position 0
LAI_climat     = LAI_climat.transpose(1,0,2)
LST_climat     = LST_climat.transpose(1,0,2)
PFTfrac_climat = PFTfrac_climat.transpose(1,0,2)

# 3.3 calculate veg fraction (sum of all PFTs)
VEGfrac_climat = np.sum(PFTfrac_climat,axis=0)

#####################################################################################
# 4. Litterfall component
##########################
# 
print 'Litterfall Calculations'
# 4.1 Calculate local litterfall
LocLit_C_ECP, Cveg_ECP = J_SC_F.JULES_LOCAL_LITTERFALL( LAI_mean, LST_mean,  \
                                                        get_Cv=True,with_phenol=False)
LocLit_C_ECP = LocLit_C_ECP*ICEmask
Cveg_ECP     = Cveg_ECP*ICEmask
#
# 4.2 Calculate the total litterfall per PFT
Lit_C_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, LocLit_C_ECP, Cveg_ECP, \
                                           Per_PFT=True )
Lit_C_ECP=Lit_C_ECP*ICEmask
#
# 4.3 Calculate the total litterfall for all PFTs
Lit_C_total_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, LocLit_C_ECP, Cveg_ECP )
Lit_C_total_ECP = Lit_C_total_ECP*ICEmask
#
# 4.4 Split the total litter fall into DPM and RPM components
Lit_SCpools_ECP = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_ECP)
Lit_SCpools_ECP = Lit_SCpools_ECP*ICEmask

#####################################################################################
# 5. Soil Respiration Factor component
#######################################
print 'Soil Respiration Calculations'
# As solving for Soil Carbon we only want the factor of the Soil 
#  respiration function at this point. Simple hack is to call my soil respiration
#  function with Soil Carbon = 1
#  Soil_Resp_Fact= kappa * F(Ts) * F(SM) * F(Vf)
# 5.1 Create Dummy soil carbon array
Carbon_Ones=np.ones([4,nMONTHs,nlandpoints])
#
# The input data all have 4 soil layers where the functions only require the surface layer,
# hence we zero index the first dimension of the soil parameters
# 5.2a Soil_Resp_factor for Q10 temperature function
print 'Soil_Resp_Fact_Q10t'
Soil_Resp_Fact_Q10t = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_climat,     \
                                                     SM_sthf_climat, SMwilt_climat,  \
                                                     VEGfrac_climat, Tfunc='Q10',    \
                                                     OUTPUT_opt='pools'  )
# multiply by ICEmask
Soil_Resp_Fact_Q10t = Soil_Resp_Fact_Q10t * ICEmask
#
# 5.2b Soil_Resp_factor for RothC temperature function
print 'Soil_Resp_Fact_RothCt'
Soil_Resp_Fact_RothCt = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones, T_soil_climat,    \
                                                       SM_sthf_climat, SMwilt_climat, \
                                                       VEGfrac_climat, Tfunc='RothC', \
                                                       OUTPUT_opt='pools')
# multiply by ICEmask
Soil_Resp_Fact_RothCt = Soil_Resp_Fact_RothCt * ICEmask
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
                                   T_soil_climat,SM_sthf_climat,SMwilt_climat, \
                                   VEGfrac_climat, Tfunc='Q10',OUTPUT_opt='pools',  \
                                   kappa=kappa_s_1SC )
Soil_Resp_Fact_singlepool_Q10t = Soil_Resp_Fact_singlepool_Q10t.squeeze() * ICEmask
# 5.3b RothC temperature Function:
print 'Soil_Resp_Fact_singlepool_RothCt'
Soil_Resp_Fact_singlepool_RothCt = \
    J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones[0,:].reshape(temp_CarbonShape), \
                                   T_soil_climat,SM_sthf_climat,SMwilt_climat, \
                                   VEGfrac_climat, Tfunc='RothC',OUTPUT_opt='pools', \
                                    kappa=kappa_s_1SC )
Soil_Resp_Fact_singlepool_RothCt = Soil_Resp_Fact_singlepool_RothCt.squeeze() * ICEmask
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
C_dpm_Q10t   = (Lit_SCpools_ECP[0,:]       /  \
                np.sum(Soil_Resp_Fact_Q10t[0,:,:], axis=0) ) /  \
                (3600.*24.*360.)

C_dpm_RothCt = (Lit_SCpools_ECP[0,:]       /  \
                np.sum(Soil_Resp_Fact_RothCt[0,:,:], axis=0) ) /  \
                (3600.*24.*360.)

#
# 6.2 Resistant Plant Material:
#  C_rpm = (f_dpm*Lit_c)/ [Soil_resp_factor]_dpm
#      Calculate Cpool based on Q10 temp function and RothC temp function
C_rpm_Q10t   = (Lit_SCpools_ECP[1,:]       /  \
                np.sum(Soil_Resp_Fact_Q10t[1,:,:], axis=0) ) /  \
                (3600.*24.*360.)

C_rpm_RothCt = (Lit_SCpools_ECP[1,:]       /  \
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

######################################################################################
# 9. Output the data to netCDF files
#########################################################################
#
print 'Output data to netCDF'
print OUT_DIR+out_tag+'_FastSoilCarbon.nc'
# 9.1 open output file for 4 pools data.
outf=nc.Dataset(OUT_DIR+out_tag+'_FastSoilCarbon.nc','w')

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
outvar=outf.createVariable('LST','float32',('pft','x'))
outvar.units='K'
outvar.note='LST for the 5 standard JULES PFTs'
outvar[:]=LST_mean

# Monthly T_soil
outvar=outf.createVariable('T_SOIL','float32',('month','x'))
outvar.units='K'
outvar.note='Monthly Soil Temperature (Top layer 10cm)'
outvar[:]=T_soil_climat

# Monthly Soil Moisture
outvar=outf.createVariable('SM','float32',('month','x'))
outvar.units='-'
outvar.note='Monthly Soil Moisture as fraction of saturation (Top layer 10cm)'
outvar[:]=SM_sthf_climat

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
outvar=outf.createVariable('Local_Litterfall_C','float32',('pft','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Local Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=LocLit_C_ECP

# Total Litterfall per PFT
outvar=outf.createVariable('Litterfall_C','float32',('pft','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=Lit_C_ECP
    
# Total Litterfall
outvar=outf.createVariable('Total_Litterfall_C','float32',('x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon summed over all 5 standard JULES PFTs'
outvar[:]=Lit_C_total_ECP

# Lit_SCpools
outvar=outf.createVariable('Lit_SCpools','float32',('LitSCpool','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon that goes into each Carbon Pool'
outvar[:]= Lit_SCpools_ECP
    
# Cveg per PFT
outvar=outf.createVariable('c_veg','float32',('pft','x'))
outvar.units='kgC m-2'
outvar.note='Total vegetation Carbon for each of the 5 standard JULES PFTs'
outvar[:]=Cveg_ECP
    
# PFT fraction
outvar=outf.createVariable('PFT_frac','float32',('pft','x'))
outvar.units='kgC m-2'
outvar.note='Fractional cover of the 5 standard JULES PFTs'
outvar[:]=PFTfrac
    
# LAI
outvar=outf.createVariable('LAI','float32',('pft','x'))
outvar.units='m2 / m2'
outvar.note='LAI of the 5 standard JULES PFTs'
outvar[:]=LAI_mean
    
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

outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.close()



