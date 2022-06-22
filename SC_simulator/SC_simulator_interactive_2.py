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
nPOOLs=4
nPFTs=5
fill_value=-999
#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:

JULES_sources, JULES_sources_info = di_SC.jules_sources_info()

iJULES = di_SC.select_source(JULES_sources,Message='Select JULES source:')
#iJULES = 6
JULES_source=JULES_sources[iJULES]
JULES_source_info=JULES_sources_info[JULES_source]
print 'JULES file = ',JULES_source_info['file']

#LAI_source_info = di_SC.LAI_source_info(JULES_source_info)

#PFT_source_info = di_SC.PFT_source_info(JULES_source_info)

#SOIL_source_info = di_SC.SOIL_source_info(JULES_source_info)

#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
print 'Current outdir: '+OUT_DIR
temp = raw_input('Enter alternative OUT_DIR or hit return to continue: \n')

if temp!='':
    OUT_DIR=temp.copy()

out_tag = 'SCsim_J_litnresp_'+JULES_source_info['tag'] 

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
LocLit_C_J = inf.variables['lit_c'][:].squeeze()
Soil_Resp_J= inf.variables['resp_s'][:].squeeze()
CS_J       = inf.variables['cs'][:].squeeze()
frac_J  = inf.variables['frac'][:].squeeze()
C_veg_J    = inf.variables['c_veg'][:].squeeze()
inf.close()

PFTfrac_J = frac_J[:,:nPFTs,:]

ICEmask   = np.round(np.ma.masked_equal((frac_J[0,8,:]-1)*-1,0))

# Mask arrays and transpose to bring soil/PFT dimension to position 0
PFTfrac_J= PFTfrac_J*ICEmask
PFTfrac_J.fill_value=fill_value
PFTfrac_J=PFTfrac_J.transpose(1,0,2)

LocLit_C_J= LocLit_C_J*ICEmask
LocLit_C_J.fill_value=fill_value
LocLit_C_J=LocLit_C_J.transpose(1,0,2)

Soil_Resp_J= Soil_Resp_J*ICEmask
Soil_Resp_J.fill_value=fill_value
Soil_Resp_J=Soil_Resp_J.transpose(1,0,2)

CS_J=CS_J*ICEmask
CS_J.fill_value=fill_value
CS_J=CS_J.transpose(1,0,2)

C_veg_J=C_veg_J*ICEmask
C_veg_J.fill_value=fill_value
C_veg_J=C_veg_J.transpose(1,0,2)

Vegfrac = np.sum(PFTfrac_J,axis=0)

# Calculate dimension lengths
nMONTHs     = JULES_source_info['tsteps']   
nYEARs      = nMONTHs/12.
nlandpoints = LocLit_C_J.shape[2]

######################################################
# 2.9 Read in grid file for plotting and diagnostics
####################################################
gridfile = '/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
inf=nc.Dataset(gridfile,'r')
grindex=inf.variables['land_index'][:]-1
grifrac=inf.variables['land_fraction'][:]
grimask=np.ones_like(grindex)
inf.close()



#####################################################################################
# 4. Litterfall component
##########################
# 
print 'Litterfall Calculations'
# 4.1 Calculate local litterfall

# 4.2 Calculate the total litterfall per PFT
Lit_C_J = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_J, LocLit_C_J, C_veg_J, \
                                           Per_PFT=True )
Lit_C_J = Lit_C_J*ICEmask
Lit_C_J.fill_value=Lit_C_J.fill_value
#
# 4.3 Calculate the total litterfall for all PFTs
Lit_C_total_J = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac_J, LocLit_C_J, C_veg_J )
Lit_C_total_J = Lit_C_total_J*ICEmask
Lit_C_total_J.fill_value = fill_value
#
# 4.4 Split the total litter fall into DPM and RPM components
Lit_SCpools_J = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_J)
Lit_SCpools_J = Lit_SCpools_J*ICEmask
Lit_SCpools_J.fill_value=fill_value



#quit()
#####################################################################################
# 5. Soil Respiration Factor component
#######################################
print 'Soil Respiration Calculations'
# As solving for Soil Carbon we only want the factor of the Soil 
#  respiration function at this point. Simple hack is to call my soil respiration
#  function with Soil Carbon = 1
#  Soil_Resp_Fact= kappa * F(Ts) * F(SM) * F(Vf)

Soil_Resp_Fact_J = (Soil_Resp_J)/(CS_J*(1-beta_r))


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
C_dpm_Q10t   = (np.sum(Lit_SCpools_J[0,:,:], axis=0)       /  \
                np.sum(Soil_Resp_Fact_J[0,:,:], axis=0) ) /  \
                (3600.*24.*360.)


#
# 6.2 Resistant Plant Material:
#  C_rpm = (f_dpm*Lit_c)/ [Soil_resp_factor]_dpm
#      Calculate Cpool based on Q10 temp function and RothC temp function
C_rpm_Q10t   = (np.sum(Lit_SCpools_J[1,:,:], axis=0)       /  \
                       np.sum(Soil_Resp_Fact_J[1,:,:], axis=0) ) /  \
                       (3600.*24.*360.)

#
# See Eddy Equations:
# 6.3 Biomass:
#  C_bio =  ( 0.46 / ( (1/beta_r)-1 )*kappa_bio ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
# 6.4 Hummus:
#  C_hum =  ( 0.54 / ( (1/beta_r)-1 )*kappa_hum ) * [kappa_dpm*C_dpm + kappa_rpm*C_rpm]
#      Calculate Cpool based on Q10 temp function and RothC temp function
dpm_rpm_term_Q10t   = (kappa_s_RothC[0]*C_dpm_Q10t)+(kappa_s_RothC[1]*C_rpm_Q10t)
bio_factor =  0.46 / ( ( (1/beta_r)-1 )*kappa_s_RothC[2] )
hum_factor =  0.54 / (( (1/beta_r)-1 )*kappa_s_RothC[3] )
# 6.3:
C_bio_Q10t   = bio_factor*dpm_rpm_term_Q10t
# 6.4:
C_hum_Q10t   = hum_factor*dpm_rpm_term_Q10t
#
# 6.5 append to single array:
C_4pools_Q10t   = np.array( [C_dpm_Q10t,C_rpm_Q10t,C_bio_Q10t,C_hum_Q10t] )*ICEmask
#

C_4pools_Q10temp=np.ma.masked_equal(C_4pools_Q10t,0)

vmaxi = [1,20,5,50]
#plot_Cq10diags=False
if (plot_Cq10diags==True):
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.imshow(C_4pools_Q10t[i,grindex-1]*grimask,\
                   origin='bottom',cmap='YlOrBr')
        plt.colorbar()
    plt.show()
    

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


#  Soil Q10t Resp Factor
outvar=outf.createVariable('Soil_Resp_Fact_J','float32',('pool','month','x'))
outvar.units='kgC m-2 s-1'
outvar.note='Monthly Soil Respiration Factor using Q10 temperature function'
outvar[:]=Soil_Resp_Fact_J   

# Local Litterfall
outvar=outf.createVariable('Local_Litterfall_C','float32',('pft','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Local Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=LocLit_C_J  

# Total Litterfall per PFT
outvar=outf.createVariable('Litterfall_C','float32',('pft','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon for each of the 5 standard JULES PFTs'
outvar[:]=Lit_C_J  
    
# Total Litterfall
outvar=outf.createVariable('Total_Litterfall_C','float32',('month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon summed over all 5 standard JULES PFTs'
outvar[:]=Lit_C_total_J  

# Lit_SCpools
outvar=outf.createVariable('Lit_SCpools','float32',('LitSCpool','month','x'))
outvar.units='kgC m-2 year-1'
outvar.note='Total Litterfall Carbon that goes into each Carbon Pool'
outvar[:]= Lit_SCpools_J  
    
# Cveg per PFT
outvar=outf.createVariable('c_veg','float32',('pft','month','x'))
outvar.units='kgC m-2'
outvar.note='Total vegetation Carbon for each of the 5 standard JULES PFTs'
outvar[:]=C_veg_J  
    
# PFT fraction
outvar=outf.createVariable('PFT_frac','float32',('pft','month','x'))
outvar.units='kgC m-2'
outvar.note='Fractional cover of the 5 standard JULES PFTs'
outvar[:]=PFTfrac_J      
    
# Carbon 4 pools Q10t
outvar=outf.createVariable('C_4pools_Q10t','float32',('pool','x'))
outvar.units='kgC m-2'
outvar.note='Soil Carbon for the RothC pools calculated using the Q10 temperature function'
outvar[:]=C_4pools_Q10t
    

outf.author    = 'Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.JULES_RUN = JULES_source 
outf.note      = 'Calculations done using JULES litter and soil respiration output'
outf.close()



