#!/usr/bin/python2.7
#
############################################################
# 
# Program: JULES_with_SCequi.py
# Author: E. Comyn-Platt, edwcom@ceh.ac.uk
# Purpose: To perform a multi-stage spin up of JULES
#          incorporating SC simulator method
###########################################################

import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
import JULES_SC_functions as J_SC_F
import netCDF4 as nc
import netcdftime as nctime


##############################################################################
# 0. Read in parsed arguments
##################################
#
run = 'CHESS_Koven' 
out_dir = '/prj/GREENHOUSE/SC_simulator/JULES_output/CHESS_KovenMethod/'
spin1_outfile = out_dir+'JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod_spin2.monthly_mean.nc' 
spin1_dumpfile = out_dir+'JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod_spin2.dump.19610101.0.nc'
spin2_startdump_file = out_dir+'JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod_spin2.Equilibrium_FromBasics.nc'
gridfile='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'


##############################################################################
# 1. Define directories and files etc.
##################################
#  FILENAMES NOW PARSED IN
print 'Performing SCequilibrium.py'
print 'spin1_outfile = '+spin1_outfile
print 'spin1_dumpfile = '+spin1_dumpfile
print 'spin2_startdump_file = '+spin2_startdump_file

##############################################################################
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
print 'beta_r = ',beta_r
#
dz_soil = 0.1   # top soil layer is 10cm  
#
nPFTs=5
PFT_names=['BL','NL','C3','shrub','crop']
nPOOLs=4
POOL_names=['DPM','RPM','BIO','HUM']
if 'CHESS' in run:
    gamma_v=[0.005,0.007,0.20,0.05,0.20]
    alpha_dr=[0.25,0.25,0.67,0.33,0.67]
elif 'WFDEI' in run:
    gamma_v=[0.010,0.004,0.10,0.10,0.05]
    alpha_dr=[0.25,0.25,0.67,0.67,0.33]
else:
    print 'Unknown run, setting gamma_v to default [0.010,0.004,0.10,0.10,0.05]'
    gamma_v=[0.010,0.004,0.10,0.10,0.05]
    alpha_dr=[0.25,0.25,0.67,0.67,0.33]


##############################################################################
# 3. Read Grid File for plotting
grinf=nc.Dataset(gridfile,'r')
grindex=grinf.variables['index_2D'][:]
lats_2d=grinf.variables['lats_2D'][:]
lons_2d=grinf.variables['lons_2D'][:]
grinf.close()
grimask=np.ones_like(grindex)


##############################################################################
# 4. Calculate equi CS
#     Using the Eleanor/Alberto/Eddy method calculate the equilibrium
#       soil carbon based on the soil resp and litter functions
#     
#     Output our first estimate of the equilibrium CS:
#       CS_2 = [CS_dpm_2,CS_rpm_2,CS_bio_2,CS_hum_2]
#
#     RPM and DPM should have reasonable equilibrium estimates now
#
##################################
#os.system('pwd')
print 'Reading Data'
# 4.1 Read in data from Spin 1 run
# Remove first year to avoid funny initial litterfall values
# squeeze and transpose soil/pft dim to position 0
inf      = nc.Dataset(spin1_outfile,'r')
time=nctime.num2date(inf.variables['time'][11:23],     \
                     units=inf.variables['time'].units,\
                     calendar='standard')
loclit_c = inf.variables['lit_c'][11:23,:].squeeze().transpose(1,0,2)
c_veg_X  = inf.variables['c_veg'][10:23,:].squeeze().transpose(1,0,2) # Reading additional month 
                                                                      #  for alternate lit_c calc
resp_s   = inf.variables['resp_s'][11:23,:].squeeze().transpose(1,0,2)
resp_s_diag   = inf.variables['resp_s_diag'][-12:,:].squeeze().transpose(1,0,2)
cs       = inf.variables['cs'][11:23,:].squeeze().transpose(1,0,2)
frac     = inf.variables['frac'][11:23,:].squeeze().transpose(1,0,2)
lai      = inf.variables['lai'][11:23,:].squeeze().transpose(1,0,2)
lai_bal  = inf.variables['lai_bal'][11:23,:].squeeze().transpose(1,0,2)
tstar    = inf.variables['tstar'][11:23,:].squeeze().transpose(1,0,2)
npp      = inf.variables['npp'][11:23,:].squeeze().transpose(1,0,2)
inf.close()

# Remove early month for normal calculations
c_veg=c_veg_X[:,:-1,:]

# Calculate the change in cveg from month to month
delta_cveg=c_veg_X[:,1:,:]-c_veg_X[:,:-1,:]

# Calculate lit_c based on npp not used in change in cveg
#  Have to convert npp to kgC month^-1
# Then convert litter to kgC per year
lit_c_npp = ((npp*30.*24.*3600.)-delta_cveg)*12.

# Extract PFT frac and then calc Veg frac
PFTfrac = frac[:nPFTs,:]
VEGfrac = np.sum(PFTfrac,axis=0)
if 'CHESS' in run:
    ICEmask = np.where(VEGfrac[0,:]<-1)[0]
else:
    ICEmask = np.where(frac[8,0,:]>0.1)[0]

#for imonth in range(6):
#    print imonth*2
#    plt.subplot(2,3,imonth+1)
#    plt.imshow(tstar[2,imonth*2,grindex]*grimask,origin='bottom')
#    plt.colorbar()
#plt.show()

#for iPFT in range(5):
#    print imonth*2
#    plt.subplot(2,3,iPFT+1)
#    plt.imshow(tstar[iPFT,2,grindex]*grimask,origin='bottom')
#    plt.colorbar()
#plt.show()


#for imonth in range(6):
#    print imonth*2
#    plt.subplot(2,3,imonth+1)
#    plt.imshow((loclit_c*PFTfrac)[3,imonth*2,grindex]*grimask,origin='bottom',vmin=0,vmax=1.5)
#    plt.colorbar()
#
#plt.show()

#print ICEmask
# 4.1 Calculate local litterfall
#      with modified PFTs [ BL, NL, C3, shrub, crop (C3) ]
LocLit_C_ECP, Cveg_ECP = \
        J_SC_F.JULES_LOCAL_LITTERFALL( lai, tstar[:5,:,:], get_Cv=True,           \
                                       with_phenol=True, LAI_b=lai_bal,           \
                                       gamma_0=[0.25,0.25,0.25,0.25,0.25],        \
                                       dT=[9.0,9.0,0.0,9.0,0.0],                  \
                                       Toff=[278.15,243.15,258.15,243.15,258.15], \
                                       a_wl=[0.65,0.65,0.005,0.10,0.005],         \
                                       b_wl=[1.667,1.667,1.667,1.667,1.667],      \
                                       sigma_l=[0.0375,0.1,0.025,0.05,0.025],     \
                                       gamma_r=[0.25,0.25,0.25,0.25,0.25],        \
                                       gamma_w=[0.01,0.01,0.20,0.05,0.20],        \
                                       gamma_p=[20.0,20.0,20.0,20.0,20.0],        \
                                       a_ws=[10.0,10.0,1.0,10.0,1.0],             \
                                       nu_sl=[0.01,0.01,0.01,0.01,0.01],          \
                                       )
#
#LocLit_C_EM, Cveg_EM = \
#         J_SC_F.JULES_LOCAL_LITTERFALL( lai, em_tstar_1D[:5,:,:], get_Cv=True,           \
#                                        with_phenol=True, LAI_b=lai_bal,           \
#                                        gamma_0=[0.25,0.25,0.25,0.25,0.25],        \
#                                        dT=[9.0,9.0,0.0,9.0,0.0],                  \
#                                        Toff=[278.15,243.15,258.15,243.15,258.15], \
#                                        a_wl=[0.65,0.65,0.005,0.10,0.005],         \
#                                        b_wl=[1.667,1.667,1.667,1.667,1.667],      \
#                                        sigma_l=[0.0375,0.1,0.025,0.05,0.025],     \
#                                        gamma_r=[0.25,0.25,0.25,0.25,0.25],        \
#                                        gamma_w=[0.01,0.01,0.20,0.05,0.20],        \
#                                        gamma_p=[20.0,20.0,20.0,20.0,20.0],        \
#                                        a_ws=[10.0,10.0,1.0,10.0,1.0],             \
#                                        nu_sl=[0.01,0.01,0.01,0.01,0.01],          \
#                                        )


# 4.2 Calculate local litterfall
print 'Litterfall Calculations'
Lit_C = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, loclit_c, c_veg, \
                                           Per_PFT=True, 
                                           gamma_v=gamma_v )


Lit_C_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, LocLit_C_ECP, Cveg_ECP, \
                                           Per_PFT=True, 
                                           gamma_v=gamma_v )


#
# 4.3 Split the total litter fall into DPM and RPM components
Lit_SCpools = J_SC_F.JULES_LITTER_to_SCpool(Lit_C,alpha_dr=alpha_dr)
#Lit_SCpools_ECP = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_ECP,alpha_dr=alpha_dr)

# 4.4 Soil Respiration Factor component
#Soil_Resp_Fact = np.zeros_like(resp_s)
#for iPOOL in range(nPOOLs):
#    Soil_Resp_Fact[iPOOL,:] = (resp_s[iPOOL,:])/(cs[iPOOL,:])  #*kappa_s_RothC[iPOOL])

Soil_Resp_Fact = resp_s/(cs*(1-beta_r))
#Soil_Resp_Fact = resp_s_diag/(cs*(1-beta_r))


#Cs_ones=np.ones_like(cs)
#Soil_Resp_Fact_2=JULES_SOIL_RESPIRATION( Cs_ones, Ts, SM, SMw, Vf, \
#                       kappa=[3.22e-7,9.65e-9,2.12e-8,6.43e-10], \
#                       Tfunc='Q10', Q10=2.0, \
#                       OUTPUT_opt='pools'  )


# 4.5 Calculate the soil carbon for each pool using the albmar/ECP equations
# DPM
C_dpm  = (np.sum(Lit_SCpools[0,:,:], axis=0)  /  \
              np.sum(Soil_Resp_Fact[0,:,:], axis=0) ) /  \
              (3600.*24.*360.)
# RPM
C_rpm  = (np.sum(Lit_SCpools[1,:,:], axis=0)       /  \
              np.sum(Soil_Resp_Fact[1,:,:], axis=0) ) /  \
              (3600.*24.*360.)

#C_dpm_ECP  = (np.sum(Lit_SCpools_ECP[0,:,:], axis=0)  /  \
#              np.sum(Soil_Resp_Fact[0,:,:], axis=0) ) /  \
#              (3600.*24.*360.)
## RPM
#C_rpm_ECP  = (np.sum(Lit_SCpools_ECP[1,:,:], axis=0)       /  \
#              np.sum(Soil_Resp_Fact[1,:,:], axis=0) ) /  \
#              (3600.*24.*360.)

C_rpm_K = np.mean(cs[1,:,:],axis=0)


    # BIO and HUM factors:
dpm_rpm_term   = (kappa_s_RothC[0]*C_dpm)+(kappa_s_RothC[1]*C_rpm)

bio_factor =  (0.46*beta_r) /( (1-beta_r)*kappa_s_RothC[2] )
hum_factor =  (0.54*beta_r) /( (1-beta_r)*kappa_s_RothC[3] )

# BIO
C_bio   = bio_factor*dpm_rpm_term
# HUM
C_hum   = hum_factor*dpm_rpm_term
#
# Mask out ICE points
cs_2 = np.array( [C_dpm,C_rpm,C_bio,C_hum] )
cs_2[:,ICEmask]=0.0

# set values less than or equal to zero to 1e-6
cs_2[cs_2<=0]=1e-6

# write new cs to spin2_startdump_file
print "Writing SC equilibrium output to: ", spin2_startdump_file
inf = nc.Dataset(spin1_dumpfile,'r')
outf= nc.Dataset(spin2_startdump_file,'w')
cs_1=inf.variables['cs'][:].squeeze()
for dim in inf.dimensions:
    outf.createDimension( str(dim),len(inf.dimensions[dim]) )
        
for var in inf.variables:
    outvar = outf.createVariable( str(var), \
                                  inf.variables[var].dtype, \
                                  inf.variables[var].dimensions )
    if str(var)=='cs':
        outvar[:]=cs_2
    else:
        outvar[:]=inf.variables[str(var)][:]

inf.close()
outf.close()



quit()



#Lit_C=loclit_c*PFTfrac
data1=np.mean(Lit_C_ECP*12,axis=1)*frac[:5,0,:]
data1name='Lit_C_ECP'
data1=np.mean(Lit_C_ECP_pheno,axis=1)*frac[:5,0,:]
data1name='Lit_C_ECP_pheno'
data2=np.mean(loclit_c,axis=1)*frac[:5,0,:]
data2name='loclit_c'
data2=np.mean(Lit_C,axis=1)*frac[:5,0,:]
data2name='Lit_C'

PFT_names=['BL','NL','C3','shrub','crop']

for iPFT in range(len(PFT_names)):
    plt.figure(figsize=[15,10])
    pft_data1=data1[iPFT,:]
    pft_data2=data2[iPFT,:]
    
    plt.subplot(1,3,1)
    plt.imshow(pft_data1[grindex]*grimask,origin='bottom')
    plt.title(data1name)
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,2)
    plt.imshow(pft_data2[grindex]*grimask,origin='bottom')
    plt.title(data2name)
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,3)
    plt.imshow((pft_data1-pft_data2)[grindex]*grimask,origin='bottom')
    plt.title(data1name+' - '+data2name)
    plt.colorbar(fraction=0.05)
    plt.suptitle(PFT_names[iPFT]+' Differences')
    plt.savefig('test_plots/Diff_'+data1name+'_'+data2name+'_'+PFT_names[iPFT]+'.png')
    plt.close()


cs_1 = np.mean(cs,axis=1)

for iPOOL in range(nPOOLs):
    plt.figure(figsize=[15,8])
    plt.subplot(1,3,1)
    plt.imshow(cs_2[iPOOL,grindex]*grimask,origin='bottom')
    plt.title('JULES-sim')
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,2)
    plt.imshow(cs_1[iPOOL,grindex]*grimask,origin='bottom')
    plt.title('JULES-Kov')
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,3)
    plt.imshow((cs_2-cs_1)[iPOOL,grindex]*grimask,origin='bottom')
    plt.title('JULES-sim - JULES-Kov')
    plt.colorbar(fraction=0.05)
    plt.suptitle('C '+POOL_names[iPOOL]+' Differences')
    plt.savefig('test_plots/CS_Diffs_'+POOL_names[iPOOL]+'.png')
    plt.close()
    


data1=np.mean(resp_s,axis=1)*3600.*24*365
data1name='resp_s'
data2=np.mean(resp_s_diag,axis=1)
data2name='resp_s_diag'
title= 'RESP_S_DIFFERENCES'

for iPOOL in range(nPOOLs):
    plt.figure(figsize=[15,8])
    plt.subplot(1,3,1)
    plt.imshow(data1[iPOOL,grindex]*grimask,origin='bottom')
    plt.title(data1name)
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,2)
    plt.imshow(data2[iPOOL,grindex]*grimask,origin='bottom')
    plt.title(data2name)
    plt.colorbar(fraction=0.05)
    plt.subplot(1,3,3)
    plt.imshow((data1-data2)[iPOOL,grindex]*grimask,origin='bottom')
    plt.title(data1name+' - '+data2name)
    plt.colorbar(fraction=0.05)
    plt.suptitle(title)
    plt.savefig('test_plots/'+title+'.png')
    plt.close()
    






plt.figure(figsize=[18,10])
 
plt.subplot(2,3,1)
plt.imshow(C_rpm_ECP[grindex]*grimask,origin='bottom',vmax=7,vmin=0)
plt.title('ECP RPM')
plt.colorbar(fraction=0.05)
plt.subplot(2,3,2)
plt.imshow(C_rpm[grindex]*grimask,origin='bottom',vmax=7,vmin=0)
plt.title('JULES-sim RPM')
plt.colorbar(fraction=0.05)
plt.subplot(2,3,3)
plt.imshow(C_rpm_K[grindex]*grimask,origin='bottom',vmax=7,vmin=0)
plt.title('JULES-Kov RPM')
plt.colorbar(fraction=0.05)

plt.subplot(2,3,4)
plt.imshow((C_rpm_ECP-C_rpm)[grindex]*grimask,origin='bottom',vmax=1,vmin=-1)
plt.title('ECP - JULES-sim')
plt.colorbar(fraction=0.05)
plt.subplot(2,3,5)
plt.imshow((C_rpm_ECP-C_rpm_K)[grindex]*grimask,origin='bottom',vmax=1,vmin=-1)
plt.title('ECP - JULES-Kov')
plt.colorbar(fraction=0.05)
plt.subplot(2,3,6)
plt.imshow((C_rpm-C_rpm_K)[grindex]*grimask,origin='bottom',vmax=1,vmin=-1)
plt.title('JULES-sim - JULES-Kov')
plt.colorbar(fraction=0.05)
plt.suptitle('C RPM Differences')
plt.savefig('RPMDiffs.png')
plt.close()


for imonth in range(6):
    print imonth*2
    plt.subplot(2,3,imonth+1)
    plt.imshow([2,imonth*2,grindex]*grimask,origin='bottom',vmin=0)
    plt.title('month = '+str(imonth*2))
    plt.colorbar()

plt.show()

plotdata=np.mean(Lit_C,axis=1)
title='Litterfall_By_PFT'
plt.figure(figsize=[18,10])
for iPFT in range(nPFTs):
    print iPFT
    plt.subplot(2,3,iPFT+1)
    plt.imshow(plotdata[iPFT,grindex]*grimask,origin='bottom',vmin=0,vmax=1)
    plt.title(PFT_names[iPFT])
    plt.colorbar()
plt.subplot(2,3,6)
plt.imshow(np.sum(plotdata,axis=0)[grindex]*grimask,origin='bottom',vmin=0,vmax=2)
plt.title('Total')
plt.colorbar()
plt.suptitle(title)
plt.savefig('test_plots/'+title+'.png')
plt.show()


plotdata=np.mean(Soil_Resp_Fact,axis=1)*3600.*24*365.
title='SoilRespFactor_By_Pool'
plt.figure(figsize=[18,10])
for iPOOL in range(nPOOLs):
    print iPOOL
    plt.subplot(2,3,iPOOL+1)
    plt.imshow(plotdata[iPOOL,grindex]*grimask,origin='bottom',vmin=0)#,vmax=1)
    plt.title(POOL_names[iPOOL])
    plt.colorbar()
plt.subplot(2,3,6)
plt.imshow(np.sum(plotdata,axis=0)[grindex]*grimask,origin='bottom',vmin=0)#,vmax=2)
plt.title('Total')
plt.colorbar()
plt.suptitle(title)
plt.savefig('test_plots/'+title+'.png')
plt.show()


plotdata=np.mean(resp_s,axis=1)*3600.*24*365.
title='RespS_By_Pool'
plt.figure(figsize=[18,10])
for iPOOL in range(nPOOLs):
    print iPOOL
    plt.subplot(2,3,iPOOL+1)
    plt.imshow(plotdata[iPOOL,grindex]*grimask,origin='bottom',vmin=0)#,vmax=1)
    plt.title(POOL_names[iPOOL])
    plt.colorbar()
plt.subplot(2,3,6)
plt.imshow(np.sum(plotdata,axis=0)[grindex]*grimask,origin='bottom',vmin=0)#,vmax=2)
plt.title('Total')
plt.colorbar()
plt.suptitle(title)
plt.savefig('test_plots/'+title+'.png')
plt.show()

plotdata=np.mean(cs,axis=1)
title='Koven_CS_By_Pool'
plt.figure(figsize=[18,10])
for iPOOL in range(nPOOLs):
    print iPOOL
    plt.subplot(2,3,iPOOL+1)
    plt.imshow(plotdata[iPOOL,grindex]*grimask,origin='bottom',vmin=0)#,vmax=1)
    plt.title(POOL_names[iPOOL])
    plt.colorbar()
plt.subplot(2,3,6)
plt.imshow(np.sum(plotdata,axis=0)[grindex]*grimask,origin='bottom',vmin=0)#,vmax=2)
plt.title('Total')
plt.colorbar()
plt.suptitle(title)
plt.savefig('test_plots/'+title+'.png')
plt.show()


plotdata=cs_2
title='ECP_CS_By_Pool'
plt.figure(figsize=[18,10])
for iPOOL in range(nPOOLs):
    print iPOOL
    plt.subplot(2,3,iPOOL+1)
    plt.imshow(plotdata[iPOOL,grindex]*grimask,origin='bottom',vmin=0)#,vmax=1)
    plt.title(POOL_names[iPOOL])
    plt.colorbar()
plt.subplot(2,3,6)
plt.imshow(np.sum(plotdata,axis=0)[grindex]*grimask,origin='bottom',vmin=0)#,vmax=2)
plt.title('Total')
plt.colorbar()
plt.suptitle(title)
plt.savefig('test_plots/'+title+'.png')
plt.show()



plt.figure(figsize=[15,10])
total_data1=np.sum(data1,axis=0)
total_data2=np.sum(data2,axis=0)
plt.subplot(1,3,1)
plt.imshow(total_data1[grindex]*grimask,origin='bottom',vmax=1.5)
plt.title(data1name)
plt.colorbar(fraction=0.05)
plt.subplot(1,3,2)
plt.imshow(total_data2[grindex]*grimask,origin='bottom',vmax=1.5)
plt.title(data2name)
plt.colorbar(fraction=0.05)
plt.subplot(1,3,3)
plt.imshow((total_data1-total_data2)[grindex]*grimask,origin='bottom',vmin=-0.5,vmax=0.5)
plt.title(data1name+' - '+data2name)
plt.colorbar(fraction=0.05)
plt.suptitle('Total Differences')
plt.savefig('test_plots/Diff_'+data1name+'_'+data2name+'_Total.png')
plt.close()

plt.figure(figsize=[15,10])
plt.subplot(1,3,1)
plt.imshow(np.mean(np.sum(Leaf_LF,axis=0),axis=0)[grindex]*grimask,origin='bottom',vmin=0,vmax=2.5)
plt.title('Leaf_LF')
plt.colorbar(fraction=0.05)
plt.subplot(1,3,2)
plt.imshow(np.mean(np.sum(Root_LF,axis=0),axis=0)[grindex]*grimask,origin='bottom',vmin=0,vmax=2.5)
plt.title('Root_LF')
plt.colorbar(fraction=0.05)
plt.subplot(1,3,3)
plt.imshow(np.mean(np.sum(Stem_LF,axis=0),axis=0)[grindex]*grimask,origin='bottom',vmin=0,vmax=2.5)
plt.title('Stem LF')
plt.colorbar(fraction=0.05)
plt.suptitle('Total LF Components')
plt.savefig('test_plots/Total_LF_Components.png')
plt.close()
