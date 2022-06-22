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
chess_BCsoil_file='/users/eow/edwcom/CHESS/'


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
    gamma_v= [0.010,0.004,0.10,0.05,0.10] # [0.005,0.007,0.20,0.05,0.20]
    alpha_dr=[0.25,0.25,0.67,0.33,0.67]   # 1.44
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
sp=180
ep=240
S_level=1
# 4.1 Read in data from Spin 1 run
# Remove first year to avoid funny initial litterfall values
# squeeze and transpose soil/pft dim to position 0
inf      = nc.Dataset(spin1_outfile,'r')
time=nctime.num2date(inf.variables['time'][sp:ep],     \
                     units=inf.variables['time'].units,\
                     calendar='standard')
loclit_c = inf.variables['lit_c'][sp:ep,:].squeeze().transpose(1,0,2)
c_veg_X  = inf.variables['c_veg'][sp-1:ep,:].squeeze().transpose(1,0,2) # Reading additional month 
                                                                      #  for alternate lit_c calc
resp_s   = inf.variables['resp_s'][sp:ep,:].squeeze().transpose(1,0,2)
resp_s_diag   = inf.variables['resp_s_diag'][sp:ep,:].squeeze().transpose(1,0,2)
cs       = inf.variables['cs'][sp:ep,:].squeeze().transpose(1,0,2)
frac     = inf.variables['frac'][sp:ep,:].squeeze().transpose(1,0,2)
#lai      = inf.variables['lai'][sp:ep,:].squeeze().transpose(1,0,2)
#lai_bal  = inf.variables['lai_bal'][sp:ep,:].squeeze().transpose(1,0,2)
#tstar    = inf.variables['tstar'][sp:ep,:].squeeze().transpose(1,0,2)
#npp      = inf.variables['npp'][sp:ep,:].squeeze().transpose(1,0,2)
t_soil    = inf.variables['t_soil'][sp:ep,S_level,:].squeeze()
sthu      = inf.variables['sthu'][sp:ep,S_level,:].squeeze()
inf.close()
nTSTEPs=len(time)

inf_d = nc.Dataset(spin1_dumpfile,'r')
SM_wilt = np.array([ inf_d.variables['sm_wilt'][S_level,:]/ \
                     inf_d.variables['sm_sat'][S_level,:]   \
                     for i in range(nTSTEPs) ] )
inf_d.close()

# Remove early month for normal calculations
c_veg=c_veg_X[:,:-1,:]

# Calculate the change in cveg from month to month
#delta_cveg=c_veg_X[:,1:,:]-c_veg_X[:,:-1,:]

# Calculate lit_c based on npp not used in change in cveg
#  Have to convert npp to kgC month^-1
# Then convert litter to kgC per year
#lit_c_npp = ((npp*30.*24.*3600.)-delta_cveg)*12.

# Extract PFT frac and then calc Veg frac
PFTfrac = frac[:nPFTs,:]
VEGfrac = np.sum(PFTfrac,axis=0)
if 'CHESS' in run:
    ICEmask = np.where(VEGfrac[0,:]<-1)[0]
else:
    ICEmask = np.where(frac[8,0,:]>0.1)[0]

# 4.2 Calculate local litterfall
#print 'Litterfall Calculations'
#Lit_C = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, loclit_c, c_veg, \
#                                           Per_PFT=True, 
#                                           gamma_v=gamma_v )

# JULES lit_c output has taken distubance and competition terms into account
# 4.3 Split the total litter fall into DPM and RPM components
#Lit_SCpools = J_SC_F.JULES_LITTER_to_SCpool(Lit_C,alpha_dr=alpha_dr)
Lit_SCpools = J_SC_F.JULES_LITTER_to_SCpool(loclit_c*PFTfrac,alpha_dr=alpha_dr)

# 4.4 Soil Respiration Factor component
Soil_Resp_Fact = resp_s/(cs*(1-beta_r))
#Soil_Resp_Fact = resp_s_diag/(cs*(1-beta_r)*3600.*24.*365.)

#Cs_ones=np.ones_like(cs)
#Soil_Resp_Fact_2=J_SC_F.JULES_SOIL_RESPIRATION( Cs_ones, t_soil, sthu, SM_wilt, VEGfrac, \
#                                         kappa=[3.22e-7,9.65e-9,2.12e-8,6.43e-10], \
#                                         Tfunc='Q10', Q10=2.0, \
#                                         OUTPUT_opt='pools'  )
#
#Soil_Resp_Fact=Soil_Resp_Fact_2


# 4.5 Calculate the soil carbon for each pool using the albmar/ECP equations
# DPM
C_dpm  = (np.sum(Lit_SCpools[0,:,:], axis=0)  /  \
              np.sum(Soil_Resp_Fact[0,:,:], axis=0) ) /  \
              (3600.*24.*365.)
#C_dpm  = np.mean(cs[0,:,:],axis=0)
# RPM
C_rpm  = (np.sum(Lit_SCpools[1,:,:], axis=0)       /  \
              np.sum(Soil_Resp_Fact[1,:,:], axis=0) ) /  \
              (3600.*24.*365.)
#C_rpm  = np.mean(cs[1,:,:],axis=0)

    # BIO and HUM factors:
dpm_rpm_term   = (kappa_s_RothC[0]*C_dpm)+(kappa_s_RothC[1]*C_rpm)

bio_factor =  (0.46*beta_r) /( 1.65*(1-beta_r)*kappa_s_RothC[2] )
hum_factor =  (0.54*beta_r) /( 1.65*(1-beta_r)*kappa_s_RothC[3] )

# BIO
C_bio   = bio_factor*dpm_rpm_term
# HUM
C_hum   = hum_factor*dpm_rpm_term
#
# Mask out ICE points
cs_2 = np.array( [C_dpm,C_rpm,C_bio,C_hum] )
#cs_2 = np.array( [np.mean(C_dpm,axis=0),np.mean(C_rpm,axis=0),\
#                  np.mean(C_bio,axis=0),np.mean(C_hum,axis=0)] )
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
diff_percent=( (cs_2-cs_1)/ ((cs_2+cs_1)/2.) ) *100.

CS_maxes = [ 0.2, 8.0, 0.5, 25.0 ]
CS_ranges=[ [-0.05,0.05],[-2,2],[-0.2,0.2],[-10,10] ]
CS_perc_ranges=[ [-10,10],[-10,10],[-10,10],[-10,10] ]
CS_perc_ranges=[ [-50,50],[-50,50],[-50,50],[-50,50] ]

for iPOOL in range(nPOOLs):
    plt.figure(figsize=[18,8])
    plt.subplot(1,4,1)
    plt.imshow(cs_2[iPOOL,grindex]*grimask,origin='bottom', \
               vmax=CS_maxes[iPOOL]) 
    plt.title('JULES-sim')
    plt.colorbar(fraction=0.05)
    plt.subplot(1,4,2)
    plt.imshow(cs_1[iPOOL,grindex]*grimask,origin='bottom', \
               vmax=CS_maxes[iPOOL])
    plt.title('JULES-Kov')
    plt.colorbar(fraction=0.05)
    plt.subplot(1,4,3)
    plt.imshow((cs_2-cs_1)[iPOOL,grindex]*grimask,origin='bottom', \
                vmax=CS_ranges[iPOOL][1], vmin=CS_ranges[iPOOL][0], \
                cmap='RdBu_r' )
    plt.title('JULES-sim - JULES-Kov')
    plt.colorbar(fraction=0.05)
    plt.subplot(1,4,4)
    plt.imshow(diff_percent[iPOOL,grindex]*grimask,origin='bottom', \
                vmax=CS_perc_ranges[iPOOL][1], vmin=CS_perc_ranges[iPOOL][0], \
                cmap='RdBu_r' )
    plt.title('JULES-sim - JULES-Kov (%)')
    plt.colorbar(fraction=0.05)
    plt.suptitle('C '+POOL_names[iPOOL]+' Differences',fontsize=30)
    plt.savefig('test_plots/FUDGE_CS_Diffs_'+POOL_names[iPOOL]+'.png', bbox_inches='tight',\
                pad_inches=0.2)
    plt.close()
        

plotdata=np.mean(loclit_c*PFTfrac,axis=1)
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


plotdata=np.mean(Soil_Resp_Fact/Soil_Resp_Fact_2,axis=1)#*3600.*24*365.
title='SoilRespFactor_ratio_By_Pool'
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

plotdata=cs_2/cs_1
title='CS_ratio'
plt.figure(figsize=[18,10])
for iPOOL in range(nPOOLs):
    print iPOOL
    plt.subplot(2,2,iPOOL+1)
    plt.imshow(plotdata[iPOOL,grindex]*grimask,origin='bottom',vmin=0,vmax=2)
    plt.title(POOL_names[iPOOL])
    plt.colorbar()

plt.suptitle(title,fontsize=30)
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
plt.suptitle(title,fontsize=30)
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
plt.suptitle(title,fontsize=30)
plt.savefig('test_plots/'+title+'.png')
plt.show()

