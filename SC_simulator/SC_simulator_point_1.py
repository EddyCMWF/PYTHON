#!/usr/bin/env python

################################################################################
# 
# Program: SC_simulator_point.py
# 
# Python Script to simulate soil carbon estimate for a point location, 
#       i.e. flux site based on the long-term equilibrium condition, 
#       i.e. dCs/dt=0 .'. litterfall=microbial_soil_respiraiton       
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
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import JULES_SC_functions as J_SC_F
import netcdftime as nctime
#import data_info_SC as di_SC
#import plot_tools as PT
#
###########################################################
def select_source(sources,Message='Available Sources:'):
    print Message 
    for i in range(len(sources)):
        print str(i)+': '+sources[i]

    iSOURCE = input('Select an option:\n')

    return iSOURCE
#####################################################################################
nPFTs=5

JULES_dir = '/prj/jules_cosmos/jules_runs/outputs/'
SITES = [ 'tadhm', \
          'rdmer', \
          'crich', \
                ]
t_ress = ['day','tstep']

#JULES_files = [ 'tadhm_operative.tstep.nc', \
#                'rdmer_operative.tstep.nc', \
#                'crich_operative.tstep.nc', \
#                ]
JULES_init_files = [ 'tadhm_initial.dump.spin1.20141014.43200.nc', \
                     'rdmer_initial.dump.spin1.20150211.1800.nc', \
                     'crich_initial.dump.spin1.20141202.54000.nc', \
                    ]

iSITE = select_source(SITES,Message='Select a site: ')
iTRES = select_source(t_ress,Message='Select a time resolution: ')

SITE=SITES[iSITE]
t_res=t_ress[iTRES]

JULES_file=JULES_dir+SITE+'_operative.'+t_res+'.nc'
JULES_init_file=JULES_dir+JULES_init_files[iSITE]

OUT_DIR = '/users/eow/edwcom/SC_simulator/output/site_estimates/'
print 'Current outdir: '+OUT_DIR

if '.day.' in JULES_file:
    nSTEPs=365
elif '.tstep.' in JULES_file:
    nSTEPs=365*48


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

print 'Reading LST, T_soil and SM_sthu from JULES output data'
print 'Opening: '+JULES_file
inf        = nc.Dataset(JULES_file,'r')
J_time     = nctime.num2date( \
                              inf.variables['time'][:nSTEPs],\
                              units=inf.variables['time'].units,\
                              calendar=inf.variables['time'].calendar )
    
LST_gb     = inf.variables['tstar_gb'][:nSTEPs,:,:].squeeze()
T_soil     = inf.variables['t_soil'][:nSTEPs,0,:,:].squeeze()
SM_sthu    = inf.variables['sthu'][:nSTEPs,0,:,:].squeeze()
LAI        = inf.variables['lai'][:nSTEPs,:,:,:].squeeze().transpose(1,0)
Lit_C_J_total = inf.variables['lit_c_mean'][:nSTEPs,:].squeeze()
C_veg_J    = inf.variables['c_veg'][:nSTEPs,:].squeeze()
soil_resp_J= inf.variables['resp_s'][:nSTEPs,:].squeeze()
cs_J       = inf.variables['cs'][:nSTEPs,:].squeeze().transpose(1,0)
NPP_gb_J   = inf.variables['npp_gb'][:nSTEPs,:].squeeze()
pc_s_J     = inf.variables['pc_s'][:nSTEPs,:].squeeze()
inf.close()
LST=np.array([LST_gb for i in range(nPFTs)])

print 'Opening: '+JULES_init_file
inf=nc.Dataset(JULES_init_file,'r')
SMwilt_sthu= inf.variables['sm_wilt'][0,:].squeeze() / \
             inf.variables['sm_sat'][0,:].squeeze()
PFTfrac    = inf.variables['frac'][:nPFTs,:].squeeze()
PFTfrac[PFTfrac < 0.1] = 0
VEGfrac = np.sum(PFTfrac)
inf.close()

nTSTEPs = len(T_soil)

#LAI=np.mean(LAI,axis=1).reshape(nPFTs,1)
#LST=np.mean(LST,axis=1).reshape(nPFTs,1)
#SM_sthu=np.array([np.mean(SM_sthu)])
#T_soil=np.array([np.mean(T_soil)])
#PFTfrac=PFTfrac.reshape(nPFTs,1)
#SMwilt_sthu=np.array([np.mean(SMwilt_sthu)])
#VEGfrac=np.array([VEGfrac])


SMwilt_sthu= np.array([SMwilt_sthu for i in range(nTSTEPs)])
PFTfrac = np.array([PFTfrac for i in range(nTSTEPs)]).transpose(1,0)
VEGfrac =  np.array([VEGfrac for i in range(nTSTEPs)])
Carbon_Ones = np.ones([4,nTSTEPs])

LAI_b=np.zeros_like(LAI)+4
LocLit_C_ECP, Cveg_ECP = \
         J_SC_F.JULES_LOCAL_LITTERFALL( LAI, LST,  \
                                        get_Cv=True,with_phenol=False, \
                                        Toff=[278.15,233.15,278.15,278.15,233.15], \
                                        dT=[9.,9.,9.,9.,9.],      \
                                        gamma_p=[15.0,20.0,20.0,20.0,20.0],  \
                                        gamma_r=[0.25,0.15,0.25,0.25,0.25],  \
                                        gamma_w=[0.005,0.005,0.20,0.20,0.05] )

LocLit_C_ECP = \
         J_SC_F.JULES_LOCAL_LITTERFALL( LAI, LST,  \
                                        LAI_b=LAI_b, \
                                        with_phenol=False, \
                                        Toff=[278.15,233.15,278.15,278.15,233.15], \
                                        dT=[9.,9.,9.,9.,9.],      \
                                        gamma_p=[15.0,20.0,20.0,20.0,20.0],  \
                                        gamma_r=[0.25,0.15,0.25,0.25,0.25],  \
                                        gamma_w=[0.005,0.005,0.20,0.20,0.05] )

Lit_C_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, LocLit_C_ECP, Cveg_ECP, \
                                           Per_PFT=True,                    \
                                           gamma_v=[0.005,0.007,0.20,0.20,0.05]  )

Lit_C_total_ECP = J_SC_F.JULES_TOTAL_LITTERFALL( PFTfrac, LocLit_C_ECP, Cveg_ECP, \
                                           gamma_v=[0.010,0.004,0.10,0.10,0.05] )

Lit_SCpools_ECP = J_SC_F.JULES_LITTER_to_SCpool(Lit_C_ECP)


Soil_Resp_Fact = J_SC_F.JULES_SOIL_RESPIRATION( Carbon_Ones,\
                                                T_soil,     \
                                                SM_sthu, SMwilt_sthu,  \
                                                VEGfrac, Tfunc='Q10',    \
                                                OUTPUT_opt='pools'  )
#

C_dpm_Q10t   = (np.sum(Lit_SCpools_ECP[0,:], axis=0)       /  \
                np.sum(Soil_Resp_Fact[0,:], axis=0) ) /  \
                (3600.*24.*360.)

C_rpm_Q10t   = (np.sum(Lit_SCpools_ECP[1,:], axis=0)       /  \
                np.sum(Soil_Resp_Fact[1,:], axis=0) ) /  \
                (3600.*24.*360.)

dpm_rpm_term_Q10t   = (kappa_s_RothC[0]*C_dpm_Q10t)+(kappa_s_RothC[1]*C_rpm_Q10t)
bio_factor =  0.46 / ( ( (1/beta_r)-1 )*kappa_s_RothC[2] )
hum_factor =  0.54 / (( (1/beta_r)-1 )*kappa_s_RothC[3] )

C_bio_Q10t   = bio_factor*dpm_rpm_term_Q10t

C_hum_Q10t   = hum_factor*dpm_rpm_term_Q10t

#
# 6.5 append to single array:
C_4pools_Q10t   = np.array( [C_dpm_Q10t,C_rpm_Q10t,C_bio_Q10t,C_hum_Q10t] )

print C_4pools_Q10t 

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
#plt.plot(NPP_gb_J[1::48])
#plt.plot((C_veg_J[5::48,2]-C_veg_J[:-5:48,2])/(1800.*24))
#ax.plot(J_time[10:],(NPP_gb_J[10:])-((C_veg_J[10:,2]-C_veg_J[9:-1,2])/(1800.*24)))
ax.plot(J_time,Lit_C_J_total/(360*24*3600),color="blue",linewidth=2)
ax.plot(J_time,Lit_C_ECP[2,:]/(360*24*3600),color="cyan",linewidth=2)
ax.plot(J_time,pc_s_J[:,2]/(360*24*3600),color="orange",linewidth=2)

ax2=ax.twinx()
ax2.plot(J_time,LST_gb-273.15,color="red")
ax2.plot(J_time,LAI[2,:]*5,color="green")
plt.plot(J_time,C_veg_J[:,2]*100,color="brown")
plt.show()

print Cveg_ECP.shape
plt.plot(J_time,C_veg_J[:,2])
plt.plot(J_time,Cveg_ECP[2,:])
plt.show()

quit()

fig=plt.figure()

ax=fig.add_subplot(1,1,1)
ax.plot(np.sum(Soil_Resp_Fact,axis=0),ls='',marker='.')

ax2=ax.twinx()
ax2.plot(soil_resp_J,ls='',marker='.',color='green')

plt.show()

C_dpm_Q10t   = (Lit_SCpools_ECP[0,:] / Soil_Resp_Fact[0,:] ) /  \
                (3600.*24.*360.)

C_rpm_Q10t   = (Lit_SCpools_ECP[1,:] / Soil_Resp_Fact[1,:] ) /  \
                (3600.*24.*360.)
