#!/bin/env python

import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import parabolic, profile, ocean_co2, data_info

l_outTprof=False   # Output temperature profile?
SCENARIO = '2deg'  #''1p81p5deg' #'1p81p5deg'  #'2deg'  #'1p5deg' #
Etminan=True
Baseline='SSP2-2.6_IMAGE'
version='vn1p2'   # For saving output driving files

l_Offset=False
l_GradientOffset=False
l_Gradient=False
l_DecayOffset=True

cmip5_runs = data_info.cmip5_runs()
cmip5_runs = [ cmip5_runs[8] ]
n_cmip5 = len(cmip5_runs)
print(cmip5_runs , n_cmip5)
#n_cmip5=3
#quit()

# RF equation Variables
q2co2=3.74                # Old CO2 RF parameter
# Etminan:
a1=-2.4e-7
b1=7.2e-4
c1=-2.1e-4

# Conversion factor for GtC to CO2 ppm
GtC_to_ppm=0.471

SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'

SSP_TAG = 'IMAGE-SSP2-26'
if SCENARIO=='1p5deg':
    # First set the 3 parameters in the temperature curves
    dt_limit=1.5
    mu_zero = 0.1
    mu_one = 0.000
elif SCENARIO=='1p81p5deg':
    # First set the 3 parameters in the temperature curves
    dt_limit=1.5
    mu_zero = 0.005
    mu_one = 0.00015
elif SCENARIO=='2deg':
    # First set the 3 parameters in the temperature curves
    dt_limit=2.0
    #mu_zero = 0.065
    mu_zero = 0.05
    mu_one = 0.00


if 'SSP' in Baseline.upper():
    PLOT_TAG='SSP2-2.6_Aerosol_'+SCENARIO
    CO2_FILE=SCENARIO_DIR+SSP_TAG+'_concs_co2.txt'
    AEROSOL_FILE=SCENARIO_DIR+SSP_TAG+'_qaerosol.txt'
    RAD_FORCING_FILE=SCENARIO_DIR+SSP_TAG+'_qnonco2.txt'
    CH4_N2O_FILE=SCENARIO_DIR+SSP_TAG+'_concs_ch4_n2o.txt'

    #RAD_FORCING_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2_smooth.txt'
    #RAD_FORCING_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2_1p5deg.dat'
    #CO2_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2_smooth.txt'
    #CH4_N2O_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o.txt'
    #AEROSOL_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_dq_aerosol_smooth.txt'
    #AEROSOL_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_qaerosol_1p5deg.dat'
    #AEROSOL_FILE=SCENARIO_DIR+'rcp2.6_dq_aerosol.txt'
elif 'RCP' in Baseline.upper():
    RAD_FORCING_FILE=SCENARIO_DIR+'rcp_qnonco2_vn1p0_2.6.txt'
    CO2_FILE=SCENARIO_DIR+'rcp2.6_concs_co2_vn1p0.txt'
    CH4_N2O_FILE=SCENARIO_DIR+'rcp2.6_concs_ch4_n2o_vn1p0.txt'
    AEROSOL_FILE=SCENARIO_DIR+'rcp2.6_dq_aerosol.txt'
    PLOT_TAG='RCP2.6'

#RCP_RAD_FORCING_FILE=SCENARIO_DIR+'rcp_qnonco2_vn1p0_2.6.txt'

   
# Also normalise the curves so that they end up at the current temperature and gradient estimate
beta=0.025                # K/yr
dt_now=0.89
yr_now=2015
end_year=1852
start_year=1849
n_yr=end_year-start_year+1

delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one)
delta_temp_global = np.append([0],delta_temp_global[:n_yr-1])
delta_temp_global = delta_temp_global.round(4)
#delta_temp_global[1] = 3.0206198E-02
#delta_temp_global[1] = 2.3598967E-03
#delta_temp_global[1] = 2.7100001E-02
delta_temp_global[1] = 0.0023808996093


yr_plot=np.arange(start_year,end_year+1)

#now_index=np.where(yr_plot==yr_now)[0][0]
#plt.plot(yr_plot,delta_temp_global)
#plt.grid(True)
#plt.show()
#quit()

# Output Tprofile if required
if l_outTprof:
    out_Tprof_filename=SCENARIO_DIR+SCENARIO+'_global_temp_anomaly_'+version+'.dat'
    outf=open(out_Tprof_filename,'w')
    for year,Temp in zip(yr_plot,delta_temp_global):
        line='%4i  %8.4f\n'%(year,Temp)
        outf.write(line)
    outf.close()

kappa_all=data_info.kappa(cmip5_runs)
lambda_l_all=data_info.lambda_l(cmip5_runs)
lambda_o_all=data_info.lambda_o(cmip5_runs)
nu_all=data_info.nu(cmip5_runs)
f_all=data_info.ocean_frac(cmip5_runs)

# Loop over the different GCMs emulated and calculate RF via dtemp_o
temp_gradient_out_all=np.zeros([n_yr,n_cmip5]) 
#  HadGEM-ES=17, CSIRO-Q=8, NOAA-2G=29
#if os.path.isfile(out_dir+'dq_all.npy'):
if False:
    #dq_all=np.load(out_dir+'dq_all.npy')
    print('nutin')
else:
    dq_all=np.zeros([n_yr,n_cmip5]) 
    for i_gcm in range(n_cmip5):
        print(i_gcm,':',cmip5_runs[i_gcm])
        print('f,nu,lambda_,lambda_o,kappa=',f_all[i_gcm],nu_all[i_gcm],lambda_l_all[i_gcm],lambda_o_all[i_gcm],kappa_all[i_gcm])
        delta_temp_ocean_yearly = delta_temp_global / (f_all[i_gcm] + nu_all[i_gcm] -(nu_all[i_gcm]*f_all[i_gcm]))
        temp_gradient_out = parabolic.parabolic(kappa_all[i_gcm], delta_temp_ocean_yearly)
        #print(temp_gradient_out)
        temp_gradient_out_all[:,i_gcm]=temp_gradient_out
        dq = np.zeros(n_yr)
        # Derive time-evolution of radiative forcing, Q
        factor=(((1.0-f_all[i_gcm])*lambda_l_all[i_gcm]*nu_all[i_gcm])/f_all[i_gcm])+lambda_o_all[i_gcm]
        print(temp_gradient_out)
        dq = (temp_gradient_out + delta_temp_ocean_yearly*factor)*f_all[i_gcm]
        for iyr in range(n_yr):
            print('year,dt_glob,dt_o,dq:')
            print(yr_plot[iyr], delta_temp_global[iyr], delta_temp_ocean_yearly[iyr],dq[iyr])
        dq_all[:,i_gcm]=dq
    #np.save(out_dir+'dq_all.npy',dq_all)
    #np.save(out_dir+'delta_temp_ocean.npy',delta_temp_ocean_yearly)

quit()

# Read in the non-CO2 RF from PIK site
dq_non_co2 = np.zeros(n_yr)
dq_non_co2_ssp = np.zeros(n_yr)

#ECP method from ssp2-2.6_image data
f_in = open(RAD_FORCING_FILE, 'r')
RF_lines=f_in.readlines()
f_in.close()

co2_f_in=open(CO2_FILE,'r')
CO2_lines=co2_f_in.readlines()
co2_f_in.close()

ch4_n2o_f_in=open(CH4_N2O_FILE,'r')
CH4_N2O_lines=ch4_n2o_f_in.readlines()
ch4_n2o_f_in.close()

aero_f_in=open(AEROSOL_FILE,'r')
AERO_lines=aero_f_in.readlines()
aero_f_in.close()

co2_ppm_ssp=np.zeros(n_yr)
ch4_ppb_ssp=np.zeros(n_yr)
n2o_ppb_ssp=np.zeros(n_yr)
ssp_year=np.zeros(n_yr)
dq_aerosol_ssp=np.zeros(n_yr)

for iyr in range(n_yr):
    split=RF_lines[iyr].split()
    dq_non_co2_ssp[iyr]=float(split[1])
        
    split=CO2_lines[iyr].split()
    ssp_year[iyr]=int(split[0])
    co2_ppm_ssp[iyr]=float(split[1])
    
    split=CH4_N2O_lines[iyr].split()
    ch4_ppb_ssp[iyr]=split[1]
    n2o_ppb_ssp[iyr]=split[2]

    split=AERO_lines[iyr].split()
    dq_aerosol_ssp[iyr]=split[1]

co2_ppm_pi=co2_ppm_ssp[0]
dq_other_ssp=dq_non_co2_ssp-dq_aerosol_ssp

if Etminan:
    Nbar=(n2o_ppb_ssp+n2o_ppb_ssp[0])/2.
    dq_co2_ssp=    np.log(co2_ppm_ssp/co2_ppm_pi)      \
                 * (  (a1*((co2_ppm_ssp-co2_ppm_pi)**2))      \
                    + (b1*np.abs(co2_ppm_ssp-co2_ppm_pi))     \
                    + (c1*Nbar) + 5.36 )
else:
    dq_co2_ssp= np.log(co2_ppm_ssp/co2_ppm_pi) * (q2co2/np.log(2.))

now_index=np.where(yr_plot==yr_now)[0][0]

#emiss_invert=np.zeros([n_yr,n_cmip5])     # Holds the final backed-out emissions (GtC/yr)
#  HadGEM-ES=17, CSIRO-Q=8, NOAA-2G=29
# Just Offset
if (os.path.isfile(out_dir+'co2_ppm.npy'))&(l_Offset):
    print('Loading Offset Data Radiative Forcing Data')
    dq_co2_all=np.load(out_dir+'dq_co2.npy')
    dq_aerosol_all=np.load(out_dir+'dq_aerosol.npy')
    co2_ppm_all=np.load(out_dir+'co2_ppm.npy')
elif (l_Offset):
    print('Calculating Offset Data Radiative Forcing Data')
    #array for all temp gradients and rfs
    dq_other_all=np.zeros([n_yr,n_cmip5])
    dq_co2_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_offset_all=np.zeros([n_cmip5])
    co2_ppm_all=np.zeros([n_yr,n_cmip5]) 
    for i_gcm in range(n_cmip5):
        # Now back out CO2 concentrations and non_co2 radiative forcings
        dq=np.copy(dq_all[:,i_gcm])        
        #first, historical record, non aerosol prescribed, hence dq_aerosol=dq-dq_co2_ssp-dq_other
        dq_aerosol=np.zeros(n_yr)
        dq_aerosol[:now_index+1]=dq[:now_index+1]-dq_co2_ssp[:now_index+1]-dq_other_ssp[:now_index+1]

        # store offset/scale_factor/grad_SF for yr_now
        dq_aerosol_offset           = dq_aerosol[now_index]-dq_aerosol_ssp[now_index]
        dq_aerosol_offset_all[i_gcm]= np.copy(dq_aerosol_offset)
        # for remaining period dq_non_co2 is the ssp value - offset
        dq_aerosol[now_index:]=dq_aerosol_ssp[now_index:] + dq_aerosol_offset
        dq_aerosol_all[:,i_gcm]=np.copy(dq_aerosol)
        
        # Can now Calculate CO2 from dq and dq_non_co2, 
        # we do whole time period to check conservation
        dq_co2 = dq - dq_aerosol - dq_other_ssp
        dq_co2_all[:,i_gcm]=np.copy(dq_co2)

        if Etminan:
            co2_ppm_pi=co2_ppm_ssp[0]
            co2_ppm=np.copy(co2_ppm_ssp)
            for iyr in range(n_yr):
                Nbar=(n2o_ppb_ssp[iyr]+n2o_ppb_ssp[0])*0.5
                # iterate to find CO2 concentration from radiative forcing.
                co2_OF = np.copy(co2_ppm[iyr-1])
                itercnt=0
                co2_iter=999.
                while ( (np.abs(co2_iter-co2_OF)>0.001) & (itercnt<=1000) ):
                    co2_iter = co2_OF.copy()
                    denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     \
                            + (b1*(co2_iter-co2_ppm_pi))           \
                            + (c1*Nbar)+5.36)
                    
                    co2_OF = co2_ppm_pi * np.exp( dq_co2[iyr]/denom )
                    itercnt+=1
                co2_ppm[iyr]=np.copy(co2_OF)
                    
        else:
            co2_ppm = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2)/q2co2)

        co2_ppm_all[:,i_gcm]=np.copy(co2_ppm)
    
    np.save(out_dir+'dq_co2.npy',dq_co2_all)
    np.save(out_dir+'dq_aerosol.npy',dq_aerosol_all)
    np.save(out_dir+'co2_ppm.npy',co2_ppm_all)

# Gradient Offset
if (os.path.isfile(out_dir+'co2_ppm_GOS.npy'))&(l_GradientOffset):
    print('Loading Gradient-Offset Data Radiative Forcing Data')
    dq_co2_GOS_all=np.load(out_dir+'dq_co2_GOS.npy')
    dq_aerosol_GOS_all=np.load(out_dir+'dq_aerosol_GOS.npy')
    co2_ppm_GOS_all=np.load(out_dir+'co2_ppm_GOS.npy')
elif (l_GradientOffset): 
    print('Calculating Gradient-Offset Data Radiative Forcing Data')
    # For Gradient Offset
    dq_other_GOS_all=np.zeros([n_yr,n_cmip5])
    dq_co2_GOS_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_GOS_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_offset_GOS_all=np.zeros([n_cmip5])
    co2_ppm_GOS_all=np.zeros([n_yr,n_cmip5])
    for i_gcm in range(n_cmip5):
        # Now back out CO2 concentrations and non_co2 radiative forcings
        dq=np.copy(dq_all[:,i_gcm])        
        #first, historical record, non aerosol prescribed, hence dq_aerosol=dq-dq_co2_ssp-dq_other
        dq_aerosol_GOS=np.zeros(n_yr)
        dq_aerosol_GOS[:now_index+1]=dq[:now_index+1]-dq_co2_ssp[:now_index+1]-dq_other_ssp[:now_index+1]

        # Or Gradient Offset
        dqim_dt = dq_aerosol_GOS[now_index]-dq_aerosol_GOS[now_index-1]
        dqssp_dt = dq_aerosol_ssp[now_index]-dq_aerosol_ssp[now_index-1]
        dq_aerosol_GOS_c = (dqim_dt/dqssp_dt)
        dq_aerosol_GOS_d = dq_aerosol_GOS[now_index] - dq_aerosol_ssp[now_index]*dq_aerosol_GOS_c
        dq_aerosol_GOS[now_index:] = (dq_aerosol_ssp[now_index:]*dq_aerosol_GOS_c) + dq_aerosol_GOS_d 
        dq_aerosol_GOS_all[:,i_gcm]=np.copy(dq_aerosol_GOS)
       
        # Can now Calculate CO2 from dq and dq_non_co2, 
        # we do whole time period to check conservation
        dq_co2_GOS = dq - dq_aerosol_GOS - dq_other_ssp
        dq_co2_GOS_all[:,i_gcm]=np.copy(dq_co2_GOS)

        if Etminan:
            co2_ppm_pi=co2_ppm_ssp[0]
            co2_ppm_GOS=np.copy(co2_ppm_ssp)
            for iyr in range(n_yr):
                Nbar=(n2o_ppb_ssp[iyr]+n2o_ppb_ssp[0])*0.5
                # iterate to find CO2 concentration from radiative forcing.
                co2_GOS = np.copy(co2_ppm_GOS[iyr-1])
                itercnt=0
                co2_iter=999.
                while ( (np.abs(co2_iter-co2_GOS)>0.001) & (itercnt<=1000) ):
                    co2_iter = co2_GOS.copy()
                    denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     \
                            + (b1*(co2_iter-co2_ppm_pi))           \
                            + (c1*Nbar)+5.36)
                    co2_GOS = co2_ppm_pi * np.exp( dq_co2_GOS[iyr]/denom )
                    itercnt+=1
                co2_ppm_GOS[iyr]=np.copy(co2_GOS)
        else:
            co2_ppm_GOS = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2_GOS)/q2co2)

        co2_ppm_GOS_all[:,i_gcm]=np.copy(co2_ppm_GOS)
    
    np.save(out_dir+'dq_co2_GOS.npy',dq_co2_GOS_all)
    np.save(out_dir+'dq_aerosol_GOS.npy',dq_aerosol_GOS_all)
    np.save(out_dir+'co2_ppm_GOS.npy',co2_ppm_GOS_all)

if (os.path.isfile(out_dir+'co2_ppm_G.npy'))&(l_Gradient):
#if False:
    print('Loading Gradient Data Radiative Forcing Data')
    dq_co2_G_all=np.load(out_dir+'dq_co2_G.npy')
    dq_aerosol_G_all=np.load(out_dir+'dq_aerosol_G.npy')
    co2_ppm_G_all=np.load(out_dir+'co2_ppm_G.npy')
elif (l_Gradient):
    print('Calculating Gradient Data Radiative Forcing Data')
    # For Gradient
    dq_other_G_all=np.zeros([n_yr,n_cmip5])
    dq_co2_G_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_G_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_offset_G_all=np.zeros([n_cmip5])
    co2_ppm_G_all=np.zeros([n_yr,n_cmip5])

    for i_gcm in range(n_cmip5):
        # Now back out CO2 concentrations and non_co2 radiative forcings
        dq=np.copy(dq_all[:,i_gcm])        
        #first, historical record, non aerosol prescribed, hence dq_aerosol=dq-dq_co2_ssp-dq_other
        dq_aerosol_G=np.zeros(n_yr)
        dq_aerosol_G[:now_index+1]=dq[:now_index+1]-dq_co2_ssp[:now_index+1]-dq_other_ssp[:now_index+1]

        # Gradient and timed Offset (CHG method)
        dqim_dt = dq_aerosol_G[now_index]-dq_aerosol_G[now_index-1]
        dqssp_dt = dq_aerosol_ssp[now_index]-dq_aerosol_ssp[now_index-1]
        
        #dq_aerosol_G_m = dq_aerosol_G[now_index]/dq_aerosol_ssp[now_index]
        #dq_aerosol_G_c = dqim_dt - (dqssp_dt*dq_aerosol_G_m)
        #dq_aerosol_G[now_index:] = (dq_aerosol_ssp[now_index:]*dq_aerosol_G_m) +  \
        #                            (dq_aerosol_G_c*(yr_plot[now_index:]-yr_now))
        
        dq_aerosol_G_m = dq_aerosol_G[now_index]/dq_aerosol_ssp[now_index]
        dq_aerosol_G_c = (dqim_dt - (dqssp_dt*dq_aerosol_G_m))*1/dq_aerosol_ssp[now_index]
        dq_aerosol_G[now_index:] = (dq_aerosol_G_m+(dq_aerosol_G_c*(yr_plot[now_index:]-yr_now))) \
                                    * dq_aerosol_ssp[now_index:]

        dq_aerosol_G_all[:,i_gcm]=np.copy(dq_aerosol_G)

        dq_co2_G = dq - dq_aerosol_G - dq_other_ssp
        dq_co2_G_all[:,i_gcm]=np.copy(dq_co2_G)
        #pdb.set_trace()
        if Etminan:
            co2_ppm_pi=co2_ppm_ssp[0]
            co2_ppm_G=np.copy(co2_ppm_ssp)
            for iyr in range(n_yr):
                Nbar=(n2o_ppb_ssp[iyr]+n2o_ppb_ssp[0])*0.5
                # iterate to find CO2 concentration from radiative forcing.
                co2_G = np.copy(co2_ppm_G[iyr-1])
                itercnt=0
                co2_iter=999.
                while ( (np.abs(co2_iter-co2_G)>0.001) & (itercnt<=1000) ):
                    co2_iter = co2_G.copy()
                    denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     \
                            + (b1*(co2_iter-co2_ppm_pi))           \
                            + (c1*Nbar)+5.36)
                    co2_G = co2_ppm_pi * np.exp( dq_co2_G[iyr]/denom )
                    itercnt+=1
                co2_ppm_G[iyr]=np.copy(co2_G)
        else:
            co2_ppm_G = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2_G)/q2co2)
        
        co2_ppm_G_all[:,i_gcm]=np.copy(co2_ppm_G)
    
    np.save(out_dir+'dq_co2_G.npy',dq_co2_G_all)
    np.save(out_dir+'dq_aerosol_G.npy',dq_aerosol_G_all)
    np.save(out_dir+'co2_ppm_G.npy',co2_ppm_G_all)
        
# Decay Offset
# fit to dq_aerosol=dq_aerosol_SSP + Ae^-Bt
if (os.path.isfile(out_dir+'co2_ppm_DOS.npy'))&(l_DecayOffset):
    print('Loading Decaying Offset Data Radiative Forcing Data')
    dq_co2_DOS_all=np.load(out_dir+'dq_co2_DOS.npy')
    dq_aerosol_DOS_all=np.load(out_dir+'dq_aerosol_DOS.npy')
    co2_ppm_DOS_all=np.load(out_dir+'co2_ppm_DOS.npy')
elif (l_DecayOffset):
    print('Calculating Decaying Offset Data Radiative Forcing Data')
    # For Decaying Offset
    dq_other_DOS_all=np.zeros([n_yr,n_cmip5])
    dq_co2_DOS_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_DOS_all=np.zeros([n_yr,n_cmip5]) 
    dq_aerosol_offsetA_DOS_all=np.zeros([n_cmip5])
    dq_aerosol_offsetB_DOS_all=np.zeros([n_cmip5])
    co2_ppm_DOS_all=np.zeros([n_yr,n_cmip5])
    for i_gcm in range(n_cmip5):
        #print(i_gcm)
        # Now back out CO2 concentrations and non_co2 radiative forcings
        dq=np.copy(dq_all[:,i_gcm])        
        #first, historical record, non aerosol prescribed, hence dq_aerosol=dq-dq_co2_ssp-dq_other
        dq_aerosol_DOS=np.zeros(n_yr)
        dq_aerosol_DOS[:now_index+1]=dq[:now_index+1]-dq_co2_ssp[:now_index+1]-dq_other_ssp[:now_index+1]

        # Or Gradient Offset
        dqim_dt = dq_aerosol_DOS[now_index]-dq_aerosol_DOS[now_index-1]
        dqssp_dt = dq_aerosol_ssp[now_index]-dq_aerosol_ssp[now_index-1]
        t = yr_plot[now_index:]-yr_now

        #dq_aerosol_DOS_A = dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index]
        #dq_aerosol_DOS_B=np.abs((dqssp_dt-dqim_dt) / dq_aerosol_DOS_A)
        #dq_aerosol_DOS_B = (dqssp_dt-dqim_dt) / dq_aerosol_DOS_A
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:]   \
        #         + ( dq_aerosol_DOS_A*np.exp(-dq_aerosol_DOS_B*(t)) )
        
        #dq_aerosol_DOS_A=(dq_aerosol_DOS[now_index]/dq_aerosol_ssp[now_index])-1.
        #dq_aerosol_DOS_B=( (1./(dq_aerosol_DOS_A*dq_aerosol_ssp[now_index])) \
        #      *  ( ((1+dq_aerosol_DOS_A)*dqssp_dt) - dqim_dt ) ) \
        #      + (dqssp_dt/dq_aerosol_ssp[now_index])
        #m_t = 1 + ( dq_aerosol_DOS_A*np.exp(-dq_aerosol_DOS_B*(t)) )
        #dq_aerosol_DOS[now_index:] =dq_aerosol_ssp[now_index:]+ m_t *  dq_aerosol_ssp[now_index:] 

        #beta1=0.045
        #beta1= np.abs((dqim_dt/dqssp_dt)-1)
        #dq_aerosol_DOS_B =  (dqim_dt + (beta1*dq_aerosol_DOS[now_index])) \
        #                        / (dqssp_dt + (beta1*dq_aerosol_ssp[now_index]))
        #dq_aerosol_DOS_A = ( dq_aerosol_DOS[now_index] \
        #                          - (dq_aerosol_DOS_B*dq_aerosol_ssp[now_index]) ) \
        #                         * np.exp(1.)
        #c_t = dq_aerosol_DOS_A * np.exp( -beta1*(t-1.)
        #dq_aerosol_DOS[now_index:] = (dq_aerosol_DOS_B*dq_aerosol_ssp[now_index:]) + c_t
        
        ##  c_t = A + B*t*e^(-beta*t)
        #beta1=0.05
        #dq_aerosol_DOS_A = dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index]
        #dq_aerosol_DOS_B = dqim_dt - dqssp_dt
        #c_t = dq_aerosol_DOS_A+(dq_aerosol_DOS_B*t*np.exp(-beta1*t))
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + c_t
        
        #  c_t = A + B*e^(-beta*t) 
        #beta1=0.1
        #dq_aerosol_DOS_B = (-1./beta1)*(dqim_dt - dqssp_dt)
        #dq_aerosol_DOS_A = dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index]-dq_aerosol_DOS_B
        #c_t = dq_aerosol_DOS_A+(dq_aerosol_DOS_B*np.exp(-beta1*t))
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + c_t
       
        # c_t = A(*t+1)*e^(-B*t)
        #dq_aerosol_DOS_A = dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index]
        #dq_aerosol_DOS_B = 1. - ( (dqim_dt - dqssp_dt)/dq_aerosol_DOS_A )
        #c_t = dq_aerosol_DOS_A * ( t+1. ) * np.exp(-dq_aerosol_DOS_B*t)
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + c_t
        
        ## c_t = A/(t+1)+B*e^(-beta*t)
        #beta1=0.001
        #dq_aerosol_DOS_B = (1./(1.-beta1))  *          \
        #        ( dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index]+dqim_dt - dqssp_dt )
        #dq_aerosol_DOS_A = dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index] -dq_aerosol_DOS_B 
        #c_t = (dq_aerosol_DOS_A /( t+1. )) +(dq_aerosol_DOS_B*np.exp(-beta1*t))
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + c_t
       
        # Q = Qssp + c_t*Qssp
        # c_t= A + B*e^(-beta*t)
        #beta1=0.001
        #dq_aerosol_DOS_B = (1./(beta1*dq_aerosol_ssp[now_index]))  *          \
        #        ( (dqssp_dt*dq_aerosol_DOS[now_index]/dq_aerosol_ssp[now_index]) - dqim_dt  )
        #dq_aerosol_DOS_A = (dq_aerosol_DOS[now_index]/dq_aerosol_ssp[now_index]) - 1 - dq_aerosol_DOS_B 
        #c_t = ( dq_aerosol_DOS_A + ( dq_aerosol_DOS_B*np.exp(-beta1*t)) )
        #dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + (c_t*dq_aerosol_ssp[now_index:] )

        # Q = Qssp + c_t*Qssp
        # c_t= A + B*t*e^(-beta*t)
        beta1=0.01
        dq_aerosol_DOS_A = (dq_aerosol_DOS[now_index]-dq_aerosol_ssp[now_index])/dq_aerosol_ssp[now_index]
        dq_aerosol_DOS_B = (1./dq_aerosol_ssp[now_index])  *          \
                ( dqim_dt-dqssp_dt -(dq_aerosol_DOS_A*dqssp_dt)  )
        c_t = ( dq_aerosol_DOS_A + ( dq_aerosol_DOS_B*t*np.exp(-beta1*t)) )
        dq_aerosol_DOS[now_index:] = dq_aerosol_ssp[now_index:] + (c_t*dq_aerosol_ssp[now_index:] )
        

        dq_aerosol_DOS_all[:,i_gcm]=np.copy(dq_aerosol_DOS)
        
        dq_aerosol_offsetA_DOS_all[i_gcm]=dq_aerosol_DOS_A
        dq_aerosol_offsetB_DOS_all[i_gcm]=dq_aerosol_DOS_B

        # Can now Calculate CO2 from dq and dq_non_co2, 
        # we do whole time period to check conservation
        dq_co2_DOS = dq - dq_aerosol_DOS - dq_other_ssp
        dq_co2_DOS_all[:,i_gcm]=np.copy(dq_co2_DOS)

        if Etminan:
            co2_ppm_pi=co2_ppm_ssp[0]
            co2_ppm_DOS=np.copy(co2_ppm_ssp)
            for iyr in range(n_yr):
                Nbar=(n2o_ppb_ssp[iyr]+n2o_ppb_ssp[0])*0.5
                # iterate to find CO2 concentration from radiative forcing.
                co2_DOS = np.copy(co2_ppm_DOS[iyr-1])
                itercnt=0
                co2_iter=999.
                while ( (np.abs(co2_iter-co2_DOS)>0.001) & (itercnt<=1000) ):
                    co2_iter = co2_DOS.copy()
                    denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     \
                            + (b1*(co2_iter-co2_ppm_pi))           \
                            + (c1*Nbar)+5.36)
                    co2_DOS = co2_ppm_pi * np.exp( dq_co2_DOS[iyr]/denom )
                    itercnt+=1
                co2_ppm_DOS[iyr]=np.copy(co2_DOS)
        else:
            co2_ppm_DOS = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2_DOS)/q2co2)

        co2_ppm_DOS_all[:,i_gcm]=np.copy(co2_ppm_DOS)
    
    np.save(out_dir+'dq_co2_DOS.npy',dq_co2_DOS_all)
    np.save(out_dir+'dq_aerosol_DOS.npy',dq_aerosol_DOS_all)
    np.save(out_dir+'co2_ppm_DOS.npy',co2_ppm_DOS_all)

        
if (os.path.isfile(out_dir+'d_ocean_atmos.npy'))&(l_Offset):
    print('Loading Emissions - Offset')
    d_ocean_atmos_all=np.load(out_dir+'d_ocean_atmos.npy')
    co2_emissions_all=np.load(out_dir+'co2_emissions.npy')
elif (l_Offset):
    d_ocean_atmos_all=np.zeros([n_yr,n_cmip5]) 
    co2_emissions_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions - Offset')
    delta_temp_ocean_yearly_all = np.array([delta_temp_global / \
                                   (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)
    d_ocean_atmos_all     = ocean_co2.ocean_co2_parallel(co2_ppm_all,
                                                         delta_temp_ocean_yearly_all)
    d_ocean_atmos_GtC_all = d_ocean_atmos_all/GtC_to_ppm
    delta_co2_ppm_all     = np.zeros_like(co2_ppm_all)
    delta_co2_ppm_all[1:] = co2_ppm_all[1:,:]-co2_ppm_all[:-1,:]
    delta_co2_GtC_all     = delta_co2_ppm_all/GtC_to_ppm
    co2_emissions_all     = (delta_co2_GtC_all-d_ocean_atmos_GtC_all)/0.75
    np.save(out_dir+'d_ocean_atmos.npy',d_ocean_atmos_all)
    np.save(out_dir+'co2_emissions.npy',co2_emissions_all)

if (os.path.isfile(out_dir+'d_ocean_atmos_GOS.npy'))&(l_GradientOffset):
    print('Loading Emissions - GradientOffset')
    d_ocean_atmos_GOS_all=np.load(out_dir+'d_ocean_atmos_GOS.npy')
    co2_emissions_GOS_all=np.load(out_dir+'co2_emissions_GOS.npy')
elif (l_GradientOffset):
    d_ocean_atmos_GOS_all=np.zeros([n_yr,n_cmip5]) 
    co2_emissions_GOS_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions - GradientOffset')
    delta_temp_ocean_yearly_all = np.array([delta_temp_global / \
                                   (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)
    d_ocean_atmos_GOS_all     = ocean_co2.ocean_co2_parallel(co2_ppm_GOS_all,
                                                         delta_temp_ocean_yearly_all)
    d_ocean_atmos_GOS_GtC_all = d_ocean_atmos_GOS_all/GtC_to_ppm
    delta_co2_GOS_ppm_all     = np.zeros_like(co2_ppm_GOS_all)
    delta_co2_GOS_ppm_all[1:] = co2_ppm_GOS_all[1:,:]-co2_ppm_GOS_all[:-1,:]
    delta_co2_GOS_GtC_all     = delta_co2_GOS_ppm_all/GtC_to_ppm
    co2_emissions_GOS_all     = (delta_co2_GOS_GtC_all-d_ocean_atmos_GOS_GtC_all)/0.75
    np.save(out_dir+'d_ocean_atmos_GOS.npy',d_ocean_atmos_GOS_all)
    np.save(out_dir+'co2_emissions_GOS.npy',co2_emissions_GOS_all)

if (os.path.isfile(out_dir+'d_ocean_atmos_G.npy'))&(l_Gradient):
    print('Loading Emissions - Gradient')
    d_ocean_atmos_G_all=np.load(out_dir+'d_ocean_atmos_G.npy')
    co2_emissions_G_all=np.load(out_dir+'co2_emissions_G.npy')
elif (l_Gradient):
    d_ocean_atmos_G_all=np.zeros([n_yr,n_cmip5]) 
    co2_emissions_G_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions - Gradient')
    delta_temp_ocean_yearly_all = np.array([delta_temp_global / \
                                   (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)
    d_ocean_atmos_G_all     = ocean_co2.ocean_co2_parallel(co2_ppm_G_all,
                                                             delta_temp_ocean_yearly_all)
    d_ocean_atmos_G_GtC_all = d_ocean_atmos_G_all/GtC_to_ppm
    delta_co2_G_ppm_all     = np.zeros_like(co2_ppm_G_all)
    delta_co2_G_ppm_all[1:] = co2_ppm_G_all[1:]-co2_ppm_G_all[:-1]
    delta_co2_G_GtC_all     = delta_co2_G_ppm_all/GtC_to_ppm
    co2_emissions_G_all     = (delta_co2_G_GtC_all-d_ocean_atmos_G_GtC_all)/0.75
    np.save(out_dir+'d_ocean_atmos_G.npy',d_ocean_atmos_G_all)
    np.save(out_dir+'co2_emissions_G.npy',co2_emissions_G_all)

if (os.path.isfile(out_dir+'d_ocean_atmos_DOS.npy'))&(l_DecayOffset):
    print('Loading Emissions - DecayOffset')
    d_ocean_atmos_DOS_all=np.load(out_dir+'d_ocean_atmos_DOS.npy')
    co2_emissions_DOS_all=np.load(out_dir+'co2_emissions_DOS.npy')
elif (l_DecayOffset):
    d_ocean_atmos_DOS_all=np.zeros([n_yr,n_cmip5]) 
    co2_emissions_DOS_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions - DecayingOffset')
    delta_temp_ocean_yearly_all = np.array([delta_temp_global / \
                                   (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)
    d_ocean_atmos_DOS_all     = ocean_co2.ocean_co2_parallel(co2_ppm_DOS_all,
                                                         delta_temp_ocean_yearly_all)
    d_ocean_atmos_DOS_GtC_all = d_ocean_atmos_DOS_all/GtC_to_ppm
    delta_co2_DOS_ppm_all     = np.zeros_like(co2_ppm_DOS_all)
    delta_co2_DOS_ppm_all[1:] = co2_ppm_DOS_all[1:,:]-co2_ppm_DOS_all[:-1,:]
    delta_co2_DOS_GtC_all     = delta_co2_DOS_ppm_all/GtC_to_ppm
    co2_emissions_DOS_all     = (delta_co2_DOS_GtC_all-d_ocean_atmos_DOS_GtC_all)/0.75
    np.save(out_dir+'d_ocean_atmos_DOS.npy',d_ocean_atmos_DOS_all)
    np.save(out_dir+'co2_emissions_DOS.npy',co2_emissions_DOS_all)


METHOD_DICT={}
if l_Offset:
    METHOD_DICT['Offset']={'dq_co2':dq_co2_all,        
                           'dq_aerosol':dq_aerosol_all,
                           'dq_nonco2':dq_aerosol_all+dq_non_co2_ssp-dq_aerosol_ssp,
                           'co2_ppm':co2_ppm_all,      
                           'co2_emissions':co2_emissions_all,      
                           'colour':'b'}

if l_GradientOffset:
    METHOD_DICT['ScaleFactor_Offset']={'dq_co2':dq_co2_GOS_all,        
                                       'dq_aerosol':dq_aerosol_GOS_all,
                                       'dq_nonco2':dq_aerosol_GOS_all+dq_non_co2_ssp-dq_aerosol_ssp,
                                       'co2_ppm':co2_ppm_GOS_all,      
                                       'co2_emissions':co2_emissions_GOS_all,      
                                       'colour':'green'}

if l_Gradient:
    METHOD_DICT['ScaleFactor']={'dq_co2':dq_co2_G_all,        
                                'dq_aerosol':dq_aerosol_G_all,
                                'dq_nonco2':dq_aerosol_G_all+dq_non_co2_ssp-dq_aerosol_ssp,
                                'co2_ppm':co2_ppm_G_all,      
                                'co2_emissions':co2_emissions_G_all,      
                                'colour':'darkorange'}

if l_DecayOffset:
    dq_nonco2_DOS_all = np.array( [ dq_non_co2_ssp-dq_aerosol_ssp+dq_aerosol_DOS_all[:,i_gcm] 
                                    for i_gcm in range(n_cmip5) ] ).transpose(1,0)
    METHOD_DICT['Decaying_Offset']={'dq_co2':dq_co2_DOS_all,
                                    'dq_aerosol':dq_aerosol_DOS_all,
                                    'dq_nonco2':dq_nonco2_DOS_all,
                                    'co2_ppm':co2_ppm_DOS_all,
                                    'co2_emissions':co2_emissions_DOS_all,
                                    'colour':'orchid'}

plot_index=range(n_cmip5)
#plot_index=np.where(np.max(dq_all,axis=0)>2.5)[0]

FONTSIZE=35

for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,28],ncols=1,nrows=5)
    colour=METHOD_DICT[method]['colour']
    # Plot Temperature Profile 
    ax=AXES[0]
    ax.plot(ssp_year,delta_temp_global,c='g',lw=3)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Global Temperature Anomaly',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    # Plot total RF in top plot
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(ssp_year,dq_all[:,i_gcm],c=colour)
    ax.plot(ssp_year,dq_non_co2_ssp+dq_co2_ssp,c='r',lw=3)
    ax.plot(ssp_year,delta_temp_global,c='k',lw=3,ls=':')
    ax.set_title('Total Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[2]
    for i_gcm in plot_index:   #range(n_cmip5):
        #ax.plot(yr_plot,METHOD_DICT[method]['dq_nonco2'][:,i_gcm],c=colour,ls=':',lw=0.8)
        ax.plot(yr_plot,METHOD_DICT[method]['dq_aerosol'][:,i_gcm],c=colour)
    ax.plot(yr_plot,dq_aerosol_ssp, c='r',lw=3)
    ax.plot(yr_plot,dq_other_ssp, c='r',lw=3,ls=':')
    ax.plot(yr_plot,dq_non_co2_ssp, c='r',lw=3,ls='--')
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Aerosol Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fourth row
    ax=AXES[3]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_ppm'][:,i_gcm],c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Atmospheric CO2',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Concentration (ppmv)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fifth row
    ax=AXES[4]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_emissions'][:,i_gcm],c=colour)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Anthropogenic CO2 Emissions',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Emission (GtC yr$^{-1}$)',fontsize=FONTSIZE/2.)
    
    FIG.savefig(out_dir+method+'_RFbreakdown.png',bbox_inches='tight')
    FIG.savefig(out_dir+method+'_RFbreakdown.eps',bbox_inches='tight')
    plt.close()

print(np.max(dq_all[now_index,:]),np.min(dq_all[now_index,:]))
print(np.max(dq_all[now_index,:])-np.min(dq_all[now_index,:]))
FONTSIZE=30.
for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,20],ncols=1,nrows=4)
    colour=METHOD_DICT[method]['colour']
    # Plot total RF in top plot
    ax=AXES[0]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(ssp_year,dq_all[:,i_gcm],c=colour)
    ax.plot(ssp_year,dq_non_co2_ssp+dq_co2_ssp,c='r',lw=3)
    #ax.plot(ssp_year,delta_temp_global,c='k',lw=3,ls=':')
    ax.set_ylabel('Total Radiative Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.grid(True)
    ax.text(yr_plot[2],((ax.get_ylim()[1]-ax.get_ylim()[0])*0.95)+ax.get_ylim()[0],
            '(a)', fontsize=FONTSIZE,verticalalignment='top')

    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_aerosol'][:,i_gcm],c=colour)
    ax.plot(yr_plot,dq_aerosol_ssp, c='r',lw=3)
    ax.plot(yr_plot,dq_other_ssp, c='r',lw=3,ls=':')
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_ylabel('Non-CO$_2$ Radiative Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.text(yr_plot[2],((ax.get_ylim()[1]-ax.get_ylim()[0])*0.95)+ax.get_ylim()[0],
            '(b)', fontsize=FONTSIZE,verticalalignment='top')
    
    # Plot CO2 radiative forcing on third row
    ax=AXES[2]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_co2'][:,i_gcm],c=colour)
    ax.plot(yr_plot,dq_co2_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_ylabel('CO$_2$ Radiative Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.text(yr_plot[2],((ax.get_ylim()[1]-ax.get_ylim()[0])*0.95)+ax.get_ylim()[0],
            '(c)', fontsize=FONTSIZE,verticalalignment='top')
    
    # Plot CO2 concentration on fourth row
    ax=AXES[3]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_ppm'][:,i_gcm],c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_ylabel('Atmospheric CO$_2$ (ppmv)',fontsize=FONTSIZE/2.)
    ax.text(yr_plot[2],((ax.get_ylim()[1]-ax.get_ylim()[0])*0.95)+ax.get_ylim()[0],
            '(d)',fontsize=FONTSIZE,verticalalignment='top')

    #FIG.suptitle(PLOT_TAG+', Imogen vs '+PLOT_TAG+' '+method,fontsize=FONTSIZE)
    FIG.savefig(out_dir+method+'_RFbreakdown_forPaper.png',bbox_inches='tight')
    FIG.savefig(out_dir+method+'_RFbreakdown_forPaper.eps',bbox_inches='tight')
    plt.close()


