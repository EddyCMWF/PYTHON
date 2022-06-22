#!/bin/env python3

import pdb
import numpy as np
import sys,os
import matplotlib.pyplot as plt
from imogen import data_info,ocean_co2,response,imogen_ebm

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

def dqCO2(CO2,CO2_pi,q2co2=3.74):
    dq_co2=np.log(CO2/CO2_pi) * (q2co2/np.log(2.))
    return dq_co2

GtC_to_ppm = 0.471

l_LOADPREVIOUS='-loadprevious' in sys.argv

C_emissions_file=optional_argparse('-c_emissions_file','./IMOGEN_MED_CO2_Emissions_Param.npy')

start_CO2 = float(optional_argparse('-start_CO2','284.0'))

start_year = int(optional_argparse('-start_year','1850'))
end_year = int(optional_argparse('-end_year','2100'))
now_year = int(optional_argparse('-now_year','2015'))
now_temp = float(optional_argparse('-now_temp','0.89'))
n_yr = end_year-start_year+1
yr_plot=np.arange(start_year,end_year+1)
npde=20

cmip5_runs=data_info.cmip5_runs()
n_cmip5=len(cmip5_runs)

C_emissions=np.load(C_emissions_file)
for i_gcm in range(n_cmip5):
    C_emissions[:10,i_gcm] = np.interp(np.arange(10),np.array([0,9]),C_emissions[[0,9],i_gcm])

kappa=data_info.kappa(cmip5_runs)
lambda_l=data_info.lambda_l(cmip5_runs)
lambda_o=data_info.lambda_o(cmip5_runs)
nu=data_info.nu(cmip5_runs)
f=data_info.ocean_frac(cmip5_runs)

dq_nonCO2_file=optional_argparse('-dq_nonCO2_file','dq_nonCO2.npy')
dq_nonCO2=np.load(dq_nonCO2_file)

n_o_levels = 255

#Set up arrays for output:
dq_CO2       = np.zeros([n_yr,n_cmip5])
atmos_CO2    = np.zeros([n_yr,n_cmip5])
d_atmos_CO2  = np.zeros([n_yr,n_cmip5])
d_atmos_C    = np.zeros([n_yr,n_cmip5])
ocean_CO2    = np.zeros([n_yr,n_cmip5])
d_ocean_CO2  = np.zeros([n_yr,n_cmip5])
d_ocean_C    = np.zeros([n_yr,n_cmip5])

dtemp_global = np.zeros([n_yr,n_cmip5])
#dtemp_land   = np.zeros([n_yr,n_cmip5])
dtemp_ocean  = np.zeros([n_yr,n_cmip5,n_o_levels])

# Our toy land surface 25% model:
d_land_C = C_emissions*0.25

# Set first year values:
atmos_CO2[0,:]=start_CO2

ocean_CO2[0,:]=start_CO2
dq_CO2[0,:]=dqCO2(atmos_CO2[0,:],start_CO2)

# atmoospheric CO2 changes can be calculated for the first year based on zero ocean uptake:
d_atmos_C[0,:]=C_emissions[0,:]-d_land_C[0,:]
d_atmos_CO2[0,:]=d_atmos_C[0,:]*GtC_to_ppm

dq = dq_CO2[0,:]+dq_nonCO2[0]
dtemp_ocean[0,:,:] = imogen_ebm.forward_ebm_year(dq, dtemp_ocean[0,:,:], 
                                                 kappa,f,lambda_l,lambda_o,nu)

fa_ocean = np.zeros([(n_yr+1)*npde,n_cmip5])
# Get the Green's function for use "under integral".
rs=response.response(npde,n_yr+1)

#ocean_area=3.627e14  # m2
#t_ocean_init=289.28
if l_LOADPREVIOUS:
  dtemp_global=np.load('dtemp_global.npy')
  atmos_CO2=np.load('atmos_CO2.npy')
else:
  for iyr in range(1,n_yr):
    year=start_year+iyr
    
    atmos_CO2_year = atmos_CO2[iyr-1] + (0.75*GtC_to_ppm*C_emissions[iyr,:])
    for i_gcm in range(n_cmip5):
        d_ocean_CO2[iyr,i_gcm],fa_ocean[:,i_gcm] = ocean_co2.ocean_co2_ECP(  iyr,
                                                    #atmos_CO2[iyr-1,i_gcm], ocean_CO2[iyr-1,i_gcm], 
                                                    atmos_CO2_year[i_gcm], ocean_CO2[iyr-1,i_gcm], 
                                                    start_CO2,dtemp_ocean[iyr-1,i_gcm,0],
                                                    d_atmos_CO2[iyr-1,i_gcm],
                                                    fa_ocean[:,i_gcm],rs,ncallyr=npde,
                                                    )
        #d_ocean_CO2[iyr,i_gcm] = ocean_co2.ocean_co2_ECP(  iyr,
        #                                            atmos_CO2[iyr-1,i_gcm], ocean_CO2[iyr-1,i_gcm], 
        #                                            start_CO2,dtemp_ocean[iyr-1,i_gcm,0],
        #                                            d_atmos_CO2[iyr-1,i_gcm],
        #                                            fa_ocean[:,i_gcm],rs,ncallyr=npde,
        #                                            )
    ocean_CO2[iyr,:] = ocean_CO2[iyr-1,:] + d_ocean_CO2[iyr,:]
    d_ocean_C[iyr,:] = d_ocean_CO2[iyr,:]/GtC_to_ppm
    
    d_atmos_C[iyr,:] = C_emissions[iyr,:] - d_land_C[iyr,:] + d_ocean_C[iyr,:]
    
    d_atmos_CO2[iyr,:] = d_atmos_C[iyr,:]*GtC_to_ppm

    atmos_CO2[iyr,:] = atmos_CO2[iyr-1,:] + d_atmos_CO2[iyr,:]
    
    dq_CO2[iyr,:]=dqCO2(atmos_CO2[iyr,:],start_CO2)

    dq = dq_CO2[iyr,:]+dq_nonCO2[iyr,:]

    dtemp_ocean[iyr,:,:] = imogen_ebm.forward_ebm_year(dq, dtemp_ocean[iyr-1,:,:], 
                                                            kappa,f,lambda_l,lambda_o,nu)

    dtemp_global[iyr,:] = dtemp_ocean[iyr,:,0] * (f + (1.0-f)*nu)
    # dtemp_land[iyr,:]   = 

    print(year, atmos_CO2[iyr,0], dtemp_global[iyr,0])


#pdb.set_trace()
                                        
np.save('dtemp_global.npy',dtemp_global)
np.save('atmos_CO2.npy',atmos_CO2)
print('C_emissions.shape=',C_emissions.shape)
print('atmos_CO2.shape=',atmos_CO2.shape)
print('dtemp_global.shape=',dtemp_global.shape)
print('yr_plot.shape=',yr_plot.shape)

fig,axes = plt.subplots(ncols=1,nrows=3,figsize=[15,10])

ax=axes[0]
for i_gcm in range(n_cmip5):
    ax.plot(yr_plot,C_emissions[:,i_gcm])
ax.set_ylabel('C emissions (GtC)')
ax.set_ylim([ax.get_ylim()[0],ax.get_ylim()[1]+1])
ax.grid(True)
#Now line:
ax.plot([now_year,now_year],ax.get_ylim(),c='k',lw=2)

ax=axes[1]
for i_gcm in range(n_cmip5):
    ax.plot(yr_plot[:-1],atmos_CO2[:,i_gcm])
ax.set_ylabel('Atmospheric CO2 (ppmv)')
ax.grid(True)
#Now line:
ax.plot([now_year,now_year],ax.get_ylim(),c='k',lw=2)


ax=axes[2]
for i_gcm in range(n_cmip5):
    ax.plot(yr_plot[:-1],dtemp_global[:,i_gcm])
ax.set_ylabel('Temperature Anomaly (K)')
ax.grid(True)
#Now line:
ax.plot([now_year,now_year],ax.get_ylim(),c='k',lw=2)
#Now temperature line:
ax.plot([now_year-5,now_year+5],[now_temp,now_temp],c='k',lw=2)

fig.savefig('Prescribed_Emissions.png')

plt.close()  #show()

