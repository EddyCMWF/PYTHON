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

CO2_file=optional_argparse('-co2_file','/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/SSP2-2.6_IMAGE_concs_co2_vn2p0.txt')
Q_NONCO2_file=optional_argparse('-qnonco2_file','/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/SSP2-2.6_IMAGE_qnonco2_vn2p0.txt')

start_CO2 = float(optional_argparse('-start_CO2','284.7250'))

start_year = int(optional_argparse('-start_year','1850'))
end_year = int(optional_argparse('-end_year','2100'))
now_year = int(optional_argparse('-now_year','2015'))
now_temp = float(optional_argparse('-now_temp','0.89'))
n_yr = end_year-start_year+1
yr_plot=np.arange(start_year,end_year+1)
npde=20

cmip5_runs = data_info.cmip5_runs()
cmip5_runs = [ cmip5_runs[8] ]
n_cmip5 = len(cmip5_runs)
print(cmip5_runs , n_cmip5)
#n_cmip5=3

CO2_lines=open(CO2_file).readlines()
CO2_year,CO2=[],[]
for line in CO2_lines:     
    split=line.split()
    CO2_year.append(int(split[0]))
    CO2.append(float(split[1]))
CO2_year,CO2 = np.array(CO2_year),np.array(CO2)

Q_NONCO2_lines=open(Q_NONCO2_file).readlines()
Q_NONCO2_year,Q_NONCO2=[],[]
for line in Q_NONCO2_lines:     
    split=line.split()
    Q_NONCO2_year.append(int(split[0]))
    Q_NONCO2.append(float(split[1]))
Q_NONCO2_year,Q_NONCO2 = np.array(Q_NONCO2_year),np.array(Q_NONCO2)
#Q_NONCO2[0]=1.0254369E-02

kappa=data_info.kappa(cmip5_runs)
lambda_l=data_info.lambda_l(cmip5_runs)
lambda_o=data_info.lambda_o(cmip5_runs)
nu=data_info.nu(cmip5_runs)
f=data_info.ocean_frac(cmip5_runs)

n_o_levels = 255

#Set up arrays for output:
dq_CO2=dqCO2(CO2,start_CO2)

dtemp_global = np.zeros([n_yr,n_cmip5])
#dtemp_land   = np.zeros([n_yr,n_cmip5])
dtemp_ocean  = np.zeros([n_yr,n_cmip5,n_o_levels])

# Set first year values:


# atmoospheric CO2 changes can be calculated for the first year based on zero ocean uptake:

dq = dq_CO2[0]+Q_NONCO2[0]
if n_cmip5==1:
    dq=np.array([dq])
print(dq)
#print(dtemp_ocean.shape)
iyr=0
print('dq,dq_CO2,dq_nonCO2',dq,dq_CO2[iyr],Q_NONCO2[iyr])
print('kappa,f,lambda_l,lambda_o,nu',kappa,f,lambda_l,lambda_o,nu)
dtemp_ocean[0,:,:] = imogen_ebm.forward_ebm_year(dq, dtemp_ocean[0,:,:], 
                                                 kappa,f,lambda_l,lambda_o,nu)
dtemp_global[0,:] = dtemp_ocean[0,:,0] * (f + (1.0-f)*nu)
print('dtemp_ocean',dtemp_ocean[0,0,:10])
print('dtemp_global',dtemp_global[iyr,0])

fa_ocean = np.zeros([(n_yr+1)*npde,n_cmip5])
# Get the Green's function for use "under integral".
rs=response.response(npde,n_yr+1)

#ocean_area=3.627e14  # m2
#t_ocean_init=289.28
if l_LOADPREVIOUS:
    print('Nothing')
else:
  for iyr in range(1,n_yr):
    year=start_year+iyr
    print(year)
    
    dq = dq_CO2[iyr]+Q_NONCO2[iyr]
    if n_cmip5==1:
        dq=np.array([dq])
     
    print('dq,dq_CO2,dq_nonCO2',dq,dq_CO2[iyr],Q_NONCO2[iyr])

    dtemp_ocean[iyr,:,:] = imogen_ebm.forward_ebm_year(dq, dtemp_ocean[iyr-1,:,:], 
                                                            kappa,f,lambda_l,lambda_o,nu)
    print('dtemp_ocean',dtemp_ocean[iyr,0,:10])
    dtemp_global[iyr,:] = dtemp_ocean[iyr,:,0] * (f + (1.0-f)*nu)
    print('dtemp_global',dtemp_global[iyr,0])
    quit()


#pdb.set_trace()
quit()
np.save('dtemp_global_test.npy',dtemp_global)
np.save('atmos_CO2_test.npy',atmos_CO2)
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

