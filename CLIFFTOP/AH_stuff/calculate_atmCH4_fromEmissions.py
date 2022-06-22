#!/bin/env python2.7

import pdb
import numpy as np
import matplotlib.pyplot as plt
from imogen import data_info,delQ



data_dir='./' 
runid='BL'
scenario='1p5deg'
iscenario=2

data_start_year=1850
data_end_year=2100
nYEARs=data_end_year-data_start_year+1


atmCO2_in=np.load(data_dir+runid+'atmCO2.npy') # ppbv
atmCO2=atmCO2_in[:,iscenario,:nYEARs]

fCH4_startyear=2000

fCH4_file_1 = data_dir+'methane_emissions_SSP2_RCP19_1p5degTgC.dat'
fCH4_in_1=np.loadtxt(fCH4_file_1)
fCH4_sf_1 = 180./fCH4_in_1[0]
fCH4_file_2 = data_dir+'methane_emissions_SSP2_SPA2_RCP26_1p5degTgC.dat'
fCH4_in_2=np.loadtxt(fCH4_file_2)
fCH4_sf_2 = 180./fCH4_in_2[0]

fCH4_TgCH4_1 = np.zeros_like(atmCO2)+180.
fCH4_TgCH4_1[:,fCH4_startyear-data_start_year:] = fCH4_in_1*fCH4_sf_1
fCH4_TgCH4_1[:,-1]=fCH4_TgCH4_1[:,-2]
fCH4_TgCH4_2 = np.zeros_like(atmCO2)+180.
fCH4_TgCH4_2[:,fCH4_startyear-data_start_year:] = fCH4_in_2*fCH4_sf_2
fCH4_TgCH4_2[:,-1]=fCH4_TgCH4_2[:,-2]

scenario_dir='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/'
ch4_n2o_file=scenario_dir+'SSP2-2.6_IMAGE_concs_ch4_n2o.txt'
co2_file=scenario_dir+'SSP2-2.6_IMAGE_concs_co2.txt'
qnonco2_file=scenario_dir+'SSP2-2.6_IMAGE_qnonco2.txt'

ch4_n2o,ch4_years=data_info.read_scenario_file(ch4_n2o_file,ndata=2,return_years=True)
co2_int_ppmv=data_info.read_scenario_file(co2_file)[0]
qnonco2_int=data_info.read_scenario_file(qnonco2_file)[0]

ch4_years=ch4_years[:nYEARs]
ch4_int_ppbv=ch4_n2o[0][:nYEARs]
n2o_int_ppbv=ch4_n2o[1][:nYEARs]

tau_ch4_ref=8.4
ch4_ppbv_ref=1751.02
ref_year=2000
atmos_gain_ch4_init_TgCH4 =180.  #   # TgCH4 per year
TgCH4_to_GtC = (1e9/1e12)*( 12.01/16.04 )
CONV=0.471

print(ch4_int_ppbv.shape)
irefyear=ref_year-data_start_year
atmCH4_new_1=np.zeros_like(atmCO2)+ch4_int_ppbv
atmCH4_new_2=np.zeros_like(atmCO2)+ch4_int_ppbv

for iyear in range(data_end_year-ref_year+1):
    year = iyear+ref_year
    iiyear = iyear + irefyear
    #print(iyear)
    ch4_ppbv_itm = np.copy(atmCH4_new_1[...,iiyear-1])  # previous time-step ch4_ppbv
    #land_gain_ch4 = -fCH4_TgCH4[...,iiyear]    # TgCH4 per year *1e12
    atmos_gain_ch4_TgCH4 = fCH4_TgCH4_1[...,iiyear]  # TgCH4 per year *1e12
    d_land_atmos_ch4 = ( atmos_gain_ch4_TgCH4 - atmos_gain_ch4_init_TgCH4 ) * CONV * 1e3 * TgCH4_to_GtC
    tau_ch4 = tau_ch4_ref * np.exp( 0.28 * np.log(ch4_ppbv_itm/ch4_ppbv_ref) )
    ch4_ppbv = ch4_ppbv_itm + d_land_atmos_ch4              \
            - ((ch4_ppbv_itm - ch4_int_ppbv[iiyear]) / tau_ch4)\
            + (ch4_int_ppbv[iiyear] - ch4_int_ppbv[iiyear-1])
    atmCH4_new_1[...,iiyear]=np.copy(ch4_ppbv)

    ch4_ppbv_itm = np.copy(atmCH4_new_2[...,iiyear-1])  # previous time-step ch4_ppbv
    #land_gain_ch4 = -fCH4_TgCH4[...,iiyear]    # TgCH4 per year *1e12
    atmos_gain_ch4_TgCH4 = fCH4_TgCH4_2[...,iiyear]  # TgCH4 per year *1e12
    d_land_atmos_ch4 = ( atmos_gain_ch4_TgCH4 - atmos_gain_ch4_init_TgCH4 ) * CONV * 1e3 * TgCH4_to_GtC
    tau_ch4 = tau_ch4_ref * np.exp( 0.28 * np.log(ch4_ppbv_itm/ch4_ppbv_ref) )
    ch4_ppbv = ch4_ppbv_itm + d_land_atmos_ch4              \
            - ((ch4_ppbv_itm - ch4_int_ppbv[iiyear]) / tau_ch4)\
            + (ch4_int_ppbv[iiyear] - ch4_int_ppbv[iiyear-1])
    atmCH4_new_2[...,iiyear]=np.copy(ch4_ppbv)

# Etminan:
a1=-2.4e-7
b1=7.2e-4
c1=-2.1e-4
co2_ppm_pi=co2_int_ppmv[0]

delQ_ch4      =delQ.etminan_CH4(ch4_int_ppbv,n2o_int_ppbv[:])
delQ_ch4_new_1=delQ.etminan_CH4(atmCH4_new_1,n2o_int_ppbv[:])
delQ_ch4_new_2=delQ.etminan_CH4(atmCH4_new_2,n2o_int_ppbv[:])

delQ_co2 = delQ.etminan_CO2(atmCO2,n2o_int_ppbv[:],co2_ppm_0=co2_ppm_pi,n2o_ppb_0=n2o_int_ppbv[0])

delQ_co2_new_1 = delQ_co2 + delQ_ch4 - delQ_ch4_new_1
delQ_co2_new_2 = delQ_co2 + delQ_ch4 - delQ_ch4_new_2

atmCO2_new_1     = np.zeros_like(atmCO2)
atmCO2_new_2     = np.zeros_like(atmCO2)

nGCMs = atmCO2.shape[0]
#nSCENARIOs = atmCO2.shape[1]
for igcm in range(nGCMs):
  #for iscen in range(nSCENARIOs)
    for iyear in range(nYEARs):
        #co2_ppmv_new=delQ.etminan_CO2_inverse(delQ_co2_new[...,iyear],
        #                                      atmCO2[...,iyear],n2o_int_ppbv[iyear])
        delQ = delQ_co2_new_1[igcm,iyear]
        co2_iter=999.
        co2=atmCO2[igcm,iyear]
        Nbar=(n2o_int_ppbv[iyear]-n2o_int_ppbv[0])*0.5
        itercnt=0
        while ( (np.abs(co2_iter-co2)>0.001) & (itercnt<=1000) ):
            co2_iter=co2.copy()
            denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     
                     + (b1*(co2_iter-co2_ppm_pi))           
                     + (c1*Nbar)+5.36 )
            co2 = co2_ppm_pi * np.exp( delQ/denom )
            itercnt+=1
        atmCO2_new_1[igcm,iyear]=np.copy(co2)

        delQ = delQ_co2_new_2[igcm,iyear]
        co2_iter = 999.
        co2 = atmCO2[igcm,iyear]
        Nbar = (n2o_int_ppbv[iyear]-n2o_int_ppbv[0])*0.5
        itercnt=0
        while ( (np.abs(co2_iter-co2)>0.001) & (itercnt<=1000) ):
            co2_iter=co2.copy()
            denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     
                     + (b1*(co2_iter-co2_ppm_pi))           
                     + (c1*Nbar)+5.36 )
            co2 = co2_ppm_pi * np.exp( delQ/denom )
            itercnt+=1
        atmCO2_new_2[igcm,iyear]=np.copy(co2)

delta_atmCO2=atmCO2_new_1 - atmCO2_new_2

np.savetxt(data_dir+runid+'delta_atmCO2.dat',delta_atmCO2)
np.savetxt(data_dir+runid+'atmCO2_1.dat',atmCO2_new_1)
np.savetxt(data_dir+runid+'atmCO2_2.dat',atmCO2_new_2)

np.savetxt(data_dir+runid+'atmCH4_1.dat',atmCH4_new_1)
np.savetxt(data_dir+runid+'atmCH4_2.dat',atmCH4_new_2)

#pdb.set_trace()
fig,axes = plt.subplots(ncols=1,nrows=2,figsize=(8,10))

axes[0].plot(ch4_years,ch4_int_ppbv,c='g')
axes[0].plot(ch4_years,atmCH4_new_1.transpose(),c='b')
axes[0].plot(ch4_years,atmCH4_new_2.transpose(),c='r',ls=':')
axes[0].set_ylabel('Atmopsheric CH$_4$ (ppbv)')

axes[1].plot(ch4_years,atmCO2.transpose(),c='g')
axes[1].plot(ch4_years,atmCO2_new_1.transpose(),c='b')
axes[1].plot(ch4_years,atmCO2_new_2.transpose(),c='r',ls=':')
axes[1].set_ylabel('Atmopsheric CO$_2$ (ppmv)')


fig.savefig('atmCH4_atmCO2_profiles.png')
plt.show()


