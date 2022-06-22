#!/bin/env python2.7

import matplotlib.pyplot as plt
import pickle,os
import numpy as np
from imogen import data_info

cmip5_runs=data_info.cmip5_runs()
GCMs=data_info.GCMs()
nGCMs=len(GCMs)

start_year= 1850
end_year  = 2100
out_years = np.arange(start_year,end_year+1)
nyears    = len(out_years)

pickle_dir='/users/global/chg/imogen/build/pickled_yearly/'

fill_value=-9999.
Tprofiles=np.zeros([nGCMs,nyears])+fill_value
#ProfNames=np.zeros([nGCMs,nyears])+fill_value

for igcm in range(nGCMs):
    gcm=GCMs[igcm]
    infile=pickle_dir+gcm+'_historical_rcp26.pickle'

    if not os.path.isfile(infile):
        continue
    
    timeseries_historical_lower = pickle.load( open(infile, 'rb'))
    
    i_ens=0

    # preI backgorund removal
    tas_control_est=np.mean(timeseries_historical_lower[i_ens][4][0:40])
    tas_delta = timeseries_historical_lower[i_ens][4].squeeze() - tas_control_est
    
    igcm_len = len(tas_delta)
    igcm_start = timeseries_historical_lower[i_ens][7]
    igcm_end = igcm_start+igcm_len-1
    igcm_years = np.arange(igcm_start,igcm_end+1)
    nigcm_years = len(igcm_years)

    if igcm_start>=start_year:
        igcm_start_index = 0
        out_start_index  = np.where(out_years==igcm_start)[0][0]
    else:
        igcm_start_index = np.where(igcm_years==start_year)[0][0]
        out_start_index  = 0

    if igcm_end<=end_year:
        igcm_end_index = nigcm_years
        out_end_index  = np.where(out_years==igcm_end)[0][0]+1
    else:
        igcm_end_index = np.where(igcm_years==end_year)[0][0]+1
        out_end_index  = nyears

    
    Tprofiles[igcm,out_start_index:out_end_index] = tas_delta[igcm_start_index:igcm_end_index]

Tprofiles = np.ma.masked_equal(Tprofiles,fill_value)

plt.plot(out_years,Tprofiles.transpose(1,0))
plt.show()

np.savetxt('CMIP5_tprofile_data.dat',Tprofiles)





