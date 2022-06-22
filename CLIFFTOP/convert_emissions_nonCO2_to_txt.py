#!/bin/env python3

import numpy as np
from imogen import data_info
import matplotlib.pyplot as plt

scenario='2deg'
outdir='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/C_emissions/'
Emissions_infile='./IMOGEN_MED_CO2_Emissions_Param.npy'
nonCO2_infile='./dq_nonCO2.npy'

start_year,end_year=1850,2100
years=np.arange(start_year,end_year+1)
n_yrs=end_year-start_year+1

cmip5_runs=data_info.GCMs()
n_cmip5=len(cmip5_runs)

C_emissions=np.load(Emissions_infile)
for i_gcm in range(n_cmip5):
    C_emissions[:10,i_gcm] = np.interp(np.arange(10),np.array([0,9]),C_emissions[[0,9],i_gcm])

qnonco2=np.load(nonCO2_infile)
print(qnonco2.shape)

for i_gcm in range(n_cmip5):
    outf=open(outdir+cmip5_runs[i_gcm]+'_'+scenario+'_Cemissions.txt','w')
    for iyr in range(n_yrs):
        outf.write('%4i %10.5f\n' % (years[iyr],C_emissions[iyr,i_gcm]) )
    outf.close()
    
    outf_nonCO2=open(outdir+cmip5_runs[i_gcm]+'_'+scenario+'_qnonco2.txt','w')
    for iyr in range(n_yrs):
        outf_nonCO2.write('%4i %10.5f\n' % (years[iyr],qnonco2[iyr,i_gcm]) )
    outf_nonCO2.close()



