#!/bin/env python

import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import parabolic, profile, ocean_co2, data_info

l_outTprof=True    # Output temperature profile?
SCENARIO = '1p81p5deg' #'1p81p5deg'  #'2deg'  #'1p5deg' #
SCENARIOs = ['1p5deg','1p81p5deg','2deg'] #'1p81p5deg'  #'2deg'  #'1p5deg' #
Etminan=True
Baseline='SSP2-2.6_IMAGE'
version='vn2p0'   # For saving output driving files

SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'

beta=0.025                # K/yr
dt_now=0.89
yr_now=2015
end_year=2500
start_year=1850
n_yr=end_year-start_year+1

for SCENARIO in SCENARIOs:
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
        mu_zero = 0.05
        mu_one = 0.00

    delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one)
    delta_temp_global = delta_temp_global[:n_yr]
    yr_plot=np.arange(start_year,end_year+1)

    now_index=np.where(yr_plot==yr_now)[0][0]
    plt.plot(yr_plot,delta_temp_global)
    plt.grid(True)
    plt.show()
#quit()

    # Output Tprofile if required
    if l_outTprof:
        out_Tprof_filename=SCENARIO_DIR+SCENARIO+'_global_temp_anomaly_'+version+'.dat'
        print( out_Tprof_filename)
        outf=open(out_Tprof_filename,'w')
        for year,Temp in zip(yr_plot,delta_temp_global):
            line='%4i  %8.4f\n'%(year,Temp)
            outf.write(line)
        outf.close()

