#!/bin/env python3

import numpy as np
from imogen import ocean_co2, profile,data_info

out_dir='/users/eow/edwcom/CLIFFTOP/plots/toy_jules/SSP2-2.6_neweqns_1p5deg/'

cmip5_runs = data_info.cmip5_runs()
n_cmip5=len(cmip5_runs)

dt_limit=1.5
mu_zero = 0.1
mu_one = 0.000
beta=0.025                # K/yr
dt_now=0.89
yr_now=2015

nu_all = np.zeros(n_cmip5)
f_all = np.zeros(n_cmip5)

file_nu = open('/users/global/chg/imogen/build/imogen_vals/nu.dat')
i_line=0
for line in file_nu:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        nu_all[i_line] = np.float(in_vals[0])
    else:
        print( 'nu file not lining up OK') ; sys.exit()
    i_line = i_line + 1

file_frac = open('/users/global/chg/imogen/build/imogen_vals/ocean_frac.dat')
i_line=0
for line in file_frac:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        f_all[i_line] = np.float(in_vals[0])
    else:
        print('frac file not lining up OK') ; sys.exit()
    i_line = i_line + 1

co2_ppm_all=np.load(out_dir+'co2_ppm.npy')
print(co2_ppm_all.shape)

delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one)

delta_temp_ocean_yearly_all = np.array([delta_temp_global / \
                                        (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])


d_ocean_atmos_all=ocean_co2.ocean_co2_parallel(co2_ppm_all,delta_temp_ocean_yearly_all)


