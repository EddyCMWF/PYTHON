# This is the top-level routine for the IMOGEN-EBM inversion, to get first estimate of emissions.
# C. Huntingford (13th March 2017)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys,os
import csv, math
sys.path.insert(0, '/users/global/chg/CLUES/imogen_co2/code/subroutines')
import parabolic, profile, ocean_co2

# First set the 3 parameters in the temperature curves
dt_limit=2.0
mu_zero = 0.05
mu_one = 0.0

# Also normalise the curves so that they end up at the current temperature and gradient estimate
beta=0.025                # K/yr
dt_now=0.89

delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one)
n_yr = len(delta_temp_global)
yr_plot = np.zeros(n_yr)

for i in range(0,n_yr):
    yr_plot[i]=1850.0 + np.float(i)

cmip5_runs = [['BCC','bcc-csm1-1'],
              ['BCC','bcc-csm1-1-m'],
              ['BNU','BNU-ESM'],
              ['CCCma','CanESM2'],
              ['CMCC','CMCC-CMS'],
              ['CNRM-CERFACS','CNRM-CM5'],
              ['CSIRO-BOM','ACCESS1-0'],
              ['CSIRO-BOM','ACCESS1-3'],
              ['CSIRO-QCCCE','CSIRO-Mk3-6-0'],
              ['INM','inmcm4'],
              ['IPSL','IPSL-CM5A-LR'],
              ['IPSL','IPSL-CM5A-MR'],
              ['IPSL','IPSL-CM5B-LR'],
              ['MIROC','MIROC5'],
              ['MIROC','MIROC-ESM'],
              ['MIROC','MIROC-ESM-CHEM'],
              ['MOHC','HadGEM2-CC'],
              ['MOHC','HadGEM2-ES'],
              ['MPI-M','MPI-ESM-LR'],
              ['MPI-M','MPI-ESM-MR'],
              ['MRI','MRI-CGCM3'],
              ['NASA-GISS','GISS-E2-H'],
              ['NASA-GISS','GISS-E2-H-CC'],
              ['NASA-GISS','GISS-E2-R'],
              ['NASA-GISS','GISS-E2-R-CC'],
              ['NCAR','CCSM4'],
              ['NCC','NorESM1-M'],
              ['NCC','NorESM1-ME'],
              ['NOAA-GFDL','GFDL-CM3'],
              ['NOAA-GFDL','GFDL-ESM2G'],
              ['NOAA-GFDL','GFDL-ESM2M'],
              ['NSF-DOE-NCAR','CESM1-BGC'],
              ['NSF-DOE-NCAR','CESM1-CAM5'],
              ['NSF-DOE-NCAR','CESM1-WACCM']]
n_cmip5 = len(cmip5_runs)
kappa_all=np.zeros(n_cmip5) ; lambda_l_all=np.zeros(n_cmip5) ; lambda_o_all=np.zeros(n_cmip5)
nu_all = np.zeros(n_cmip5) ; f_all = np.zeros(n_cmip5)

# Read in the EBM parameters
file_kappa = open('/users/global/chg/imogen/build/imogen_vals/kappa.dat')
i_line=0
for line in file_kappa:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        kappa_all[i_line] = np.float(in_vals[0])
    else:
        print 'kappa file not lining up OK' ; sys.exit()
    i_line = i_line + 1 

file_lambda_o = open('/users/global/chg/imogen/build/imogen_vals/lambda_o.dat')
i_line=0
for line in file_lambda_o:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        lambda_o_all[i_line] = np.float(in_vals[0])
    else:
        print 'lambda_o file not lining up OK' ; sys.exit()
    i_line = i_line + 1 

file_lambda_l = open('/users/global/chg/imogen/build/imogen_vals/lambda_l.dat')
i_line=0
for line in file_lambda_l:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        lambda_l_all[i_line] = np.float(in_vals[0])
    else:
        print 'lambda_l file not lining up OK' ; sys.exit()
    i_line = i_line + 1 

file_nu = open('/users/global/chg/imogen/build/imogen_vals/nu.dat')
i_line=0
for line in file_nu:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        nu_all[i_line] = np.float(in_vals[0])
    else:
        print 'nu file not lining up OK' ; sys.exit()
    i_line = i_line + 1 

file_frac = open('/users/global/chg/imogen/build/imogen_vals/ocean_frac.dat')
i_line=0
for line in file_frac:
    in_vals = line.split()
    if cmip5_runs[i_line][0] == in_vals[1] and cmip5_runs[i_line][1] == in_vals[2] :
        f_all[i_line] = np.float(in_vals[0])
    else:
        print 'frac file not lining up OK' ; sys.exit()
    i_line = i_line + 1 

# Loop over the different GCMs emulated
for i_gcm in range(0, n_cmip5):
#for i_gcm in range(0, 3):
    delta_temp_ocean_yearly = delta_temp_global / (f_all[i_gcm] + (1.0-f_all[i_gcm])*nu_all[i_gcm])
    temp_gradient_out = parabolic.parabolic(kappa_all[i_gcm], delta_temp_ocean_yearly)

    n_yr = temp_gradient_out.shape[0]

    if i_gcm == 0:
        temp_gradient_out_all=np.zeros([n_yr,n_cmip5])                  # Holds all CO2 concentrations 
    temp_gradient_out_all[:,i_gcm]=temp_gradient_out

    dq = np.zeros(n_yr)
    # Derive time-evolution of radiative forcing, Q
    factor = (((1.0 - f_all[i_gcm])*lambda_l_all[i_gcm]*nu_all[i_gcm])/f_all[i_gcm]) + lambda_o_all[i_gcm]
    dq = (temp_gradient_out + delta_temp_ocean_yearly*factor)*f_all[i_gcm]
    if i_gcm == 0:
        dq_all=np.zeros([n_yr,n_cmip5])                  # Holds all CO2 concentrations 
    dq_all[:,i_gcm]=dq

    # Read in the non-CO2 RF from PIK site
    if i_gcm == 0:
        dq_non_co2 = np.zeros(n_yr)
        f_in = open('RCP3PD_MIDYEAR_RADFORCING_mpeg.txt', 'r')

        # Ignore headers
        for i in range(0,59):
            header=f_in.readline()
        # Ignore years 1765 to 1849
        for i in range(0,85):
            header=f_in.readline()

        for iyr in range(0, 2500-1850+1): 
            nums_in_string=f_in.readline()
            dq_non_co2[iyr]= np.float(nums_in_string.split()[1]) - np.float(nums_in_string.split()[8]) 
        f_in.close()

        for iyr in range(651, n_yr):
            dq_non_co2[iyr]=dq_non_co2[650]


        dq_smooth=True
        if dq_smooth:
            dq_non_co2_temp=np.copy(dq_non_co2)
            mean_len_dq = 31
            n_offset_dq=np.int(math.floor(np.float(mean_len_dq)/2.0))

            for i_yr in range(0, n_yr-mean_len_dq+1):
                dq_non_co2_temp[i_yr+n_offset_dq]=np.mean(dq_non_co2[i_yr: i_yr+mean_len_dq])
            for i_yr in range(0,n_offset_dq):
                dq_non_co2_temp[i_yr] = dq_non_co2_temp[n_offset_dq]
            for i_yr in range(0,n_offset_dq):
                dq_non_co2_temp[n_yr-n_offset_dq+i_yr] = dq_non_co2_temp[n_yr-n_offset_dq-1]

            dq_non_co2=np.copy(dq_non_co2_temp)

    # Now back out CO2 concentrations, 
    dq_co2 = dq - dq_non_co2
    q2co2 = 3.74
    co2_ppm_pi = 286.085

    co2_ppm = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2)/q2co2)
    if i_gcm == 0:
        co2_ppm_all=np.zeros([n_yr,n_cmip5])                  # Holds all CO2 concentrations 
    co2_ppm_all[:,i_gcm]=co2_ppm

    # Now calculate the components of the global carbon cycle
    # First calculate amount that has gone in to the oceans.
    fa_ocean=np.zeros(20000) ; nfarray = 20000 ; t_ocean_init = 289.28
    ocean_area = 3.627E14                          # (m2)
    conv = 0.471                                   # From ppm to GtC

    if i_gcm == 0:
        emiss_invert=np.zeros([n_yr,n_cmip5])                  # Holds the final backed-out emissions (GtC/yr)

    for i_yr in range(0, n_yr):
        if i_yr ==0:
            co2_atmos_change_ppmv = 0.0
        else:
            co2_atmos_change_ppmv = co2_ppm[i_yr]-co2_ppm[i_yr-1]
#    d_ocean_atmos = ocean_co2.ocean_co2(j, co2_atmos_ppmv, co2_atmos_ppmv_init, t_ocean_new[0],
#                 fa_ocean,ocean_area,co2_atmos_change_ppmv,nyr,t_ocean_init,nfarray)
        d_ocean_atmos = ocean_co2.ocean_co2(i_yr, co2_ppm[i_yr], co2_ppm_pi, delta_temp_ocean_yearly[i_yr], 
                          fa_ocean,ocean_area,co2_atmos_change_ppmv,n_yr,t_ocean_init,nfarray)

        d_ocean_atmos_gtc = d_ocean_atmos/conv                  # From ppm change back to GtC in to the ocean

        # Close the carbon cycle. [dCO2 atmos = E - 0.25E + d_ocean_atmos] (all in GtC units. 0.25 E assumed in to land) 
        co2_atmos_change_gtc = co2_atmos_change_ppmv/conv

        emiss= (co2_atmos_change_gtc - d_ocean_atmos_gtc)/0.75
        emiss_invert[i_yr, i_gcm]=emiss

    # Get running means of emissions.
    mean_len = 31
    n_offset=np.int(math.floor(np.float(mean_len)/2.0))
    if i_gcm == 0:
        emiss_invert_running=np.zeros([n_yr-mean_len+1,n_cmip5])           # Running mean backed-out emissions (GtC/yr)
    yr_plot_running=np.zeros(n_yr-mean_len+1)

    for i_yr in range(0, n_yr-mean_len+1):
        yr_plot_running[i_yr] = yr_plot[i_yr+n_offset]
        emiss_invert_running[i_yr, i_gcm]=np.mean(emiss_invert[i_yr: i_yr+mean_len, i_gcm])

np.save('output/emiss_invert_running.npy', emiss_invert_running)
np.save('output/delta_temp_global.npy', delta_temp_global)
np.save('output/dq_all.npy', dq_all)
np.save('output/co2_ppm_all.npy', co2_ppm_all)

sys.exit()
