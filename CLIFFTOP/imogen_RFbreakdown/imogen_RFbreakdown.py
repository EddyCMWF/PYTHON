#!/bin/env python

import ipdb
import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import parabolic, profile, ocean_co2, data_info, delQ, imogen_ebm
from imogen import CMIP_CCFBs as CCFBs
#CCFB = CCFBs.CarbonClimateFeedbacks()

# Directories:
SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/toy_imogen_experiment/'
out_dir='/prj/CLIFFTOP/IMOGEN_RFbeakdown/plots/toy_imogen/'

# Configuration options:
SCENARIO ='1p5deg'  # '2deg'  #'1p5deg' #'1p81p5deg'
Etminan=True
Baseline='IMAGE-SSP2-26'   #'RCP2.6'
version='v1.0'   # For saving output driving files

ocean_model = 'Joos'
land_model = 'CCFB'
ccfb_version = 'CMIP5'

# RF treatement options:
l_Offset=False
l_ScaleFactor=False
l_OffsetScaleFactor=False
l_DecayOffset=True
l_Tprof_Adjust = False
l_Tprof_Adjust_DOS = False

l_outTprof=False   # Output temperature profile?

# Conversion factor for GtC to CO2 ppm
GtC_to_ppm=0.471
q2co2=3.74                # Old CO2 RF parameter
epsilon_stable = 0.01  # value of epsilon at end_year (W/m2), used to calculate gamma
gamma = 0.01   # efolding rate for TprfAdDOC

# Options for time period and for calculating the temperature pathways.
beta=0.025                # K/yr at present day
dt_now=0.89
yr_now=2015
yr_beta=1995              # year to apply the beta gradient
end_year=2100


# User defined functions:
def co2ppm_to_datmosGtC(co2_ppm):
    delta_co2_ppm     = np.zeros_like(co2_ppm)
    # yearly change in atmospheric co2
    delta_co2_ppm[:-1] = co2_ppm[1:,...]-co2_ppm[:-1,...]
    delta_co2_ppm[-1]  = delta_co2_ppm[-2]
    d_atmos_GtC        = delta_co2_ppm/GtC_to_ppm
    return d_atmos_GtC 

def co2nT_to_deltaOceanGtC(co2_ppm, d_atmos_GtC, 
                           delta_temp_global, delta_temp_ocean,
                           ocean_model='Joos', ccfb_version='CMIP5' ):
    if ocean_model == 'Joos':
        # Calculate Ocean Uptake as a function of T and co2 history (Joos 1996)
        d_ocean_GtC = ocean_co2.ocean_co2_parallel(co2_ppm,
                                           delta_temp_ocean) / GtC_to_ppm * -1.0
        d_ocean_GtC = d_ocean_GtC[...,np.newaxis]
        ####### NOT USING CMIP OCEAN MODELS AS THEY ARE LESS PHYSICALLY BASED THAN JOOS 
    elif ocean_model == 'CCFB':
        ccfb_model= CCFBs.CarbonClimateFeedbacks( co2_ppm-co2_ppm[0,:],
                                               delta_temp_global, cmip=ccfb_version )
        ccfb_models = ccfb_model.params.models
        n_ccfb_models =  ccfb_model.params.nmodels
        Cuptake_ocean  = ccfb_model.delCocean_Accum_Ensemble()
        d_ocean_GtC_list = [ np.zeros_like(co2_ppm) for imod in range(n_ccfb_models) ]
        for imod in range(n_ccfb_models):
            model = ccfb_models[imod]
            d_ocean_GtC_list[imod][:-1,:] = Cuptake_ocean[model][1:,:] - Cuptake_ocean[model][:-1,:]
            # Replicate the last value
            d_ocean_GtC_list[imod][-1,:]  = d_ocean_GtC_list[imod][-2,:]
        d_ocean_GtC = np.array(d_ocean_GtC_list).transpose(1,2,0)
    else:
        print('Unrecognised ocean model, setting to half of atmospheric increase')
        d_ocean_GtC= d_atmos_GtC/2.0
        d_ocean_GtC= d_ocean_GtC[...,np.newaxis]
    return d_ocean_GtC

def co2nT_to_deltaLandGtC(co2_ppm, d_atmos_GtC, d_ocean_GtC, delta_temp_global, 
                            land_model='CCFB', ccfb_version='CMIP5' ): 
    # Land uptake and CO2 emissions together: 
    if land_model == 'CCFB':
        ccfb_model= CCFBs.CarbonClimateFeedbacks( co2_ppm-co2_ppm[0,:],
                                               delta_temp_global, cmip=ccfb_version )
        Cuptake_land   = ccfb_model.delCland_Accum_Ensemble()
        ccfb_models   = ccfb_model.params.models
        n_ccfb_models =  ccfb_model.params.nmodels
        d_land_GtC_list = [ np.zeros_like(co2_ppm) for i in range(n_ccfb_models) ]
        for imod in range(n_ccfb_models):
            model = ccfb_models[imod]
            d_land_GtC_list[imod][:-1,:] = Cuptake_land[model][1:,:] - Cuptake_land[model][:-1,:]
            # Replicate the last value
            d_land_GtC_list[imod][-1,:]  = d_land_GtC_list[imod][-2,:]
        d_land_GtC = np.array(d_land_GtC_list).transpose(1,2,0)
    else:
        # 25% uptake by land = 1/3 of Atmos+Ocean:
        d_land_GtC = (d_atmos_GtC+d_ocean_GtC)/3.0
        d_land_GtC = d_land_GtC[...,np.newaxis]

    return d_land_GtC


# Get CMIP5 run names
cmip5_runs = data_info.cmip5_runs()
n_cmip5 = len(cmip5_runs)
# Read in the EBM parameters
##imogen_ebm = data_info.imogen_ebm_params(cmip5_runs)  # Easier to read one by one
kappa_all=data_info.kappa(cmip5_runs)
lambda_l_all=data_info.lambda_l(cmip5_runs)
lambda_o_all=data_info.lambda_o(cmip5_runs)
nu_all=data_info.nu(cmip5_runs)
f_all=data_info.ocean_frac(cmip5_runs)

# Dictionary to store scenario options
SCENARIOS={'1p5deg':   {'dt_limit':1.5,'mu_zero':0.08,'mu_one':0.0,
                        'color':'#6495ED','label':'1.5$^o$C','ls':'-',
                        'SSP_TAG':'IMAGE-SSP2-26','RCP_TAG':'RCP2.6'},
           '1p81p5deg':{'dt_limit':1.5,'mu_zero':-0.01,'mu_one':0.00087,
                        'color':'#EDBC64','label':'1.5$^o$C (overshoot)','ls':'--', 
                        'SSP_TAG':'IMAGE-SSP2-26','RCP_TAG':'RCP2.6'},
           '2deg':     {'dt_limit':2.0,'mu_zero':0.08,'mu_one':0.0,
                        'color':'#ED6495','label':'2.0$^o$C','ls':':', 
                        'SSP_TAG':'IMAGE-SSP2-26','RCP_TAG':'RCP2.6'},
           '3deg':     {'dt_limit':3.0,'mu_zero':0.065,'mu_one':0.0,
                        'color':'#6F36C8','label':'3.0$^o$C','ls':'-.', 
                        'SSP_TAG':'IMAGE-SSP2-26', 'RCP_TAG':'RCP2.6'},
             }
dt_limit=SCENARIOS[SCENARIO]['dt_limit']
mu_zero =SCENARIOS[SCENARIO]['mu_zero']
mu_one  =SCENARIOS[SCENARIO]['mu_one']
SSP_TAG =SCENARIOS[SCENARIO]['SSP_TAG']
RCP_TAG =SCENARIOS[SCENARIO]['RCP_TAG']

# Set up RF equation dependent filenames and plot_tag
if Etminan:
    RAD_FORCING_FILE = SCENARIO_DIR + Baseline + '_qnonco2.txt'
    PLOT_TAG = Baseline+'_neweqns_'+SCENARIO
    if 'RCP' in Baseline.upper(): 
        print('Warning, RCP data compputed with old RF equations, update file before continuing')
        quit()
else:
    RAD_FORCING_FILE=SCENARIO_DIR+Baseline+'_qnonco2_ipccRF.txt'
    PLOT_TAG=Baseline+'_oldeqns_'+SCENARIO

# Set up standard start_years and time indexes
if 'SSP2' in PLOT_TAG:
    start_year=1850
elif 'RCP' in PLOT_TAG:
    start_year=1859
n_yr=end_year-start_year+1
yr_plot=np.arange(start_year,end_year+1)
now_index=np.where(yr_plot==yr_now)[0][0]

print(PLOT_TAG)
out_dir = out_dir+PLOT_TAG+'/'
print(out_dir)
os.system('mkdir -p '+out_dir)

# Read in the non-CO2 RF
if os.path.isfile(out_dir+'dq_non_co2_ssp.npy'):
    dq_non_co2_ssp = np.load(out_dir+'dq_non_co2_ssp.npy')
else:
    dq_non_co2_ssp = np.zeros(n_yr)
    RF_lines=open(RAD_FORCING_FILE, 'r').readlines()
    for iyr in range(n_yr):
        split=RF_lines[iyr].split()
        dq_non_co2_ssp[iyr]=float(split[1])
    np.save(out_dir+'dq_non_co2_ssp.npy', dq_non_co2_ssp)

# Read in the CH4 and N2O concentrations
if os.path.isfile(out_dir+'ch4_ppb_ssp.npy'):
    ch4_ppb_ssp=np.load(out_dir+'ch4_ppb_ssp.npy')
    n2o_ppb_ssp=np.load(out_dir+'n2o_ppb_ssp.npy')
else:
    CH4_N2O_FILE = SCENARIO_DIR + Baseline + '_concs_ch4_n2o.txt'
    CH4_N2O_lines=open(CH4_N2O_FILE,'r').readlines()
    ch4_ppb_ssp=np.zeros(n_yr)
    n2o_ppb_ssp=np.zeros(n_yr)
    for iyr in range(n_yr):
        split=CH4_N2O_lines[iyr].split()
        ch4_ppb_ssp[iyr]=split[1]
        n2o_ppb_ssp[iyr]=split[2]
    np.save(out_dir+'ch4_ppb_ssp.npy', ch4_ppb_ssp)
    np.save(out_dir+'n2o_ppb_ssp.npy', n2o_ppb_ssp)

# Read in the CO2 concentration and calculate dq_CO2
if os.path.isfile(out_dir+'co2_ppm_ssp.npy'):
    co2_ppm_ssp = np.load(out_dir+'co2_ppm_ssp.npy')
    dq_co2_ssp = np.load(out_dir+'dq_co2_ssp.npy')
else:
    CO2_FILE     = SCENARIO_DIR + Baseline + '_concs_co2.txt'
    CO2_lines=open(CO2_FILE,'r').readlines()
    co2_ppm_ssp=np.zeros(n_yr)
    for iyr in range(n_yr):
        split=CO2_lines[iyr].split()
        co2_ppm_ssp[iyr]=float(split[1])
    # Calculate corresponding radiative forcing contribution
    if Etminan:
        dq_co2_ssp = delQ.etminan_CO2( co2_ppm_ssp, n2o_ppb_ssp,
                                  co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
    else:
        dq_co2_ssp = np.log(co2_ppm_ssp/co2_ppm_ssp[0]) * (q2co2/np.log(2.))
    np.save(out_dir+'co2_ppm_ssp.npy', co2_ppm_ssp)
    np.save(out_dir+'dq_co2_ssp.npy', dq_co2_ssp)


if os.path.isfile(out_dir+'delta_temp_global.npy'):
    delta_temp_global = np.load(out_dir+'delta_temp_global.npy')
else:
    delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one, yr_now=yr_now, yr_beta=yr_beta)
    delta_temp_global = delta_temp_global[:n_yr]
    np.save(out_dir+'delta_temp_global.npy', delta_temp_global)

# Output Tprofile if required
if l_outTprof:
    out_Tprof_filename=out_dir+SCENARIO+'_global_temp_anomaly_'+version+'.dat'
    outf=open(out_Tprof_filename,'w')
    for year,Temp in zip(yr_plot,delta_temp_global):
        line='%4i  %8.4f\n'%(year,Temp)
        outf.write(line)
    outf.close()

# Save delta_temp_global per gcm for later functionality
delta_temp_global_all = np.array([delta_temp_global for i_gcm in range(n_cmip5)]).transpose()
delta_temp_ocean_all = delta_temp_global_all / (f_all + (1.0-f_all)*nu_all)
np.save(out_dir+'delta_temp_ocean_all.npy',delta_temp_ocean_all)
delta_temp_land_all = delta_temp_ocean_all * nu_all
np.save(out_dir+'delta_temp_land_all.npy',delta_temp_land_all)

# Loop over the different GCMs emulated and calculate RF via dtemp_o
temp_gradient_out_all=np.zeros([n_yr,n_cmip5]) 
if os.path.isfile(out_dir+'dq_all.npy'):
    dq_all=np.load(out_dir+'dq_all.npy')
else:
    dq_all=np.zeros([n_yr,n_cmip5]) 
    # Have to transpose the arrays in and out of thie function:
    temp_gradient_out_all = parabolic.parabolic_parallel(kappa_all, 
                            delta_temp_ocean_all.transpose()).transpose()
    dq_factor = ( (1.0-f_all)*lambda_l_all*nu_all/f_all ) + lambda_o_all
    dq_all = (temp_gradient_out_all + (delta_temp_ocean_all*dq_factor)) * f_all
    np.save(out_dir+'dq_all.npy',dq_all)

# Dictionary to store the data for plotting:
METHOD_DICT = {}

if l_Offset:
#  Constant Offset Calculations
    if os.path.isfile(out_dir+'co2_ppm_Off.npy'):
        dq_co2_all_Off     = np.load(out_dir+'dq_co2_Off.npy')
        dq_non_co2_all_Off = np.load(out_dir+'dq_non_co2_Off.npy')
        epsilon_all_Off    = np.load(out_dir+'epsilon_Off.npy')
        co2_ppm_all_Off    = np.load(out_dir+'co2_ppm_Off.npy')
    else:
        print('Calculating Emissions for constant offset method')
        #array for all temp gradients and rfs
        co2_ppm_all_Off=np.zeros([n_yr,n_cmip5]) 
        
        # Start with dq_non_co2 as the residual between dq and dq_co2_ssp for the whole time-series
        dq_non_co2_all_Off = dq_all - dq_co2_ssp[...,np.newaxis]
        
        # Calculate the osffet at the now index year
        epsilon_all_Off = dq_non_co2_all_Off[now_index,:]-dq_non_co2_ssp[now_index]

        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_Off[now_index:,:] = dq_non_co2_ssp[now_index:,np.newaxis] + epsilon_all_Off
            
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_Off = dq_all - dq_non_co2_all_Off

        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_Off = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_Off[:,igcm],
                                   n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                           for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_Off = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_Off)/q2co2)
                
        np.save(out_dir+'dq_co2_Off.npy',dq_co2_all_Off)
        np.save(out_dir+'dq_non_co2_Off.npy',dq_non_co2_all_Off)
        np.save(out_dir+'epsilon_Off.npy',epsilon_all_Off)
        np.save(out_dir+'co2_ppm_Off.npy',co2_ppm_all_Off)

    if os.path.isfile(out_dir+'d_atmos_GtC_Off.npy'):
        d_land_GtC_all_Off        = np.load(out_dir+'d_land_GtC_Off.npy')
        d_ocean_GtC_all_Off       = np.load(out_dir+'d_ocean_GtC_Off.npy')
        co2_emissions_GtC_all_Off = np.load(out_dir+'co2_emissions_GtC_Off.npy')
        d_atmos_GtC_all_Off       = np.load(out_dir+'d_atmos_GtC_Off.npy')
    else:
        print('Calculating Emissions for constant offset method')
        d_atmos_GtC_all_Off = co2ppm_to_datmosGtC(co2_ppm_all_Off)[...,np.newaxis]
         
        d_ocean_GtC_all_Off = co2nT_to_deltaOceanGtC(co2_ppm_all_Off, d_atmos_GtC_all_Off,
                                                     delta_temp_global_all,delta_temp_ocean_all,
                                                     ocean_model=ocean_model)
        
        d_land_GtC_all_Off = co2nT_to_deltaLandGtC(co2_ppm_all_Off, d_atmos_GtC_all_Off, d_ocean_GtC_all_Off, 
                                                    delta_temp_global_all, land_model=land_model )
        
        co2_emissions_GtC_all_Off = d_land_GtC_all_Off + d_atmos_GtC_all_Off + d_ocean_GtC_all_Off
        
        np.save(out_dir+'d_ocean_GtC_Off.npy',d_ocean_GtC_all_Off)
        np.save(out_dir+'d_land_GtC_Off.npy',d_land_GtC_all_Off)
        np.save(out_dir+'co2_emissions_GtC_Off.npy',co2_emissions_GtC_all_Off)
        np.save(out_dir+'d_atmos_GtC_Off.npy',d_atmos_GtC_all_Off)
        
    METHOD_DICT['Offset']= {'dtemp_global':delta_temp_global_all,
                            'dq':dq_all,
                            'dq_co2':dq_co2_all_Off,
                            'dq_non_co2':dq_non_co2_all_Off,
                            'co2_ppm':co2_ppm_all_Off,
                            'd_co2_GtC':d_atmos_GtC_all_Off,
                            'co2_emissions':co2_emissions_GtC_all_Off,
                            'ocean_uptake':d_ocean_GtC_all_Off,
                            'land_uptake':d_land_GtC_all_Off,
                            'cum_colour':'#a6cee3',    # Shades of blue for offset
                            'colour':'#1f78b4',
                            }



if l_ScaleFactor:
#  Constant Offset Calculations
    if os.path.isfile(out_dir+'co2_ppm_SF.npy'):
        dq_co2_all_SF       = np.load(out_dir+'dq_co2_SF.npy')
        dq_non_co2_all_SF   = np.load(out_dir+'dq_non_co2_SF.npy')
        epsilon_all_SF      = np.load(out_dir+'epsilon_SF.npy')
        co2_ppm_all_SF      = np.load(out_dir+'co2_ppm_SF.npy')
    else:
        print('Calculating Emissions for constant scale factor method')
        #array for all temp gradients and rfs
        co2_ppm_all_SF=np.zeros([n_yr,n_cmip5]) 
        
        # Start with dq_non_co2 as the residual between dq and dq_co2_ssp for the whole time-series
        dq_non_co2_all_SF = dq_all - dq_co2_ssp[...,np.newaxis]
        
        # Calculate the epsilon scale factor at the now index year
        # epsilon = (dq_all - dq_co2_ssp)/dq_non_co2_ssp   - 1 = (dq_non_co2_all_SF/dq_non_co2_ssp) - 1
        epsilon_all_SF = (dq_non_co2_all_SF[now_index,:]/dq_non_co2_ssp[now_index])-1

        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_SF[now_index:,:] = dq_non_co2_ssp[now_index:,np.newaxis] * (epsilon_all_SF+1.)
            
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_SF = dq_all - dq_non_co2_all_SF

        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_SF = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_SF[:,igcm],
                                   n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                           for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_SF = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_SF)/q2co2)
                
        np.save(out_dir+'dq_co2_SF.npy',dq_co2_all_SF)
        np.save(out_dir+'dq_non_co2_SF.npy',dq_non_co2_all_SF)
        np.save(out_dir+'epsilon_SF.npy',epsilon_all_SF)
        np.save(out_dir+'co2_ppm_SF.npy',co2_ppm_all_SF)

    if os.path.isfile(out_dir+'d_atmos_GtC_SF.npy'):
        d_land_GtC_all_SF        = np.load(out_dir+'d_land_GtC_SF.npy')
        d_ocean_GtC_all_SF       = np.load(out_dir+'d_ocean_GtC_SF.npy')
        co2_emissions_GtC_all_SF = np.load(out_dir+'co2_emissions_GtC_SF.npy')
        d_atmos_GtC_all_SF       = np.load(out_dir+'d_atmos_GtC_SF.npy')
    else:
#if True:
        print('Calculating Emissions for constant scale factor method')
        d_atmos_GtC_all_SF = co2ppm_to_datmosGtC(co2_ppm_all_SF)[...,np.newaxis]
    
        d_ocean_GtC_all_SF = co2nT_to_deltaOceanGtC(co2_ppm_all_SF, d_atmos_GtC_all_SF,
                                                     delta_temp_global_all,delta_temp_ocean_all,
                                                     ocean_model=ocean_model)

        d_land_GtC_all_SF = co2nT_to_deltaLandGtC(co2_ppm_all_SF, d_atmos_GtC_all_SF, d_ocean_GtC_all_SF, 
                                                    delta_temp_global_all, land_model=land_model )

        co2_emissions_GtC_all_SF = d_land_GtC_all_SF + d_atmos_GtC_all_SF + d_ocean_GtC_all_SF
        
        np.save(out_dir+'d_ocean_GtC_SF.npy',d_ocean_GtC_all_SF)
        np.save(out_dir+'d_land_GtC_SF.npy',d_land_GtC_all_SF)
        np.save(out_dir+'co2_emissions_GtC_SF.npy',co2_emissions_GtC_all_SF)
        np.save(out_dir+'d_atmos_GtC_SF.npy',d_atmos_GtC_all_SF)

    # Store data in plotting dictionary:
    METHOD_DICT['ScaleFactor']= {'dtemp_global':delta_temp_global_all,
                                 'dq':dq_all,
                                 'dq_co2':dq_co2_all_SF,
                                 'dq_non_co2':dq_non_co2_all_SF,
                                 'co2_ppm':co2_ppm_all_SF,
                                 'd_co2_GtC':d_atmos_GtC_all_SF,
                                 'co2_emissions':co2_emissions_GtC_all_SF,
                                 'ocean_uptake':d_ocean_GtC_all_SF,
                                 'land_uptake':d_land_GtC_all_SF,
                                 'cum_colour':'#fb9a99',    #  Shades of Orange for scale factor
                                 'colour':'#e31a1c',
                                }
#ipdb.set_trace()
#  Constant Offset and Scale Factor Calculations
#   eps = eps_0 + eps_1*delQ_nonCO2
if l_OffsetScaleFactor:
    if os.path.isfile(out_dir+'co2_ppm_OffSF.npy'):
        dq_co2_all_OffSF     = np.load(out_dir+'dq_co2_OffSF.npy')
        dq_non_co2_all_OffSF = np.load(out_dir+'dq_non_co2_OffSF.npy')
        epsilon_0_all_OffSF  = np.load(out_dir+'epsilon_0_OffSF.npy')
        epsilon_1_all_OffSF  = np.load(out_dir+'epsilon_1_OffSF.npy')
        co2_ppm_all_OffSF    = np.load(out_dir+'co2_ppm_OffSF.npy')
    else:
        print('Calculating Emissions for constant offset and scale factor method')
        # Start with dq_non_co2 as the residual between dq and dq_co2_ssp for the whole time-series
        dq_non_co2_all_OffSF = dq_all - dq_co2_ssp[...,np.newaxis]
        
        # Calculate the gradient od dq, dq_CO2 and dq_nonCO2 at the now index year
        dDelq_dt       = dq_all[now_index,:] - dq_all[now_index-1,:] 
        dDelqCO2_dt    = (dq_co2_ssp[now_index] - dq_co2_ssp[now_index-1])[...,np.newaxis]
        dDelqNonCO2_dt = dq_non_co2_all_OffSF[now_index,:] - dq_non_co2_all_OffSF[now_index-1,:] 
        dDelqNonCO2_dt_ssp = (dq_non_co2_ssp[now_index] - dq_non_co2_ssp[now_index-1])[...,np.newaxis]
        
        # calculate epsilon_1 = [ ( dDelQ_T/dt - dDelQCO2_ssp/dt) / dDelQnonCO2_ssp/dt ] - 1
        #                     = [dDelqNonCO2_dt/dDelqNonCO2_dt_ssp] - 1
        #epsilon_1_all_OffSF = ((dDelq_dt-dDelqCO2_dt)/dDelqNonCO2_dt_ssp) - 1
        epsilon_1_all_OffSF = (dDelqNonCO2_dt/dDelqNonCO2_dt_ssp) - 1

        # calculate espilon_0 = DelQ_T - DelQCO2_ssp - (1+epsilon_1)*delQnonCO2_ssp
        epsilon_0_all_OffSF = dq_all[now_index,:] - dq_co2_ssp[now_index][...,np.newaxis]  \
                - ( (1+epsilon_1_all_OffSF) * dq_non_co2_ssp[now_index,np.newaxis])

        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_OffSF[now_index:,:] =  dq_non_co2_ssp[now_index:,np.newaxis] * (epsilon_1_all_OffSF+1.) \
                                            + epsilon_0_all_OffSF
            
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_OffSF = dq_all - dq_non_co2_all_OffSF

        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_OffSF = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_OffSF[:,igcm],
                                   n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                           for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_OffSF = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_OffSF)/q2co2)
                
        np.save(out_dir+'dq_co2_OffSF.npy',dq_co2_all_OffSF)
        np.save(out_dir+'dq_non_co2_OffSF.npy',dq_non_co2_all_OffSF)
        np.save(out_dir+'epsilon_0_OffSF.npy',epsilon_0_all_OffSF)
        np.save(out_dir+'epsilon_1_OffSF.npy',epsilon_1_all_OffSF)
        np.save(out_dir+'co2_ppm_OffSF.npy',co2_ppm_all_OffSF)

    if os.path.isfile(out_dir+'d_atmos_GtC_OffSF.npy'):
        d_land_GtC_all_OffSF   = np.load(out_dir+'d_land_GtC_OffSF.npy')
        d_ocean_GtC_all_OffSF  = np.load(out_dir+'d_ocean_GtC_OffSF.npy')
        co2_emissions_GtC_all_OffSF  = np.load(out_dir+'co2_emissions_GtC_OffSF.npy')
        d_atmos_GtC_all_OffSF  = np.load(out_dir+'d_atmos_GtC_OffSF.npy')
    else:
#if True:
        print('Calculating Emissions for constant offset and scale factor method')
        d_atmos_GtC_all_OffSF = co2ppm_to_datmosGtC(co2_ppm_all_OffSF)[...,np.newaxis]
    
        d_ocean_GtC_all_OffSF = co2nT_to_deltaOceanGtC(co2_ppm_all_OffSF, d_atmos_GtC_all_OffSF,
                                                     delta_temp_global_all,delta_temp_ocean_all,
                                                     ocean_model=ocean_model)

        d_land_GtC_all_OffSF = co2nT_to_deltaLandGtC(co2_ppm_all_OffSF, d_atmos_GtC_all_OffSF, d_ocean_GtC_all_OffSF, 
                                                    delta_temp_global_all, land_model=land_model )

        co2_emissions_GtC_all_OffSF = d_land_GtC_all_OffSF + d_atmos_GtC_all_OffSF + d_ocean_GtC_all_OffSF
        
        np.save(out_dir+'d_ocean_GtC_OffSF.npy',d_ocean_GtC_all_OffSF)
        np.save(out_dir+'d_land_GtC_OffSF.npy',d_land_GtC_all_OffSF)
        np.save(out_dir+'co2_emissions_GtC_OffSF.npy',co2_emissions_GtC_all_OffSF)
        np.save(out_dir+'d_atmos_GtC_OffSF.npy',d_atmos_GtC_all_OffSF)

    METHOD_DICT['OffsetScaleFactor']= {'dtemp_global':delta_temp_global_all,
                                       'dq':dq_all,
                                       'dq_co2':dq_co2_all_OffSF,
                                       'dq_non_co2':dq_non_co2_all_OffSF,
                                       'co2_ppm':co2_ppm_all_OffSF,
                                       'd_co2_GtC':d_atmos_GtC_all_OffSF,
                                       'co2_emissions':co2_emissions_GtC_all_OffSF,
                                       'ocean_uptake':d_ocean_GtC_all_OffSF,
                                       'land_uptake':d_land_GtC_all_OffSF,
                                       'cum_colour':'#b2df8a',    # Shades of Green for Offset + Scale Factor
                                       'colour':'#33a02c',
                                       }
#ipdb.set_trace()
# Time decaying offset term and one over t term
#  eps =  eps_0 * exp[ eps_1*t - gamma*t^2 ]
if l_DecayOffset:
    if os.path.isfile(out_dir+'co2_ppm_DecOff.npy'):
        dq_co2_all_DecOff     = np.load(out_dir+'dq_co2_DecOff.npy')
        dq_non_co2_all_DecOff = np.load(out_dir+'dq_non_co2_DecOff.npy')
        epsilon_0_all_DecOff  = np.load(out_dir+'epsilon_0_DecOff.npy')
        epsilon_1_all_DecOff  = np.load(out_dir+'epsilon_1_DecOff.npy')
        co2_ppm_all_DecOff    = np.load(out_dir+'co2_ppm_DecOff.npy')
    else:
        print('Calculating Emissions for decaying offset method')
        # Start with dq_non_co2 as the residual between dq and dq_co2_ssp for the whole time-series
        dq_non_co2_all_DecOff = dq_all - dq_co2_ssp[...,np.newaxis]
        
        # Calculate the gradient od dq, dq_CO2 and dq_nonCO2 at the now index year
        dDelq_dt       = dq_all[now_index,:] - dq_all[now_index-1,:] 
        dDelqCO2_dt    = (dq_co2_ssp[now_index] - dq_co2_ssp[now_index-1])[...,np.newaxis]
        #dDelqNonCO2_dt = dq_non_co2_all_DecOff[now_index,:] - dq_non_co2_all_DecOff[now_index-1,:] 
        dDelqNonCO2_dt_ssp = (dq_non_co2_ssp[now_index] - dq_non_co2_ssp[now_index-1])[...,np.newaxis]

        # Calcualtate the residuals of the gradient and the absolute at present day stage:
        dDelQ_resid = dDelq_dt - dDelqCO2_dt - dDelqNonCO2_dt_ssp
        dQ_resid = dq_all[now_index,:] - dq_co2_ssp[now_index,np.newaxis] -  dq_non_co2_ssp[now_index,np.newaxis]

        epsilon_0_all_DecOff = dQ_resid[np.newaxis,...]

        epsilon_1_all_DecOff = (dDelQ_resid / dQ_resid)[np.newaxis,...]
        
        epsilon_t0    = epsilon_0_all_DecOff 
        #d_epsilon_t0  = epsilon_1_all_DecOff*epsilon_0_all_DecOff
        # sign of epsilon_t0 (offset)
        eps_t0_sign = epsilon_t0/np.abs(epsilon_t0)
        
        # calculate gamma as the minimum value that satisfies the following conditions:
        # 1. stability end_year, i.e. epsilon_t < threshold @ t=end_year
        # 2. epsilon and d2epsilon/dt2 must have opposite signs, such that the system behaves as an oscilator,
        #     i.e. returning to zero.
        #
        # 1:
        eps_stab = epsilon_stable * eps_t0_sign  # epsilon at stability with the correct sign
        t_stab   = end_year-yr_now               # time to stability
        gamma_all_DecOff = (  epsilon_1_all_DecOff*t_stab 
                           - np.log(eps_stab/epsilon_0_all_DecOff))/(t_stab**2)
        
        # 2: Condition statisfied when gamma > epsilon_1**2 / 2  (see paper/notes)
        index_2 = gamma_all_DecOff < (epsilon_1_all_DecOff**2.)/2
        gamma_all_DecOff[index_2] = (epsilon_1_all_DecOff[index_2]**2.)/2 
        
        # 3: Condition that gamma > -epsilon_1:
        index_3 = np.abs(gamma_all_DecOff)<np.abs(epsilon_1_all_DecOff)
        gamma_all_DecOff[index_3] = epsilon_1_all_DecOff[index_3]*-1.
        ipdb.set_trace()

        # Calculate the epsilon term as a function of time
        t = (yr_plot-yr_now)[...,np.newaxis]
        epsilon_t_all_DecOff = epsilon_0_all_DecOff * np.exp( (epsilon_1_all_DecOff*t) 
                                                            - (gamma_all_DecOff * t**2 ) )
        #epsilon_t_all_DecOff[:now_index,:] = 0.0
        
        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_DecOff[now_index:,:] =  dq_non_co2_ssp[now_index:,np.newaxis] + epsilon_t_all_DecOff[now_index:,:] 
        
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_DecOff = dq_all - dq_non_co2_all_DecOff

        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_DecOff = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_DecOff[:,igcm],
                                   n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                           for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_DecOff = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_DecOff)/q2co2)
        np.save(out_dir+'dq_co2_DecOff.npy',dq_co2_all_DecOff)
        np.save(out_dir+'dq_non_co2_DecOff.npy',dq_non_co2_all_DecOff)
        np.save(out_dir+'epsilon_0_DecOff.npy',epsilon_0_all_DecOff)
        np.save(out_dir+'epsilon_1_DecOff.npy',epsilon_1_all_DecOff)
        np.save(out_dir+'co2_ppm_DecOff.npy',co2_ppm_all_DecOff)

    if os.path.isfile(out_dir+'d_atmos_GtC_DecOff.npy'):
        d_land_GtC_all_DecOff   = np.load(out_dir+'d_land_GtC_DecOff.npy')
        d_ocean_GtC_all_DecOff  = np.load(out_dir+'d_ocean_GtC_DecOff.npy')
        co2_emissions_GtC_all_DecOff  = np.load(out_dir+'co2_emissions_GtC_DecOff.npy')
        d_atmos_GtC_all_DecOff  = np.load(out_dir+'d_atmos_GtC_DecOff.npy')
    else:
        print('Calculating Emissions for decay offset method')
        d_atmos_GtC_all_DecOff = co2ppm_to_datmosGtC(co2_ppm_all_DecOff)[...,np.newaxis]
        ipdb.set_trace()
    
        d_ocean_GtC_all_DecOff = co2nT_to_deltaOceanGtC(co2_ppm_all_DecOff, d_atmos_GtC_all_DecOff,
                                                     delta_temp_global_all,delta_temp_ocean_all,
                                                     ocean_model=ocean_model)

        d_land_GtC_all_DecOff = co2nT_to_deltaLandGtC(co2_ppm_all_DecOff, d_atmos_GtC_all_DecOff, d_ocean_GtC_all_DecOff, 
                                                    delta_temp_global_all, land_model=land_model )

        co2_emissions_GtC_all_DecOff = d_land_GtC_all_DecOff + d_atmos_GtC_all_DecOff + d_ocean_GtC_all_DecOff
        ipdb.set_trace()
        
        np.save(out_dir+'d_ocean_GtC_DecOff.npy',d_ocean_GtC_all_DecOff)
        np.save(out_dir+'d_land_GtC_DecOff.npy',d_land_GtC_all_DecOff)
        np.save(out_dir+'co2_emissions_GtC_DecOff.npy',co2_emissions_GtC_all_DecOff)
        np.save(out_dir+'d_atmos_GtC_DecOff.npy',d_atmos_GtC_all_DecOff)

    METHOD_DICT['DecayOffset']= {'dtemp_global':delta_temp_global_all,
                                 'dq':dq_all,
                                 'dq_co2':dq_co2_all_DecOff,
                                 'dq_non_co2':dq_non_co2_all_DecOff,
                                 'co2_ppm':co2_ppm_all_DecOff,
                                 'd_co2_GtC':d_atmos_GtC_all_DecOff,
                                 'co2_emissions':co2_emissions_GtC_all_DecOff,
                                 'ocean_uptake':d_ocean_GtC_all_DecOff,
                                 'land_uptake':d_land_GtC_all_DecOff,
                                 'cum_colour':'#fdbf6f',    # Shades of red for decay offset
                                 'colour':'#ff7f00',
                                 }


#ipdb.set_trace()
# T profile adjust method. Gradient is corrected with the Tprofile adjust + an offset to account for the absoulte difference
#  2 offset possibilities, constant offset = epsilon_all_TprfAd or time-decaying (efolding) offset = epsilon_all_TprfAdDOS
if l_Tprof_Adjust:
    if os.path.isfile(out_dir+'co2_ppm_TprfAd.npy'):
        delta_temp_global_all_TprfAd = np.load(out_dir+'delta_temp_global_TprfAd.npy')
        delta_temp_ocean_all_TprfAd  = np.load(out_dir+'delta_temp_ocean_TprfAd.npy')
        delta_temp_land_all_TprfAd   = np.load(out_dir+'delta_temp_land_TprfAd.npy')
        dq_all_TprfAd                = np.load(out_dir+'dq_TprfAd.npy')
        dq_co2_all_TprfAd            = np.load(out_dir+'dq_co2_TprfAd.npy')
        dq_non_co2_all_TprfAd        = np.load(out_dir+'dq_non_co2_TprfAd.npy')
        epsilon_all_TprfAd           = np.load(out_dir+'epsilon_TprfAd.npy')
        co2_ppm_all_TprfAd           = np.load(out_dir+'co2_ppm_TprfAd.npy')

    else:
        # First step calcualte the new trpfoiles using the gradient of the inverted RF
        # Assume that the ratio of the delQ and delT gradients are constant for the range we are adjusting for each GCM,
        #  denote with zeta:
        dDelq_dt_inv = dq_all[now_index,:] - dq_all[now_index-1,:]
        dDelT_dt_inv = delta_temp_global[now_index] - delta_temp_global[now_index-1]
        zeta_all = dDelT_dt_inv / dDelq_dt_inv
        # Now the beta of our new t-profiles is zeta multiplied by the prescribed Qco2 and Qnonco2 pathways:
        dDelq_dt_for = (dq_co2_ssp[now_index]-dq_co2_ssp[now_index-1]) + (dq_non_co2_ssp[now_index]-dq_non_co2_ssp[now_index-1])
        beta_all = zeta_all * (dDelq_dt_for)
        
        # Calculate the GCM specific tprofiles with gradients that match the Q profiles:
        delta_temp_global_all_TprfAd = np.array([profile.profile(beta_all[igcm], dt_now,dt_limit, 
                                                                 mu_zero, mu_one, 
                                                                 yr_now=yr_now, yr_beta=yr_beta)[:n_yr]
                                                  for igcm in range(n_cmip5) ]).transpose()
        np.save(out_dir+'delta_temp_global_TprfAd.npy', delta_temp_global_all_TprfAd)
        
        # Save the GCM specific land and ocean tprofiles
        delta_temp_ocean_all_TprfAd = delta_temp_global_all_TprfAd / (f_all + (1.0-f_all)*nu_all)
        np.save(out_dir+'delta_temp_ocean_TprfAd.npy', delta_temp_ocean_all_TprfAd)
        delta_temp_land_all_TprfAd = delta_temp_ocean_all_TprfAd * nu_all
        np.save(out_dir+'delta_temp_land_TprfAd.npy', delta_temp_land_all_TprfAd)

        # Now recalculate dq_all based on new t-profiles:
        dq_all_TprfAd = np.zeros_like(dq_all)
        # Have to transpose the arrays in and out of thie function:
        temp_gradient_out_all_TprfAd = parabolic.parabolic_parallel(kappa_all,delta_temp_ocean_all_TprfAd.transpose()).transpose()
        dq_factor = ( (1.0-f_all)*lambda_l_all*nu_all/f_all ) + lambda_o_all
        dq_all_TprfAd = (temp_gradient_out_all_TprfAd + (delta_temp_ocean_all_TprfAd*dq_factor)) * f_all
        np.save(out_dir+'dq_TprfAd.npy',dq_all_TprfAd)
        
        # Calculate the absolute offset at t=0
        # Historical dq_non_co2 is the residual between dq and dq_co2_ssp 
        dq_non_co2_all_TprfAd = dq_all_TprfAd - dq_co2_ssp[...,np.newaxis]
        
        # CONSTANT OFFSET:
        # Calculate the offset at the now index year
        epsilon_all_TprfAd = dq_non_co2_all_TprfAd[now_index,:]-dq_non_co2_ssp[now_index]
        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_TprfAd[now_index:,:] = dq_non_co2_ssp[now_index:,np.newaxis] + epsilon_all_TprfAd
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_TprfAd = dq_all_TprfAd - dq_non_co2_all_TprfAd
        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_TprfAd = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_TprfAd[:,igcm],
                                   n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                           for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_TprfAd = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_TprfAd)/q2co2)
        np.save(out_dir+'dq_co2_TprfAd.npy',dq_co2_all_TprfAd)
        np.save(out_dir+'dq_non_co2_TprfAd.npy',dq_non_co2_all_TprfAd)
        np.save(out_dir+'epsilon_TprfAd.npy',epsilon_all_TprfAd)
        np.save(out_dir+'co2_ppm_TprfAd.npy',co2_ppm_all_TprfAd)
    
    # Calculate emissions for the constant offset
    if os.path.isfile(out_dir+'d_atmos_GtC_TprfAd.npy'):
        d_land_GtC_all_TprfAd   = np.load(out_dir+'d_land_GtC_TprfAd.npy')
        d_ocean_GtC_all_TprfAd  = np.load(out_dir+'d_ocean_GtC_TprfAd.npy')
        co2_emissions_GtC_all_TprfAd  = np.load(out_dir+'co2_emissions_GtC_TprfAd.npy')
        d_atmos_GtC_all_TprfAd  = np.load(out_dir+'d_atmos_GtC_TprfAd.npy')
    else:
        print('Calculating Emissions for t-profile adjust with constant offset method')
        d_atmos_GtC_all_TprfAd = co2ppm_to_datmosGtC(co2_ppm_all_TprfAd)[...,np.newaxis]
        d_ocean_GtC_all_TprfAd = co2nT_to_deltaOceanGtC(co2_ppm_all_TprfAd, d_atmos_GtC_all_TprfAd,
                                                         delta_temp_global_all_TprfAd,delta_temp_ocean_all_TprfAd,
                                                         ocean_model=ocean_model)
        d_land_GtC_all_TprfAd = co2nT_to_deltaLandGtC(co2_ppm_all_TprfAd, d_atmos_GtC_all_TprfAd, d_ocean_GtC_all_TprfAd, 
                                                    delta_temp_global_all_TprfAd, land_model=land_model )
        co2_emissions_GtC_all_TprfAd = d_land_GtC_all_TprfAd + d_atmos_GtC_all_TprfAd + d_ocean_GtC_all_TprfAd

        np.save(out_dir+'d_ocean_GtC_TprfAd.npy',d_ocean_GtC_all_TprfAd)
        np.save(out_dir+'d_land_GtC_TprfAd.npy',d_land_GtC_all_TprfAd)
        np.save(out_dir+'co2_emissions_GtC_TprfAd.npy',co2_emissions_GtC_all_TprfAd)
        np.save(out_dir+'d_atmos_GtC_TprfAd.npy',d_atmos_GtC_all_TprfAd)

    METHOD_DICT['TprofAdjust']= {'dtemp_global':delta_temp_global_all_TprfAd,
                                 'dq':dq_all_TprfAd,
                                 'dq_co2':dq_co2_all_TprfAd,
                                 'dq_non_co2':dq_non_co2_all_TprfAd,
                                 'co2_ppm':co2_ppm_all_TprfAd,
                                 'd_co2_GtC':d_atmos_GtC_all_TprfAd,
                                 'co2_emissions':co2_emissions_GtC_all_TprfAd,
                                 'ocean_uptake':d_ocean_GtC_all_TprfAd,
                                 'land_uptake':d_land_GtC_all_TprfAd,
                                 'cum_colour':'#cab2d6',            # Shades of purple for Tprofil adjust with constant offset
                                 'colour':'#6a3d9a',
                                 }

# T profile adjust method. Gradient is corrected with the Tprofile adjust + an offset to account for the absoulte difference
#  offset possibility 2: time-decaying (efolding) offset = epsilon_t_all_TprfAdDOS
if l_Tprof_Adjust_DOS:
    if os.path.isfile(out_dir+'co2_ppm_TprfAdDOS.npy'):
        delta_temp_global_all_TprfAdDOS = np.load(out_dir+'delta_temp_global_TprfAdDOS.npy')
        delta_temp_ocean_all_TprfAdDOS  = np.load(out_dir+'delta_temp_ocean_TprfAdDOS.npy')
        delta_temp_land_all_TprfAdDOS   = np.load(out_dir+'delta_temp_land_TprfAdDOS.npy')
        dq_all_TprfAdDOS                = np.load(out_dir+'dq_TprfAdDOS.npy')
        dq_co2_all_TprfAdDOS            = np.load(out_dir+'dq_co2_TprfAdDOS.npy')
        dq_non_co2_all_TprfAdDOS        = np.load(out_dir+'dq_non_co2_TprfAdDOS.npy')
        epsilon_t_all_TprfAdDOS         = np.load(out_dir+'epsilon_t_TprfAdDOS.npy')
        co2_ppm_all_TprfAdDOS           = np.load(out_dir+'co2_ppm_TprfAdDOS.npy')
    else:         
        # Copy temp global and dq from the constant offset method.
        delta_temp_global_all_TprfAdDOS = delta_temp_global_all_TprfAd
        delta_temp_ocean_all_TprfAdDOS  = delta_temp_ocean_all_TprfAd
        delta_temp_land_all_TprfAdDOS   = delta_temp_land_all_TprfAd
        dq_all_TprfAdDOS = dq_all_TprfAd 

        # DECAYING OFFSET:
        dq_all_TprfAdDOS = dq_all_TprfAd
        # Historical dq_non_co2 is the residual between dq and dq_co2_ssp 
        dq_non_co2_all_TprfAdDOS = dq_all_TprfAdDOS - dq_co2_ssp[...,np.newaxis]
        # Calculate time decyaing espilon = A*e^(-gamma*t)
        epsilon_t_all_TprfAdDOS = epsilon_all_TprfAd[np.newaxis,...] * np.exp(-gamma*(yr_plot-yr_now))[...,np.newaxis]
        # for period after now_index, dq_non_co2 is the ssp value - offset
        dq_non_co2_all_TprfAdDOS[now_index:,:] = dq_non_co2_ssp[now_index:,np.newaxis] + epsilon_t_all_TprfAdDOS[now_index:,:]
        # Can now Calculate CO2 from dq and dq_non_co2, we do whole time period to check conservation
        dq_co2_all_TprfAdDOS = dq_all_TprfAdDOS - dq_non_co2_all_TprfAdDOS
        # Calculate the co2 due to the RF time series:
        if Etminan:
            # Iterative solution, so must be done gcm by gcm:
            co2_ppm_all_TprfAdDOS = np.array( [ delQ.etminan_CO2_inverse_series(dq_co2_all_TprfAdDOS[:,igcm],
                                               n2o_ppb_ssp, co2_ppm_0=co2_ppm_ssp[0], n2o_ppb_0=n2o_ppb_ssp[0])
                                               for igcm in range(n_cmip5) ]).transpose()
        else:
            co2_ppm_all_TprfAdDOS = co2_ppm_ssp[0] * np.exp((np.log(2.0)*dq_co2_all_TprfAdDOS)/q2co2)
        
        np.save(out_dir+'delta_temp_global_TprfAdDOS.npy', delta_temp_global_all_TprfAdDOS)
        np.save(out_dir+'delta_temp_ocean_TprfAdDOS.npy', delta_temp_ocean_all_TprfAdDOS)
        np.save(out_dir+'delta_temp_land_TprfAdDOS.npy', delta_temp_land_all_TprfAdDOS)
        np.save(out_dir+'dq_TprfAdDOS.npy',dq_all_TprfAdDOS)
        np.save(out_dir+'dq_co2_TprfAdDOS.npy',dq_co2_all_TprfAdDOS)
        np.save(out_dir+'dq_non_co2_TprfAdDOS.npy',dq_non_co2_all_TprfAdDOS)
        np.save(out_dir+'epsilon_t_TprfAdDOS.npy',epsilon_t_all_TprfAdDOS)
        np.save(out_dir+'co2_ppm_TprfAdDOS.npy',co2_ppm_all_TprfAdDOS)
    
    # Calculate emissions for the decaying offset
    if os.path.isfile(out_dir+'d_atmos_GtC_TprfAdDOS.npy'):
        d_land_GtC_all_TprfAdDOS   = np.load(out_dir+'d_land_GtC_TprfAdDOS.npy')
        d_ocean_GtC_all_TprfAdDOS  = np.load(out_dir+'d_ocean_GtC_TprfAdDOS.npy')
        co2_emissions_GtC_all_TprfAdDOS  = np.load(out_dir+'co2_emissions_GtC_TprfAdDOS.npy')
        d_atmos_GtC_all_TprfAdDOS  = np.load(out_dir+'d_atmos_GtC_TprfAdDOS.npy')
    else:
        print('Calculating Emissions for t-profile adjust with decaying offset method')
        d_atmos_GtC_all_TprfAdDOS = co2ppm_to_datmosGtC(co2_ppm_all_TprfAdDOS)[...,np.newaxis]
        d_ocean_GtC_all_TprfAdDOS = co2nT_to_deltaOceanGtC(co2_ppm_all_TprfAdDOS, d_atmos_GtC_all_TprfAdDOS,
                                                         delta_temp_global_all_TprfAdDOS, delta_temp_ocean_all_TprfAdDOS,
                                                         ocean_model=ocean_model)
        d_land_GtC_all_TprfAdDOS = co2nT_to_deltaLandGtC(co2_ppm_all_TprfAdDOS, d_atmos_GtC_all_TprfAdDOS, d_ocean_GtC_all_TprfAdDOS, 
                                                         delta_temp_global_all_TprfAdDOS, land_model=land_model )
        co2_emissions_GtC_all_TprfAdDOS = d_land_GtC_all_TprfAdDOS + d_atmos_GtC_all_TprfAdDOS + d_ocean_GtC_all_TprfAdDOS

        np.save(out_dir+'d_ocean_GtC_TprfAdDOS.npy', d_ocean_GtC_all_TprfAdDOS)
        np.save(out_dir+'d_land_GtC_TprfAdDOS.npy', d_land_GtC_all_TprfAdDOS)
        np.save(out_dir+'co2_emissions_GtC_TprfAdDOS.npy', co2_emissions_GtC_all_TprfAdDOS)
        np.save(out_dir+'d_atmos_GtC_TprfAdDOS.npy', d_atmos_GtC_all_TprfAdDOS)
    
    METHOD_DICT['TprofAdjustDecOff']= {'dtemp_global':delta_temp_global_all_TprfAdDOS,
                                       'dq':dq_all_TprfAdDOS,
                                       'dq_co2':dq_co2_all_TprfAdDOS,
                                       'dq_non_co2':dq_non_co2_all_TprfAdDOS,
                                       'co2_ppm':co2_ppm_all_TprfAdDOS,
                                       'd_co2_GtC':d_atmos_GtC_all_TprfAdDOS,
                                       'co2_emissions':co2_emissions_GtC_all_TprfAdDOS,
                                       'ocean_uptake':d_ocean_GtC_all_TprfAdDOS,
                                       'land_uptake':d_land_GtC_all_TprfAdDOS,
                                       'cum_colour':'#c88a68',            # Shades of brown for Tprofil adjust with decay offset
                                       'colour':'#b15928',
                                       }


PLOT_dir = out_dir+'plots/'
os.system('mkdir -p '+PLOT_dir)
plot_index=range(n_cmip5)
#plot_index=np.where(np.max(dq_all,axis=0)>2.5)[0]

xlim=[1950,2100]
FONTSIZE=35
SSP_colour = 'r'
Dtemp_OP_colour = 'k'

for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,28],ncols=1,nrows=5)
    colour=METHOD_DICT[method]['colour']
    # Plot Temperature Profile 
    ax=AXES[0]
    ax.plot(yr_plot,METHOD_DICT[method]['dtemp_global'],c=colour,lw=3)
    ax.set_title('Global Temperature Anomaly',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    #ax.grid(True)
    ax.set_xlim(xlim)

    # Plot total RF in top plot
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq'][:,i_gcm],c=colour) #,alpha=0.2)
    ax.plot(yr_plot,dq_non_co2_ssp+dq_co2_ssp, c=SSP_colour, lw=3, ls=':')
    #ax.plot(yr_plot,METHOD_DICT[method]['dtemp_global'],c='k',lw=3,ls=':')
    ax.set_title('Total Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    #ax.grid(True)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[2]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_non_co2'][:,i_gcm],c=colour) #,alpha=0.2)
    ax.plot(yr_plot,dq_non_co2_ssp, c=SSP_colour, lw=3, ls=':')
    #ax.grid(True)
    ax.set_title('Non-CO2 Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Plot CO2 concentration on fourth row
    ax=AXES[3]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_ppm'][:,i_gcm],c=colour) #,alpha=0.2)
    ax.plot(yr_plot,co2_ppm_ssp, c=SSP_colour, lw=3, ls=':')
    #ax.grid(True)
    ax.set_title('Atmospheric CO2',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Concentration (ppmv)',fontsize=FONTSIZE/2.)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Plot CO2 concentration on fifth row
    ax=AXES[4]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_emissions'][:,i_gcm],c=colour,alpha=0.3)
    #ax.grid(True)
    ax.set_title('Anthropogenic CO2 Emissions',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Emission (GtC yr$^{-1}$)',fontsize=FONTSIZE/2.)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    FIG.savefig(PLOT_dir+method+'_RFbreakdown.png',bbox_inches='tight')
    FIG.savefig(PLOT_dir+method+'_RFbreakdown.eps',bbox_inches='tight')
    plt.close()

for method in METHOD_DICT.keys():
    print(method)
    FIG,AXES=plt.subplots(figsize=[14,28],ncols=1,nrows=4)
    colour=METHOD_DICT[method]['colour']
    cum_colour=METHOD_DICT[method]['cum_colour']
    # Plot Temperature Profile 
    ax=AXES[0]
    ax2=ax.twinx()
    ax2.plot(yr_plot, METHOD_DICT[method]['dtemp_global'], c=Dtemp_OP_colour, lw=3)
    for i_gcm in plot_index:  
        plot_data     = METHOD_DICT[method]['co2_ppm'][:,i_gcm]
        ax.plot(yr_plot,plot_data,c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c=SSP_colour,lw=3)
    ax.set_title('Atmospheric CO2 Concetration',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Atmospheric CO2 (ppmv)',fontsize=FONTSIZE/2.)
    ax2.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    #ax.grid(True)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Plot CO2 concentration and cumulative uptake on fourth row
    ax=AXES[1]
    ax2=ax.twinx()
    for i_gcm in plot_index:  
        plot_data     = METHOD_DICT[method]['d_co2_GtC'][:,i_gcm]
        plot_cum_data = np.cumsum((plot_data)/GtC_to_ppm)
        ax.plot(yr_plot,plot_data,c=colour)
        ax2.plot(yr_plot,plot_cum_data,c=cum_colour)
    ax.set_title('Atmospheric Carbon Uptake',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Uptake (GtC)',fontsize=FONTSIZE/2.,color=colour)
    ax2.set_ylabel('Cumalative Uptake (GtC)',fontsize=FONTSIZE/2.,color=cum_colour)
    #ax.grid(True)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Plot Ocean Uptake and Cumulative Uptake:
    ax=AXES[2]
    ax2=ax.twinx()
    for i_gcm in plot_index:
        plot_data = METHOD_DICT[method]['ocean_uptake'][:,i_gcm]*-1./GtC_to_ppm
        plot_cum_data = np.cumsum(plot_data)
        ax.plot(yr_plot,plot_data,c=colour)
        ax2.plot(yr_plot,plot_cum_data,c=cum_colour)
    ax.set_title('Ocean Carbon Uptake',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Uptake (GtC)',fontsize=FONTSIZE/2.,color=colour)
    ax2.set_ylabel('Cumulative Uptake (GtC)',fontsize=FONTSIZE/2.,color=cum_colour)
    #ax.grid(True)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Plot Land Uptake
    ax=AXES[3]
    ax2=ax.twinx()
    for i_gcm in plot_index:
        plot_data = METHOD_DICT[method]['co2_emissions'][:,i_gcm]*0.25/GtC_to_ppm
        plot_cum_data = np.cumsum(plot_data,axis=0)
        dat_lin=ax.plot(yr_plot,plot_data,c=colour)
        cum_lin=ax2.plot(yr_plot,plot_cum_data,c=cum_colour)
    ax.grid(True)
    ax.set_title('Land Carbon Uptake',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Uptake (GtC)',fontsize=FONTSIZE/2.,color=colour)
    ax2.set_ylabel('Cumulative Uptake (GtC)',fontsize=FONTSIZE/2.,color=cum_colour)
    #ax.grid(True)
    ylim = ax.get_ylim()
    ax.plot([yr_now,yr_now],ylim,c='k',lw=1.3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    FIG.savefig(PLOT_dir+method+'_UptakeBreakdown.png',bbox_inches='tight')
    FIG.savefig(PLOT_dir+method+'_UptakeBreakdown.eps',bbox_inches='tight')
    plt.close()


quit()
    
    # LEGEND FOR UPTAKE PLOTS:
    #temp_han1,temp_lab1=ax.get_legend_handles_labels()
    #temp_han2,temp_lab2=ax2.get_legend_handles_labels()
    #print(temp_han1)
    #handles= [ dat_lin,cum_lin ]
    #labels = [ 'Annual','Cumlative' ]
    #FIG.legend( handles, labels, loc=8, fontsize=20, ncol=2 )
            



#print(np.max(dq_all[now_index,:]),np.min(dq_all[now_index,:]))
#print(np.max(dq_all[now_index,:])-np.min(dq_all[now_index,:]))
FONTSIZE=30.
for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,20],ncols=1,nrows=4)
    colour=METHOD_DICT[method]['colour']
    # Plot total RF in top plot
    ax=AXES[0]
    for i_gcm in plot_index:   #range(n_cmip5):
        #ax.plot(yr_plot,dq_all[:,i_gcm],c=colour)
        ax.plot(yr_plot,METHOD_DICT[method]['dq'],c='k',lw=3)
    ax.plot(yr_plot,dq_non_co2_ssp+dq_co2_ssp,c='r',lw=3)
    ax.set_ylabel('Total Radiative Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.grid(True)
    ax.text(yr_plot[2],((ax.get_ylim()[1]-ax.get_ylim()[0])*0.95)+ax.get_ylim()[0],
            '(a)', fontsize=FONTSIZE,verticalalignment='top')
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_non_co2'][:,i_gcm],c=colour)
    ax.plot(yr_plot,dq_non_co2_ssp, c='r',lw=3)
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
    FIG.savefig(PLOT_dir+method+'_RFbreakdown_forPaper.png',bbox_inches='tight')
    FIG.savefig(PLOT_dir+method+'_RFbreakdown_forPaper.eps',bbox_inches='tight')
    plt.close()


