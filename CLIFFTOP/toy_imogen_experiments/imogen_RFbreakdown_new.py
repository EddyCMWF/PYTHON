#!/bin/env python

import ipdb
import numpy as np
import matplotlib.pyplot as plt
import sys,os
import argparse

from imogen import parabolic, profile, ocean_co2, data_info, delQ

# Constants and convertors:
GtC_to_ppm=0.471  # Conversion factor for GtC to CO2 ppm

def myparser():
    parser = argparse.ArgumentParser(
                description="Plot the Radiative Forcing Breakdown associated with an IMOGEN configurarion" 
                )

    parser.add_argument("-T", "--T_pathway", choices=['1p5deg','1p81p5deg','2deg','3deg'], default='1p5deg',
                                help="Which temperature pathway do you want to take? default: %(default)s)")

    parser.add_argument("-B", "--Baseline_Scenario", choices=['SSP2-2.6_IMAGE'], default='SSP2-2.6_IMAGE',
                                help="Which Baseline Scenario do you want to use? default: %(default)s)")

    parser.add_argument("-E", "--etminan", choices=[True,False], default=True,
                                help="Use the Etminan RF equations? default: %(default)s)")

    parser.add_argument("-SD", "--scenario_dir", default='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/',
                                help="Directory containing the Scenario data? default: %(default)s)")

    parser.add_argument("-OD", "--out_dir", default='/prj/CLIFFTOP/IMOGEN_RFbeakdown/toy_jules/',
                                help="Out directory for plots and data? default: %(default)s)")
    
    parser.add_argument("-dq", "--l_recal_dq", choices=[True,False], default=True,
                                help="Recalculate the total RF even if file exists? default: %(default)s)")

    args = parser.parse_args()
    return(args.T_pathway, args.Baseline_Scenario, args.etminan, args.scenario_dir, args.out_dir, args.l_recal_dq)

def get_cmip5_variables():
    cmip5_runs = data_info.cmip5_runs()
    n_cmip5 = len(cmip5_runs)
    # Read in the EBM parameters
    kappa_all=data_info.kappa(cmip5_runs)
    lambda_l_all=data_info.lambda_l(cmip5_runs)
    lambda_o_all=data_info.lambda_o(cmip5_runs)
    nu_all=data_info.nu(cmip5_runs)
    f_all=data_info.ocean_frac(cmip5_runs)
    return(cmip5_runs,n_cmip5,kappa_all,lambda_l_all,lambda_o_all,nu_all,f_all)

def get_Tprof_params(SCENARIO):
    beta=0.025                # K/yr
    dt_now=0.89
    yr_now=2015
    yr_beta=1995
    end_year=2100
    if SCENARIO=='1p5deg':
        # First set the 3 parameters in the temperature curves
        dt_limit=1.5
        mu_zero = 0.08
        mu_one = 0.000
    elif SCENARIO=='1p81p5deg':
        # First set the 3 parameters in the temperature curves
        dt_limit=1.5
        mu_zero = -0.01
        mu_one = 0.00087
    elif SCENARIO=='2deg':
        # First set the 3 parameters in the temperature curves
        dt_limit=2.0
        mu_zero = 0.08
        mu_one = 0.00
    elif SCENARIO=='2deg':
        # First set the 3 parameters in the temperature curves
        dt_limit=3.0
        mu_zero = 0.065
        mu_one = 0.00
    return(beta,dt_now,yr_now,yr_beta,end_year,dt_limit,mu_zero,mu_one)

def output_parameters(out_dir,Baseline,SCENARIO,
                       mu_zero,mu_one,beta,dt_now,dt_limit,
                       start_year,yr_now,yr_beta,end_year):
    # Output paramters to file in output directory:
    outf=open(out_dir+'paramters.txt','w')
    outf.write('#'+'%20s,%10s\n'% ('Parameter','Value'))
    outf.write('%20s,%20s\n'%('Baseline Scenario',Baseline))
    outf.write('%20s,%20s\n'%('Temperature Pathway',SCENARIO))
    outf.write('%20s,%20.5f\n'%('$\mu_{0}$',mu_zero))
    outf.write('%20s,%20.5f\n'%('$\mu_{1}$',mu_one))
    outf.write('%20s,%20.5f\n'%('$\\beta}$',beta))
    outf.write('%20s,%20.5f\n'%('$T_{now}$',dt_now))
    outf.write('%20s,%20.5f\n'%('$T_{lim}$',dt_limit))
    outf.write('%20s,%20i\n'%('Start Year',start_year))
    outf.write('%20s,%20i\n'%('Now Year',yr_now))
    outf.write('%20s,%20i\n'%('Now Year',yr_beta))
    outf.write('%20s,%20i\n'%('End Year',end_year))
    outf.close()


def read_SCENARIO_data(SCENARIO_DIR,Baseline,n_yr):
    # Read in the non-CO2 RF from Baseline Template Scenario
    # Filenames:
    RAD_FORCING_FILE=SCENARIO_DIR+Baseline+'_qnonco2.txt'
    CO2_FILE=SCENARIO_DIR+Baseline+'_concs_co2.txt'
    CH4_N2O_FILE=SCENARIO_DIR+Baseline+'_concs_ch4_n2o.txt'
    
    # lines from raditative forcing file:
    RF_lines=open(RAD_FORCING_FILE, 'r').readlines()
    
    # lines from CO2 concentration file:
    CO2_lines=open(CO2_FILE,'r').readlines()
    
    # lines from CH4 and N2O concentration file:
    CH4_N2O_lines=open(CH4_N2O_FILE,'r').readlines()
    
    dq_non_co2_ssp = np.zeros(n_yr)
    co2_ppm_ssp=np.zeros(n_yr)
    ch4_ppb_ssp=np.zeros(n_yr)
    n2o_ppb_ssp=np.zeros(n_yr)
    ssp_year=np.zeros(n_yr)
    
    for iyr in range(n_yr):
        split=RF_lines[iyr].split()
        dq_non_co2_ssp[iyr]=float(split[1])
        split=CO2_lines[iyr].split()
        ssp_year[iyr]=int(split[0])
        co2_ppm_ssp[iyr]=float(split[1])
        split=CH4_N2O_lines[iyr].split()
        ch4_ppb_ssp[iyr]=split[1]
        n2o_ppb_ssp[iyr]=split[2]

    return dq_non_co2_ssp, co2_ppm_ssp, ch4_ppb_ssp, n2o_ppb_ssp, ssp_year


def calc_and_store_tprof( out_dir, SCENARIO, 
                         beta, dt_now, dt_limit, mu_zero, mu_one, 
                         n_yr,yr_now, yr_beta, yr_plot):
    delta_temp_global = profile.profile( beta, dt_now, dt_limit, mu_zero, mu_one,
                                         yr_now=yr_now, yr_beta=yr_beta)
    np.save(out_dir+'delta_temp_global.npy',delta_temp_global)
    delta_temp_global = delta_temp_global[:n_yr]
    out_Tprof_filename=out_dir+SCENARIO+'_global_temp_anomaly.dat'
    outf=open(out_Tprof_filename,'w')
    for year,Temp in zip(yr_plot,delta_temp_global):
        line='%4i  %8.4f\n'%(year,Temp)
        outf.write(line)
    outf.close()
    return delta_temp_global


def calc_and_store_dq_all(out_dir, delta_temp_global, n_yr,
                        n_cmip5, f_all, nu_all, kappa_all, lambda_l_all, lambda_o_all): 
    dq_all                  = np.zeros([n_yr,n_cmip5])
    delta_temp_ocean_yearly = np.array( 
                   [ delta_temp_global / (f_all[i_gcm] + (1.0-f_all[i_gcm])*nu_all[i_gcm]) 
                       for i_gcm in range(n_cmip5) ] )
    
    temp_gradient_out_all = parabolic.parabolic_parallel(kappa_all, delta_temp_ocean_yearly)
    
    factor=(((1.0-f_all)*lambda_l_all*nu_all)/f_all)+lambda_o_all

    dq_all = (temp_gradient_out_all.transpose() + (delta_temp_ocean_yearly.transpose()*factor))*f_all
    dq_all = dq_all.transpose()

    np.save(out_dir+'dq_all.npy',dq_all)
    return dq_all


def calc_and_store_dq_pathways_Offset(dq_all, dq_co2_ssp, dq_non_co2_ssp, 
                                        co2_ppm_ssp, n2o_ppb_ssp, now_index, n_cmip5, n_yr):
    
    # for the historical period, dq_non_co2 is the residual of dq-dq_co2
    dq_non_co2_all = dq_all - dq_co2_ssp

    # At now_yr, calcualte the offset between residual dq_non_co2
    # and the prescribed ssp dq_non_co2
    dq_non_co2_offset_all       = dq_non_co2_all[:,now_index]-dq_non_co2_ssp[now_index]
    # for period post now_yr, dq_non_co2 is the ssp value + offset
    dq_non_co2_all[:,now_index:]  = np.array( [ 
                        dq_non_co2_ssp[now_index:] + dq_non_co2_offset_all[igcm] 
                           for igcm in range(n_cmip5) ] )
     
    # now calculate the dq_CO2 by taking the difference between total and nonCO2
    #  The historical period will remain the same.
    dq_co2_all = dq_all - dq_non_co2_all
    co2_ppm_pi = co2_ppm_ssp[0]
    co2_ppm_all=np.zeros_like(dq_co2_all)
    for iyr in range(n_yr):
        for igcm in range(n_cmip5):
            co2_ppm_all[igcm,iyr] = delQ.etminan_CO2_inverse(dq_co2_all[igcm,iyr], 
                                                            co2_ppm_ssp[iyr],n2o_ppb_ssp[iyr], 
                                                            co2_ppm_0=co2_ppm_ssp[0],
                                                            n2o_ppb_0=n2o_ppb_ssp[0])
    np.save(out_dir+'dq_co2.npy',dq_co2_all)
    np.save(out_dir+'dq_non_co2.npy',dq_non_co2_all)
    np.save(out_dir+'co2_ppm.npy',co2_ppm_all)
    
    return dq_co2_all, dq_non_co2_all, co2_ppm_all


def main(SCENARIO, Baseline, Etminan, SCENARIO_DIR, out_dir, l_recal_dq):
    # Get cmip5 names and parameter values:
    cmip5_runs,n_cmip5,kappa_all,lambda_l_all,lambda_o_all,nu_all,f_all = get_cmip5_variables()
    # Get temperature pathway variables:
    beta,dt_now,yr_now,yr_beta,end_year,dt_limit,mu_zero,mu_one = get_Tprof_params(SCENARIO)
    # Start year depends on the Baseline scenario template 
    if 'SSP' in Baseline:
        start_year=1850
    elif 'RCP' in Baseline:
        start_year=1859
    n_yr=end_year-start_year+1                     # number of years
    yr_plot=np.arange(start_year,end_year+1)       # array of years
    now_index=np.where(yr_plot==yr_now)[0][0]      # index of yr_now

    # Create plot_tag and append to out_dir, create directory if necessary
    PLOT_TAG=Baseline+'_'+SCENARIO
    out_dir = out_dir+PLOT_TAG+'/'
    print(out_dir)
    os.system('mkdir -p '+out_dir)
    
    # Output the parameter set up used for reference:
    output_parameters(out_dir,Baseline,SCENARIO,
                       mu_zero,mu_one,beta,dt_now,dt_limit,
                       start_year,yr_now,yr_beta,end_year) 
    
    # Gather all the Baseline Scenario data (suffixed ssp)
    dq_non_co2_ssp, co2_ppm_ssp, ch4_ppb_ssp, n2o_ppb_ssp, ssp_year = \
            read_SCENARIO_data(SCENARIO_DIR,Baseline,n_yr)
    co2_ppm_pi=co2_ppm_ssp[0]
    ch4_ppb_pi=ch4_ppb_ssp[0]
    n2o_ppb_pi=n2o_ppb_ssp[0]
    if Etminan:
        dq_co2_ssp = delQ.etminan_CO2(co2_ppm_ssp,n2o_ppb_ssp,co2_ppm_0=co2_ppm_pi,n2o_ppb_0=n2o_ppb_pi)
    else:
        dq_co2_ssp= np.log(co2_ppm_ssp/co2_ppm_pi) * (q2co2/np.log(2.))
    
    # Create and output temperature pathway based on the parameters,
    # this may change to include variable T-pathways.
    delta_temp_global= calc_and_store_tprof( out_dir, SCENARIO, 
                                             beta, dt_now, dt_limit, mu_zero, mu_one, 
                                             n_yr,yr_now, yr_beta, yr_plot)
    #ipdb.set_trace() 
    if os.path.isfile(out_dir+'dq_all.npy')|l_recal_dq:
        dq_all = np.load(out_dir+'dq_all.npy')
        l_recal_emissions=False
    else:
        dq_all = calc_and_store_dq_all(out_dir, delta_temp_global, n_yr, 
                         n_cmip5, f_all, nu_all, kappa_all, lambda_l_all, lambda_o_all)
        l_recal_emissions=True
    
    dq_co2_all, dq_non_co2_all, co2_ppm_all = calc_and_store_dq_pathways_Offset(
                        dq_all, dq_co2_ssp, dq_non_co2_ssp, co2_ppm_ssp,n2o_ppb_ssp,
                        now_index, n_cmip5, n_yr )


    if os.path.isfile(out_dir+'delta_co2_GtC.npy')|l_recal_emissions:
        d_ocean_atmos_all=np.load(out_dir+'d_ocean_atmos.npy')
        co2_emissions_all=np.load(out_dir+'co2_emissions.npy')
        delta_co2_GtC_all=np.load(out_dir+'delta_co2_GtC.npy')
        #delta_co2_ppm_all= delta_co2_GtC_all*GtC_to_ppm
    else:
        d_atmos_co2_GtC_all, d_ocean_atmos_GtC_all, d_land_atmos_GtC_all, co2_emissions_all  \
            = calc_and_store_emissions_uptake(delta_temp_global, co2_ppm_all, f_all,nu_all, n_yr, n_cmip5 )
    

def calc_and_store_emissions_uptake( delta_temp_global, co2_ppm_all, f_all,nu_all, n_yr, n_cmip5):
    #d_ocean_atmos_all=np.zeros([n_yr,n_cmip5]) 
    #co2_emissions_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions')

    delta_temp_ocean_yearly_all = np.array([delta_temp_global / 
                                           (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    #delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)
    
    d_atmos_co2_ppm_all     = np.zeros_like(co2_ppm_all)
    d_atmos_co2_ppm_all[:,1:] = co2_ppm_all[:,1:]-co2_ppm_all[:,:-1]
    d_atmos_co2_GtC_all     = d_atmos_co2_ppm_all/GtC_to_ppm
    
    # Calculate ocean uptake using Joos (replace with Friedlingsteing?)
    d_ocean_atmos_all     = ocean_co2.ocean_co2_parallel(co2_ppm_all.transpose(),
                                         delta_temp_ocean_yearly_all.transpose()).transpose()
    d_ocean_atmos_GtC_all = d_ocean_atmos_all/GtC_to_ppm

    # Assuming 25% uptake by land means that atmsophere+ocean = 75% of emissions
    # For compatibililty with Friedlingstein, we define d_land expclicitly as a third of 
    #  atmosphere and ocean uptake (sign is from land to atmos):
    d_land_atmos_GtC_all  = -(d_atmos_co2_GtC_all-d_ocean_atmos_all)/3.0
    co2_emissions_all     = (d_atmos_co2_GtC_all-d_ocean_atmos_GtC_all)/0.75

    np.save(out_dir+'d_atmos_co2_GtC.npy',d_atmos_co2_GtC_all)
    np.save(out_dir+'d_ocean_atmos.npy',d_ocean_atmos_GtC_all)
    np.save(out_dir+'d_land_atmos.npy',d_land_atmos_GtC_all)
    np.save(out_dir+'co2_emissions.npy',co2_emissions_all)
    
    return d_atmos_co2_GtC_all, d_ocean_atmos_GtC_all, d_land_atmos_GtC_all, co2_emissions_all 

if __name__ == "__main__":
    SCENARIO, Baseline, Etminan, SCENARIO_DIR, out_dir, l_recal_dq   = myparser()
    #
    #main(SCENARIO, Baseline, Etminan, SCENARIO_DIR, out_dir, l_recal_dq)
    quit()
else:
    SCENARIO='1p5deg'
    Baseline='SSP2-2.6_IMAGE'
    Etminan=True
    SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/'
    out_dir='/prj/CLIFFTOP/IMOGEN_RFbeakdown/toy_jules/'


#break 

# PARSED:    SCENARIO ='1p5deg'  # '2deg'  #'1p5deg' #'1p81p5deg' 
# PARSED:    Etminan=True    # Always use Etminan RF equations
# PARSED:    Baseline='SSP2-2.6_IMAGE'
# PARSED:   SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/'
# PARSED:   out_dir='/prj/CLIFFTOP/IMOGEN_RFbeakdown/plots/toy_jules/'



METHOD_DICT={ 'Offset': {'dq_co2':dq_co2_all,
                         'dq_non_co2':dq_non_co2_all,
                         'co2_ppm':co2_ppm_all,
                         'd_co2_GtC':d_atmos_co2_GtC_all,
                         'co2_emissions':co2_emissions_all,
                         'ocean_uptake':d_ocean_atmos_GtC_all,
                         'land_uptake':d_land_atmos_GtC_all,
                         'cum_colour':'g',
                         'colour':'b'},               
              } 

plot_index=range(n_cmip5)
#plot_index=np.where(np.max(dq_all,axis=0)>2.5)[0]

FONTSIZE=35

for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,28],ncols=1,nrows=5)
    colour=METHOD_DICT[method]['colour']
    # Plot Temperature Profile 
    ax=AXES[0]
    ax.plot(ssp_year,delta_temp_global,c='g',lw=3)
    ax.set_title('Global Temperature Anomaly',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    # Plot total RF in top plot
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(ssp_year,dq_all[i_gcm,:],c=colour)
    ax.plot(ssp_year,dq_non_co2_ssp+dq_co2_ssp,c='r',lw=3)
    ax.plot(ssp_year,delta_temp_global,c='k',lw=3,ls=':')
    ax.set_title('Total Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[2]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_non_co2'][i_gcm,:],c=colour)
    ax.plot(yr_plot,dq_non_co2_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Non-CO2 Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fourth row
    ax=AXES[3]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_ppm'][i_gcm,:],c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Atmospheric CO2',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Concentration (ppmv)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fifth row
    ax=AXES[4]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_emissions'][i_gcm,:],c=colour)
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
    FIG.savefig(out_dir+method+'_RFbreakdown_forPaper.png',bbox_inches='tight')
    FIG.savefig(out_dir+method+'_RFbreakdown_forPaper.eps',bbox_inches='tight')
    plt.close()


for method in METHOD_DICT.keys():
    FIG,AXES=plt.subplots(figsize=[14,28],ncols=1,nrows=4)
    colour=METHOD_DICT[method]['colour']
    cum_colour=METHOD_DICT[method]['cum_colour']
    # Plot Temperature Profile 
    ax=AXES[0]
    ax2=ax.twinx()
    ax2.plot(ssp_year,delta_temp_global,c='k',lw=3)
    for i_gcm in plot_index:  
        plot_data     = METHOD_DICT[method]['co2_ppm'][:,i_gcm]
        ax.plot(yr_plot,plot_data,c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c='r',lw=3)
    ax.set_title('Atmospheric CO2 Concetration',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Atmospheric CO2 (ppmv)',fontsize=FONTSIZE/2.)
    ax2.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    # Plot CO2 concentration and cumulative uptake on fourth row
    ax=AXES[1]
    ax2=ax.twinx()
    for i_gcm in plot_index:  
        plot_data     = METHOD_DICT[method]['d_co2_GtC'][:,i_gcm]
        plot_cum_data = np.cumsum((plot_data)/GtC_to_ppm)
        ax.plot(yr_plot,plot_data,c=colour)
        ax2.plot(yr_plot,plot_cum_data,c=cum_colour)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Atmospheric Carbon Uptake',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Uptake (GtC)',fontsize=FONTSIZE/2.,color=colour)
    ax2.set_ylabel('Cumalative Uptake (GtC)',fontsize=FONTSIZE/2.,color=cum_colour)

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
    ax.grid(True)

    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot Land Uptake
    ax=AXES[3]
    ax2=ax.twinx()
    for i_gcm in plot_index:
        plot_data = METHOD_DICT[method]['co2_emissions'][:,i_gcm]*0.25/GtC_to_ppm
        plot_cum_data = np.cumsum(plot_data)
        dat_lin=ax.plot(yr_plot,plot_data,c=colour)
        cum_lin=ax2.plot(yr_plot,plot_cum_data,c=cum_colour)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Land Carbon Uptake',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Uptake (GtC)',fontsize=FONTSIZE/2.,color=colour)
    ax2.set_ylabel('Cumulative Uptake (GtC)',fontsize=FONTSIZE/2.,color=cum_colour)
    ax.grid(True)

    #pdb.set_trace() 
    #temp_han1,temp_lab1=ax.get_legend_handles_labels()
    #temp_han2,temp_lab2=ax2.get_legend_handles_labels()
    #print(temp_han1)
    #handles= [ dat_lin,cum_lin ]
    #labels = [ 'Annual','Cumlative' ]
    #FIG.legend( handles, labels, loc=8, fontsize=20, ncol=2 )
            
    FIG.savefig(out_dir+method+'_UptakeBreakdown.png',bbox_inches='tight')
    FIG.savefig(out_dir+method+'_UptakeBreakdown.eps',bbox_inches='tight')
    plt.close()

