#!/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import parabolic, profile, ocean_co2, data_info

l_outTprof=False   # Output temperature profile?
SCENARIO = '1p5deg'
Etminan=True
Baseline='SSP2-2.6_IMAGE'
version='vn1p2'   # For saving output driving files

SMOOTH_CO2 = False # To smooth the co2 profile prior to calculation of emissions

l_SSP_Tprof=True

cmip5_runs = data_info.cmip5_runs()
n_cmip5 = len(cmip5_runs)
#n_cmip5=3

# RF equation Variables
q2co2=3.74                # Old CO2 RF parameter
# Etminan:
a1=-2.4e-7
b1=7.2e-4
c1=-2.1e-4

# Conversion factor for GtC to CO2 ppm
GtC_to_ppm=0.471

SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/'

if 'SSP' in Baseline.upper():
    if Etminan:
        RAD_FORCING_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2.txt'
        PLOT_TAG='SSP2-2.6_neweqns_'+SCENARIO
    else:
        RAD_FORCING_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2_splice.txt'
        PLOT_TAG='SSP2-2.6_oldeqns_'+SCENARIO

    CO2_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2.txt'
    CH4_N2O_FILE=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o.txt'

elif 'RCP' in Baseline.upper():
    RAD_FORCING_FILE=SCENARIO_DIR+'rcp_qnonco2_vn1p0_2.6.txt'
    CO2_FILE=SCENARIO_DIR+'rcp2.6_concs_co2_vn1p0.txt'
    CH4_N2O_FILE=SCENARIO_DIR+'rcp2.6_concs_ch4_n2o_vn1p0.txt'
    PLOT_TAG='RCP2.6'

if SMOOTH_CO2:
    PLOT_TAG+='_smoothCO2'

if l_SSP_Tprof:
    PLOT_TAG+='_SSP_Tprof'

#RCP_RAD_FORCING_FILE=SCENARIO_DIR+'rcp_qnonco2_vn1p0_2.6.txt'

#PLOT_TAG+='_mod_T'

print(PLOT_TAG)
out_dir='/prj/CLIFFTOP/IMOGEN_RFbeakdown/plots/toy_jules/'+PLOT_TAG+'/'
os.system('mkdir -p '+out_dir)

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
   
# Also normalise the curves so that they end up at the current temperature and gradient estimate
beta=0.025                # K/yr
dt_now=0.89
yr_now=2015
end_year=2100
if 'SSP2-2.6' in PLOT_TAG:
    start_year=1850
elif PLOT_TAG=='RCP2.6':
    start_year=1859
n_yr=end_year-start_year+1

if l_SSP_Tprof:
    Tprof_file=SCENARIO_DIR+'SSP2-2.6_IMAGE_temp_anomaly_'+SCENARIO+'.dat'
    #Tprof_file=SCENARIO_DIR+SCENARIO+'_global_temp_anomaly.dat'
    delta_temp_global=[]
    yr_plot=[]
    for line in open(Tprof_file,'r').readlines():
        split=line.split()
        yr_plot.append(int(split[0]))
        delta_temp_global.append(float(split[1]))
    
    yr_plot=np.array(yr_plot)
    delta_temp_global=np.array(delta_temp_global)
else:
    delta_temp_global = profile.profile(beta, dt_now, dt_limit, mu_zero, mu_one)
    delta_temp_global = delta_temp_global[:n_yr]
    yr_plot=np.arange(start_year,end_year+1)

print(yr_plot)
now_index=np.where(yr_plot==yr_now)[0][0]

# Output Tprofile if required
if l_outTprof:
    out_Tprof_filename=SCENARIO_DIR+SCENARIO+'_global_temp_anomaly_'+version+'.dat'
    outf=open(out_Tprof_filename,'w')
    for year,Temp in zip(yr_plot,delta_temp_global):
        line='%4i  %8.4f\n'%(year,Temp)
        outf.write(line)
    outf.close()
    
    
kappa_all=data_info.kappa(cmip5_runs)
lambda_l_all=data_info.lambda_l(cmip5_runs)
lambda_o_all=data_info.lambda_o(cmip5_runs)
nu_all=data_info.nu(cmip5_runs)
f_all=data_info.ocean_frac(cmip5_runs)


# Loop over the different GCMs emulated and calculate RF via dtemp_o
temp_gradient_out_all=np.zeros([n_yr,n_cmip5]) 
#  HadGEM-ES=17, CSIRO-Q=8, NOAA-2G=29
if os.path.isfile(out_dir+'dq_all.npy'):
    dq_all=np.load(out_dir+'dq_all.npy')
else:
    dq_all=np.zeros([n_yr,n_cmip5]) 
    for i_gcm in range(n_cmip5):
        print(i_gcm,':',cmip5_runs[i_gcm])
        delta_temp_ocean_yearly = delta_temp_global / (f_all[i_gcm] + (1.0-f_all[i_gcm])*nu_all[i_gcm])
        temp_gradient_out = parabolic.parabolic(kappa_all[i_gcm], delta_temp_ocean_yearly)
        print(temp_gradient_out.shape[0])
        temp_gradient_out_all[:,i_gcm]=temp_gradient_out
        dq = np.zeros(n_yr)
        # Derive time-evolution of radiative forcing, Q
        factor=(((1.0-f_all[i_gcm])*lambda_l_all[i_gcm]*nu_all[i_gcm])/f_all[i_gcm])+lambda_o_all[i_gcm]
        dq = (temp_gradient_out + delta_temp_ocean_yearly*factor)*f_all[i_gcm]
        dq_all[:,i_gcm]=dq
    np.save(out_dir+'dq_all.npy',dq_all)
    #np.save(out_dir+'delta_temp_ocean.npy',delta_temp_ocean_yearly)

# Read in the non-CO2 RF from PIK site
dq_non_co2 = np.zeros(n_yr)
dq_non_co2_ssp = np.zeros(n_yr)

#ECP method from ssp2-2.6_image data
f_in = open(RAD_FORCING_FILE, 'r')
RF_lines=f_in.readlines()
f_in.close()

co2_f_in=open(CO2_FILE,'r')
CO2_lines=co2_f_in.readlines()
co2_f_in.close()

ch4_n2o_f_in=open(CH4_N2O_FILE,'r')
CH4_N2O_lines=ch4_n2o_f_in.readlines()
ch4_n2o_f_in.close()

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


co2_ppm_pi=co2_ppm_ssp[0]

if Etminan:
    print('here')
    Nbar=(n2o_ppb_ssp+n2o_ppb_ssp[0])/2.
    dq_co2_ssp=    np.log(co2_ppm_ssp/co2_ppm_pi)      \
                 * (  (a1*((co2_ppm_ssp-co2_ppm_pi)**2))      \
                    + (b1*np.abs(co2_ppm_ssp-co2_ppm_pi))     \
                    + (c1*Nbar) + 5.36 )
else:
    dq_co2_ssp= np.log(co2_ppm_ssp/co2_ppm_pi) * (q2co2/np.log(2.))

now_index=np.where(yr_plot==yr_now)[0][0]

#emiss_invert=np.zeros([n_yr,n_cmip5])     # Holds the final backed-out emissions (GtC/yr)
#  HadGEM-ES=17, CSIRO-Q=8, NOAA-2G=29
if False: #os.path.isfile(out_dir+'co2_ppm_nutin.npy'):
    dq_co2_all=np.load(out_dir+'dq_co2_nutin.npy')
    dq_non_co2_all=np.load(out_dir+'dq_non_co2_nutin.npy')
    co2_ppm_all=np.load(out_dir+'co2_ppm_nutin.npy')
    
else:
    #array for all temp gradients and rfs
    dq_non_co2_all=np.zeros([n_yr,n_cmip5])
    dq_co2_all=np.zeros([n_yr,n_cmip5]) 
    dq_non_co2_offset_all=np.zeros([n_cmip5])
    co2_ppm_all=np.zeros([n_yr,n_cmip5]) 

    for i_gcm in range(n_cmip5):
        # Now back out CO2 concentrations and non_co2 radiative forcings
        dq=np.copy(dq_all[:,i_gcm])        
        dq_non_co2=dq_non_co2_ssp 
        dq_co2 = dq-dq_non_co2

        if Etminan:
            co2_ppm_pi=co2_ppm_ssp[0]
            co2_ppm=np.copy(co2_ppm_ssp)
            for iyr in range(n_yr):
                Nbar=(n2o_ppb_ssp[iyr]+n2o_ppb_ssp[0])*0.5
                # iterate to find CO2 concentration from radiative forcing.
                co2_OF = np.copy(co2_ppm[iyr-1])
                itercnt=0
                co2_iter=999.
                while ( (np.abs(co2_iter-co2_OF)>0.001) & (itercnt<=1000) ):
                    co2_iter = co2_OF.copy()
                    denom = (  (a1*((co2_iter-co2_ppm_pi)**2.))     \
                            + (b1*(co2_iter-co2_ppm_pi))           \
                            + (c1*Nbar)+5.36)
                    
                    co2_OF = co2_ppm_pi * np.exp( dq_co2[iyr]/denom )
                    itercnt+=1
                    co2_ppm[iyr]=np.copy(co2_OF)
                    
            co2_ppm = co2_ppm_pi * np.exp((np.log(2.0)*dq_co2)/q2co2)
        
        co2_ppm_all[:,i_gcm]=np.copy(co2_ppm)
        dq_co2_all[:,i_gcm]=np.copy(dq_co2)
        dq_non_co2_all[:,i_gcm]=np.copy(dq_non_co2)
    
    np.save(out_dir+'dq_co2_nutin.npy',dq_co2_all)
    np.save(out_dir+'dq_non_co2_nutin.npy',dq_non_co2_all)
    np.save(out_dir+'co2_ppm_nutin.npy',co2_ppm_all)
        
if os.path.isfile(out_dir+'d_ocean_atmos_nutin.npy'):
    d_ocean_atmos_all=np.load(out_dir+'d_ocean_atmos_nutin.npy')
    co2_emissions_all=np.load(out_dir+'co2_emissions_nutin.npy')
else:
    d_ocean_atmos_all=np.zeros([n_yr,n_cmip5]) 
    co2_emissions_all=np.zeros([n_yr,n_cmip5]) 
    print('Calculating Emissions')

    delta_temp_ocean_yearly_all = np.array([delta_temp_global / 
                                   (f + (1.0-f)*nu) for f,nu in zip(f_all,nu_all)])
    delta_temp_ocean_yearly_all = delta_temp_ocean_yearly_all.transpose(1,0)

    d_ocean_atmos_all     = ocean_co2.ocean_co2_parallel(co2_ppm_all,
                                                         delta_temp_ocean_yearly_all)
    d_ocean_atmos_GtC_all = d_ocean_atmos_all/GtC_to_ppm
    delta_co2_ppm_all     = np.zeros_like(co2_ppm_all)
    delta_co2_ppm_all[1:] = co2_ppm_all[1:,:]-co2_ppm_all[:-1,:]
    delta_co2_GtC_all     = delta_co2_ppm_all/GtC_to_ppm
    co2_emissions_all     = (delta_co2_GtC_all-d_ocean_atmos_GtC_all)/0.75
        
    np.save(out_dir+'d_ocean_atmos_nutin.npy',d_ocean_atmos_all)
    np.save(out_dir+'co2_emissions_nutin.npy',co2_emissions_all)


fig,axes=plt.subplots(ncols=1,nrows=4,figsize=[10,15])
for i_gcm in range(n_cmip5):
    delta_temp_ocean_yearly = delta_temp_global / (f_all[i_gcm] + (1.0-f_all[i_gcm])*nu_all[i_gcm])
    axes[0].plot(yr_plot,delta_temp_ocean_yearly)
axes[0].set_title('Ocean $\Delta$Temperature')
axes[0].set_ylabel('$\Delta$Temperature (K)')
axes[0].set_ylim([0.5,1.0])

for i_gcm in range(n_cmip5):
    axes[1].plot(yr_plot,d_ocean_atmos_all[:,i_gcm]/GtC_to_ppm)
axes[1].set_title('Ocean CO$_2$ Uptake')
axes[1].set_ylabel('Carbon Uptake (GtC $^{-1}$)')
axes[1].set_ylim([-3,-1])

for i_gcm in range(n_cmip5):
    co2_ppm = co2_ppm_all[:,i_gcm]
    delta_co2_ppm = np.zeros_like(co2_ppm)
    delta_co2_ppm[1:] = co2_ppm[1:]-co2_ppm[:-1]
    delta_co2_GtC = delta_co2_ppm/GtC_to_ppm
    axes[2].plot(yr_plot,delta_co2_ppm)
axes[2].set_title('$\Delta$CO$_2$ - Atmosphere')
axes[2].set_ylabel('$\Delta$CO$_2$ (ppmv)')

for i_gcm in range(n_cmip5):
    co2_ppm = co2_ppm_all[:,i_gcm]
    axes[3].plot(yr_plot,co2_ppm)
axes[3].set_title('CO$_2$ - Atmosphere')
axes[3].set_ylabel('CO$_2$ (ppmv)')
axes[3].set_ylim([380,420])

SY,EY=2014,2020
xticklabels=[str(year) for year in range(SY,EY+1)]
print(xticklabels)
for ax in axes:
    ax.set_xlim([SY,EY])
    ax.set_xticklabels(xticklabels) 

fig.savefig(out_dir+'docean.png',bbox_inches='tight')
plt.close()




METHOD_DICT={ 'Nothing': {'dq_co2':dq_co2_all,        
                         'dq_non_co2':dq_non_co2_all,
                         'co2_ppm':co2_ppm_all,      
                         'co2_emissions':co2_emissions_all,      
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
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    ax.set_title('Global Temperature Anomaly',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('$\delta$ Temperature (K)',fontsize=FONTSIZE/2.)
    ax.grid(True)

    # Plot total RF in top plot
    ax=AXES[1]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(ssp_year,dq_all[:,i_gcm],c=colour)

    ax.plot(ssp_year,dq_non_co2_ssp+dq_co2_ssp,c='r',lw=3)
    #ax.plot(ssp_year,delta_temp_global,c='k',lw=3,ls=':')
    ax.set_title('Total Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=1.3)
    
    # Plot NON-CO2 radiative forcing on second row
    ax=AXES[2]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['dq_non_co2'][:,i_gcm],c=colour)
    ax.plot(yr_plot,dq_non_co2_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Non-CO2 Radiative Forcing',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Rad. Forcing (W m$^{-2}$)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fourth row
    ax=AXES[3]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_ppm'][:,i_gcm],c=colour)
    ax.plot(yr_plot,co2_ppm_ssp, c='r',lw=3)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Atmospheric CO2',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Concentration (ppmv)',fontsize=FONTSIZE/2.)
    
    # Plot CO2 concentration on fifth row
    ax=AXES[4]
    for i_gcm in plot_index:   #range(n_cmip5):
        ax.plot(yr_plot,METHOD_DICT[method]['co2_emissions'][:,i_gcm],c=colour)
    ax.grid(True)
    ax.plot([yr_now,yr_now],ax.get_ylim(),c='k')
    ax.set_title('Anthropogenic CO2 Emissions',fontsize=2+FONTSIZE/2.)
    ax.set_ylabel('Annual Emission (GtC yr$^{-1}$)',fontsize=FONTSIZE/2.)
    
    #FIG.suptitle(PLOT_TAG+', Imogen vs '+PLOT_TAG+' '+method,fontsize=FONTSIZE)
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


