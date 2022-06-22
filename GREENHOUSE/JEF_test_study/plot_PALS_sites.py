#!/bin/python
#
# Scripts to plot output from the JULES simulations alongside the field data.
# 
#  ./plot_PALS_sites.py [INTERACTIVE] [iSITE] [Tres] [PLOTS] [iDISPLAY]
#  ./plot_PALS_sites.py Y
#  ./plot_PALS_sites.py N 0 day YN N
# 
# 

import netCDF4 as nc
import netcdftime as nctime
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import glob, sys, os
import csv

###################################################
# Check for INTERACTIVE, add additional check for no arguments
if (len(sys.argv)>1):
    INTERACTIVE=sys.argv[1]
else:
    INTERACTIVE='Y'

if INTERACTIVE=='Y':
    iDISPLAY=raw_input('Display plots? [Y/N] ')
else:
    iDISPLAY=sys.argv[5]

########################################
# Constants and options
kgC_to_umolsCO2_factor = (1./12.)* 1e9
seconds_in_year = 60.*60.*24.*365.
SITE_colour = 'skyblue'
JULES_colour= 'green'


###################################
# Directories and filenames
PALS_dir = '/prj/GREENHOUSE/PALS_comparison/'
site_metadata_file = PALS_dir+'sites_meta.csv'
site_frac_file = PALS_dir+'sites_julesfracs.dat'
#
JULES_output_dir = PALS_dir+'output/'
#
SITE_data_dir = '/data/grp/fluxdata/PALS_sites/'
#
BASE_plot_dir = PALS_dir+'plots/'

###################################
# Read in meta data
meta_data=list(csv.reader(open(site_metadata_file,'r')))
meta_data_hdrs=meta_data.pop(0)

site_frac_dat = list(csv.reader(open(site_frac_file,'r'),delimiter='-'))
frac_list = [line[1][1:] for line in site_frac_dat]

###################################
# Select a site from the list of available sites
if (INTERACTIVE=='Y'):
    print 'Available sites: '
    print '%3s%15s%15s%10s%10s' % (' # ','Site','Country','latitude','longitude')
    for i, site_meta in zip(range(len(meta_data)),meta_data):
        print '%3s%15s%15s%10s%10s' % (i, site_meta[0].strip(), \
                                      site_meta[1].strip(), \
                                      site_meta[7].strip(), \
                                      site_meta[8].strip())
        iSITE = input('Enter flux site to plot: ')
else:
    iSITE = sys.argv[2]

for iSITE in SITES:
site_meta=meta_data[iSITE]
frac=np.array(frac_list[iSITE].split(),dtype='float')
site_name=site_meta[0].strip()
site_country=site_meta[1].strip()
site_vegtype=site_meta[2].strip()
site_canht=site_meta[3].strip()
site_measht=site_meta[4].strip()
site_lat=site_meta[7].strip()
site_lon=site_meta[8].strip()

print 'Site: '+site_name
######################################
# Select tstep or days
if (INTERACTIVE=='Y'):
    Tres_list = ['tstep','day']
    for i in range(len(Tres_list)):
        print i,': ',Tres_list[i]
    Tres = Tres_list[input('Select a temporal resolution: ')]
else:
    Tres = sys.argv[3]

######################################
# Construct JULES and site filenames:
JULES_fname=site_name+'.'+Tres+'.nc'
if Tres=='tstep':
    SITE_fname=site_name+'Fluxnet.1.4_flux.nc'
elif Tres=='day':
    SITE_fname=site_name+'Fluxnet.1.4_flux_daily.nc'

# open files:
Jinf = nc.Dataset(JULES_output_dir+JULES_fname,'r')
Sinf = nc.Dataset(SITE_data_dir+SITE_fname,'r')


##########################################################
# Append site name to BASE_plot_dir:
plot_dir=BASE_plot_dir#+site_name+'/'
if (INTERACTIVE=='Y'):
    print 'Plot output directory: '+plot_dir
    plot_dir_temp=raw_input('Enter alternative plot output directory or hit return to continue: ')
    if (len(plot_dir_temp)>0):
        plot_dir=plot_dir_temp

# Create plot_dir if doesn't already exist
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

############################################################
# Extract data and convert to same units (Site data units)
SITE_time_obj = nctime.num2date(Sinf.variables['time'][:],        \
                                units=Sinf.variables['time'].units)
SITE_NEE      = Sinf.variables['NEE'][:]
SITE_GPP      = Sinf.variables['GPP'][:]
SITE_LHF      = Sinf.variables['Qle'][:]
SITE_SHF      = Sinf.variables['Qh'][:]
SITE_Net_Annual_Cuptake = (np.mean(SITE_NEE)/kgC_to_umolsCO2_factor)*seconds_in_year*-1.


JULES_time_obj = nctime.num2date(Jinf.variables['time'][:],        \
                                units=Jinf.variables['time'].units)
JULES_NPP_GB = np.sum(Jinf.variables['npp_nuptake_out'][:].squeeze() \
                        * frac[:5],axis=1 )   \
                * kgC_to_umolsCO2_factor      
JULES_GPP_GB = np.sum( Jinf.variables['gpp'][:].squeeze() \
                        * frac[:5],axis=1 )   \
                * kgC_to_umolsCO2_factor     
JULES_resp_s = Jinf.variables['co2_soil_gb'][:].squeeze() \
                * kgC_to_umolsCO2_factor      
JULES_NEE    =  (JULES_NPP_GB - JULES_resp_s)*-1
JULES_LHF    = Jinf.variables['latent_heat'][:].squeeze() 
JULES_SHF    = np.sum( Jinf.variables['ftl'][:].squeeze() \
                           * frac, axis=1)
JULES_Net_Annual_Cuptake = (np.mean(JULES_NEE)/kgC_to_umolsCO2_factor)*seconds_in_year*-1.
JULES_Gross_Annual_Cuptake = (np.mean(JULES_GPP_GB)/kgC_to_umolsCO2_factor)*seconds_in_year

##########################################################
# Plot time-series of SHF, LHF, and NEE
PLOTS=''
if (INTERACTIVE=='Y'):
    PLOTS+=raw_input('Plot time-series? (Y/N) ')
else:
    PLOTS=sys.argv[4]

if (PLOTS[0]=='Y'):
    FIG = plt.figure(figsize=(12,12))
    
    AX = FIG.add_subplot(3,1,1)
    AX.plot(SITE_time_obj,SITE_SHF,lw=2,c=SITE_colour)
    AX.plot(JULES_time_obj,JULES_SHF,lw=1.5,c=JULES_colour)
    AX.set_ylabel('Sensible Heat Flux (W $m^{-2}$)')
    correlation = stats.pearsonr(SITE_SHF[:-2],JULES_SHF)
    AX.text(1.01,0.1,"Pearson's R = %6.2f"%(correlation[0]),\
            transform=AX.transAxes, fontsize=14)

    AX = FIG.add_subplot(3,1,2)
    AX.plot(SITE_time_obj,SITE_LHF,lw=2,c=SITE_colour)
    AX.plot(JULES_time_obj,JULES_LHF,lw=1.5,c=JULES_colour)
    AX.set_ylabel('Latent Heat Flux (W $m^{-2}$)')
    correlation = stats.pearsonr(SITE_LHF[:-2],JULES_LHF) 
    AX.text(1.01,0.1,"Pearson's R = %6.2f"%(correlation[0]),\
            transform=AX.transAxes, fontsize=14)

    AX = FIG.add_subplot(3,1,3)
    AX.plot(SITE_time_obj,SITE_NEE,label='Site',lw=2,c=SITE_colour)
    AX.plot(JULES_time_obj,JULES_NEE,label='JULES-ECOSSE-FUN',lw=1.5,c=JULES_colour)
    AX.plot(SITE_time_obj,SITE_GPP,lw=1.5,c=SITE_colour,ls=':')
    AX.plot(JULES_time_obj,JULES_GPP_GB*-1,lw=1.2,c=JULES_colour,ls=':')
    AX.set_ylabel('NEE and GPP ($\mu$mol $m^{-2} s^{-1}$)')
    correlation = stats.pearsonr(SITE_NEE[:-2],JULES_NEE)
    AX.text(1.01,0.1,"Pearson's R = %6.2f"%(correlation[0]),\
            transform=AX.transAxes, fontsize=14)
    AX.text(1.01,0.3,"JULES gross C Uptake = \n %6.3f kgC y$^{-1}$ m$^{-2}$"%\
            (JULES_Gross_Annual_Cuptake),\
            transform=AX.transAxes, fontsize=14)
    AX.text(1.01,0.5,"JULES net C Uptake = \n %6.3f kgC y$^{-1}$ m$^{-2}$"%\
            (JULES_Net_Annual_Cuptake),\
            transform=AX.transAxes, fontsize=14)
    AX.text(1.01,0.7,"SITE net C Uptake = \n %6.3f kgC y$^{-1}$ m$^{-2}$"%\
            (SITE_Net_Annual_Cuptake),\
            transform=AX.transAxes, fontsize=14)


    AX.legend( bbox_to_anchor=(0.5,-0.25),loc=8,borderaxespad=0.,ncol=2)
    
    FIG.tight_layout(rect=[0,0.05,0.8,0.96])
    FIG.suptitle(site_name+'flux tower time-series, obs and J-E-F',\
                    fontsize=24)
    FIG.text(0.78,0.92,'Country: '+site_country,fontsize=16)
    FIG.text(0.78,0.89,'PFT: '+site_vegtype,fontsize=16)
    FIG.text(0.78,0.86,'Lat: '+site_lat+'$^o$E',fontsize=16)
    FIG.text(0.78,0.83,'Lon: '+site_lon+'$^o$N',fontsize=16)
    FIG.text(0.78,0.80,'Can. hgt: '+site_canht+' m',fontsize=16)
    FIG.text(0.78,0.77,'Meas. hgt: '+site_measht+' m',fontsize=16)
    

    FIG.savefig(plot_dir+site_name+'_Flux_Comparison_Time-Series.png', \
                bbox_inches='tight')
    
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


##########################################################
# Plot scatter plot of LHF against NEE
if (INTERACTIVE=='Y'):
    PLOTS+=raw_input('Plot LHF/NEE scatter? (Y/N) ')

if (PLOTS[1]=='Y'):
    FIG = plt.figure(figsize=(12,12))

    AX = FIG.add_subplot(1,1,1)
    AX.scatter(SITE_LHF,SITE_GPP,c=SITE_colour,label='Site')
    AX.scatter(JULES_LHF,JULES_GPP_GB,c=JULES_colour,label='JULES-ECOSSE-FUN')
    
    AX.set_xlabel('Latent Heat Flux (W $m^{-2}$)')
    AX.set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)')
    AX.legend()

    FIG.tight_layout(rect=[0,0,1,0.94])
    FIG.suptitle(site_name+'('+site_country+') flux tower, LHF vs GPP',\
                    fontsize=24)
    FIG.savefig(plot_dir+site_name+'_Flux_Comparison_Time-Series.png', \
                bbox_inches='tight')
    

    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()
    


