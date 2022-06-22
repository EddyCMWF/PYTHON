#!/bin/python
#
# Scripts to plot output from the JULES simulations alongside the field data.
# 
#  ./plot_PALS_sites_hourly_pandas.py [INTERACTIVE] [iDISPLAY] [PLOTS] [SITES] [Jsources] [iTRES]
#  ./plot_PALS_sites_hourly_pandas.py Y
#  ./plot_PALS_sites_hourly_pandas.py N N YNNNNNNNNN [0,1,2,3] [0,1] 0
# 
# 
# PLOTS[0] - Time Series
# PLOTS[1] - LHF & SHF scatter
# PLOTS[2] - GPP & NEE scatter
# PLOTS[3] - GPPvLHF scatter
# PLOTS[4] - Blank
# PLOTS[5] - GPP/NEE Scatter and Taylor for all sites
# PLOTS[6] - LHF/SHF Scatter and Taylor for all sites 
# PLOTS[7] - Blank
# PLOTS[8] - Blank
# PLOTS[9] - Blank
# PLOTS[10] - Diurnal plots full sim
# PLOTS[11] - Diurnal plots by season ( DJF, MAM, JJA, SON )
#

import netCDF4 as nc
import netcdftime as nctime
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import glob, sys, os
import csv
import pandas as pd
import PlotTools.plot_tools as PT

# ECP modules
import data_info_PALS as datainfo
from maths_tools import TimeSeriesTools as TST

###################################################
# Check for INTERACTIVE, add additional check for no arguments
if (len(sys.argv)>1):
    INTERACTIVE=sys.argv[1]
else:
    INTERACTIVE='Y'

if INTERACTIVE=='Y':
    iDISPLAY=raw_input('Display plots? [Y/N] ')
else:
    iDISPLAY=sys.argv[2]

########################################
# Constants and options
kgC_to_umolsCO2_factor = (1./12.)* 1e9
seconds_in_year = 60.*60.*24.*365.
Lc_H2O = 2.501e6
SITE_colour = 'red'
JULES_colour= 'green'

#LC_names = ['Broadleaf','Needleleaf','C3 grass','C4 Grass','shrub',\
#            'urban','lake','soil','ice']
LC_names = ['Tree','Tree','Grass','Grass','shrub',\
            'urban','lake','soil','ice']
nLCs = len(LC_names)
LC_markers = [ 'x','+','v','^','D',\
               's','p','o','h']
LC_mews    = [ 2,2,0,0,0,\
               0,0,0,0]
LC_symsizes= [ 80 for i in range(nLCs) ]
LC_handles = [ plt.Line2D((0,1),(0,0), color='k', marker=LC_markers[iLC], linestyle='') for iLC in range(nLCs) ]

SEAS_names=['DJF','MAM','JJA','SON']
nSEAS=4

###################################
# Directories and filenames
PALS_dir = '/prj/GREENHOUSE/PALS_comparison/'
site_metadata_file = PALS_dir+'sites_meta_plot.csv'
site_frac_file = PALS_dir+'sites_julesfracs.dat'
#
JULES_output_dir = PALS_dir+'output/'
#
SITE_data_dir = '/data/grp/fluxdata/PALS_sites_ECP/'
#
BASE_plot_dir = PALS_dir+'plots/'

###################################
# Read in meta data
meta_data=list(csv.reader(open(site_metadata_file,'r')))
meta_data_hdrs=meta_data.pop(0)

site_frac_dat = list(csv.reader(open(site_frac_file,'r'),delimiter='-'))
frac_list = [line[1][1:] for line in site_frac_dat]

###################################
# Fetch JULES sources info
JULES_sources = datainfo.JULES_sources()

######################################################################################
#  User input/ Parsed information section
######################################################################################
##  Select Plots
PLOTS=''
if (INTERACTIVE=='Y'):
    PLOTS+=raw_input('Plot time-series? (Y/N) ')  # PLOTS[0]
    PLOTS+=raw_input('Plot LHF and SHF scatter? (Y/N) ')  # PLOTS[1] 
    PLOTS+=raw_input('Plot NEE and GPP scatter? (Y/N) ')  # PLOTS[2]
    PLOTS+=raw_input('Plot LHF vs GPP scatter? (Y/N) ')  # PLOTS[3]
    PLOTS+=raw_input('Plot CO2 flux time-series breakdown? (Y/N) ')  # PLOTS[4]
    PLOTS+=raw_input('Plot scatter and Taylor diagrams of C uptake for all sites? (Y/N) ')  # PLOTS[5]
    PLOTS+=raw_input('Plot scatter and Taylor diagrams of Heat  uptake for all sites? (Y/N) ')  # PLOTS[6]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[7]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[8]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[9]
    PLOTS+=raw_input('Plot diurnal behaviour of GPP and NEE (Y/N) ')  # PLOTS[10]
    PLOTS+=raw_input('Plot seasonal diurnal behaviour of GPP and NEE (Y/N) ')  # PLOTS[11]
    PLOTS+=raw_input('Plot scatter and Taylor diagrams of diurnal amplitudes for all site (Y/N) ')  # PLOTS[12]
    PLOTS+=raw_input('Plot seasonal scatter and Taylor diagrams of diurnal amplitudes for all site (Y/N) ')  # PLOTS[13]
    PLOTS+=raw_input('Plot diurnal behaviour of GPP, NEE and RESP (Y/N) ')  # PLOTS[14]
    PLOTS+=raw_input('Plot seasonal diurnal behaviour of GPP and NEE and RESP (Y/N) ')  # PLOTS[15]
else:
    PLOTS=sys.argv[3]

################################################################################
# Select a site from the list of available sites
if (INTERACTIVE=='Y'):
    print 'Available sites: '
    print '%3s%15s%15s%10s%10s' % (' # ','Site','Country','latitude','longitude')
    for i, site_meta in zip(range(len(meta_data)),meta_data):
        print '%3s:%15s%15s%10s%10s' % (i, site_meta[0].strip(), \
                                      site_meta[1].strip(), \
                                      site_meta[7].strip(), \
                                      site_meta[8].strip())
    
    SITES_temp = raw_input('Enter flux sites to plot seperated by commas, '+\
                           'for all sites enter ALL: ')
    if SITES_temp=='ALL':
        SITES = range(len(meta_data))
    elif ':' in SITES_temp:
        SITES_temp=SITES_temp.split(':')
        sSITE = int(SITES_temp[0])
        eSITE = int(SITES_temp[1])
        SITES = range(sSITE,eSITE+1)
        del sSITE
        del eSITE
    else:
        SITES = [int(site) for site in SITES_temp.split(',')]
    
    del SITES_temp
else:
    SITES_temp = sys.argv[4]
    if SITES_temp=='ALL':
        SITES=range(len(meta_data))
    elif ':' in SITES_temp:
        SITES_temp=SITES_temp.split(':')
        sSITE = int(SITES_temp[0])
        eSITE = int(SITES_temp[1])
        SITES = range(sSITE,eSITE+1)
        del sSITE
        del eSITE
    else:
        SITES_temp=filter(None,SITES_temp.replace('[','').replace(']','').split(','))
        SITES=[int(temp) for temp in SITES_temp]

nSITES=len(SITES)
print 'nSITES, SITES = ',nSITES,SITES

######################################################################################
# Select JULES source[s]
if (INTERACTIVE=='Y'):
    print 'Available JULES sources: '
    for i,Jsource in zip(range(len(JULES_sources)),JULES_sources):
        print '%3s:%25s' % (i,Jsource[0])

    source_temp = raw_input('Enter JULES sources to plot sepeated by commas, '+\
                            'for all sources enter ALL: ')
    if source_temp=='ALL':
        J_SOURCES = range(len(JULES_sources))
    else:
        J_SOURCES = [int(source) for source in source_temp.split(',')]

    del source_temp
else:
    J_SOURCES_temp = sys.argv[5]
    if J_SOURCES_temp=='ALL':
        J_SOURCES=range(len(JULES_sources))
    else:
        J_SOURCES_temp=filter(None,J_SOURCES_temp.replace('[','').replace(']','').split(','))
        J_SOURCES=[int(temp) for temp in J_SOURCES_temp]

nJSOURCES = len(J_SOURCES)
J_names = [JULES_sources[iSOURCE][1]['longname'] for iSOURCE in J_SOURCES]
print 'nJSOURCES, J_SOURCES = ',nJSOURCES,J_names

# construct list of inputs required for plotting later:
PS_names  = [ 'Site' ] + J_names 
PS_colours = [ SITE_colour ] + [JULES_sources[iSOURCE][1]['colour'] for iSOURCE in J_SOURCES]
nPS = len(PS_names)
PS_handles = [plt.Line2D((0,1),(0,0),color=PS_colours[iPS],lw=3) for iPS in range(nPS)] 

#######################################################################################
# T resolution
Tres_names = ['tstep','day']
Tres_Stags = ['','_daily']
if INTERACTIVE=='Y':
    print 'Available time resolutions: '
    for iTRES in range(len(Tres_names)):
        print iTRES, '- '+Tres_names[iTRES]
    iTRES = int(raw_input('Select a time resolution: '))
else:
    iTRES = int(sys.argv[6])

Tres = Tres_names[iTRES]
Tres_Stag= Tres_Stags[iTRES]

#######################################################################################
# Optional plot_dir input
plot_dir=BASE_plot_dir+Tres+'/'
if (INTERACTIVE=='Y'):
    print 'Plot output directory: '+plot_dir
    plot_dir_temp=raw_input('Enter alternative plot output directory or hit return to continue: ')
    if (len(plot_dir_temp)>0):
        plot_dir=plot_dir_temp
elif ('-plot_dir' in sys.argv):
    argloc=sys.argv.index('-plot_dir')
    plot_dir=sys.argv[argloc+1]

# Create plot_dir if doesn't already exist
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

########################################################################
# If INTERACTIVE, print resubmit command for ease next time:
if (INTERACTIVE=='Y'):
    print 'Resubmit Command: '
    #print './plot_PALS_sites_day_pandas.py N '+\
    print '\033[1;31m '+\
            __file__ + ' N ' + \
            iDISPLAY+' '+PLOTS+' '+\
            str(SITES).replace(' ','')+' '+\
            str(J_SOURCES).replace(' ','')+' '+\
            str(iTRES)+\
            '\033[0m'
    if (plot_dir!=BASE_plot_dir):
        print '-plot_dir '+plot_dir

################################################################################################
#
# Create empty lists for storing data
#####################################
SITE_lats     = []
SITE_lons     = []
SITE_PFT      = []
SITE_PFT_name = []
SITE_markers  = []

#create lists for storing data to plot for all sites
SITE_Net_Annual_Cuptake = []
SITE_Gross_Annual_Cuptake = []
SITE_Net_Annual_Evap = []
SITE_Net_Annual_Heat = []
JULES_Net_Annual_Evap = [ [] for i in range(nJSOURCES) ]
JULES_Net_Annual_Heat = [ [] for i in range(nJSOURCES) ]
JULES_Net_Annual_Cuptake = [ [] for i in range(nJSOURCES) ]
JULES_Gross_Annual_Cuptake = [ [] for i in range(nJSOURCES) ]

SITE_mean_diurnal_GPP_Amplitude = []
SITE_mean_diurnal_NEE_Amplitude = []
SITE_mean_diurnal_GPP_Amplitude_BySeas = [[] for i in range(4)]
SITE_mean_diurnal_NEE_Amplitude_BySeas = [[] for i in range(4)]
JULES_mean_diurnal_GPP_Amplitude = [ [] for i in range(nJSOURCES) ]
JULES_mean_diurnal_NEE_Amplitude = [ [] for i in range(nJSOURCES) ]
JULES_mean_diurnal_GPP_Amplitude_BySeas = [ [ [] for i in range(nJSOURCES) ] for j in range(4) ]
JULES_mean_diurnal_NEE_Amplitude_BySeas = [ [ [] for i in range(nJSOURCES) ] for j in range(4) ]

# Correlations and standard deviations for taylor plots
SHF_correlation = [ [] for i in range(nJSOURCES) ]
SHF_meanbias    = [ [] for i in range(nJSOURCES) ]
SHF_stddev      = [ [] for i in range(nJSOURCES) ]
LHF_correlation = [ [] for i in range(nJSOURCES) ]
LHF_meanbias    = [ [] for i in range(nJSOURCES) ]
LHF_stddev      = [ [] for i in range(nJSOURCES) ]
NEE_correlation = [ [] for i in range(nJSOURCES) ]
NEE_meanbias    = [ [] for i in range(nJSOURCES) ]
NEE_stddev      = [ [] for i in range(nJSOURCES) ]
GPP_correlation = [ [] for i in range(nJSOURCES) ]
GPP_meanbias    = [ [] for i in range(nJSOURCES) ]
GPP_stddev      = [ [] for i in range(nJSOURCES) ]

######################################################################################################
# Loop round selected site 
for cnt in range(nSITES):
    iSITE=SITES[cnt]
    temp_frac = frac_list[iSITE].split()
    site_meta=meta_data[iSITE]
    frac=np.array(temp_frac,dtype='float')
    site_name=site_meta[0].strip()
    site_country=site_meta[1].strip()
    site_vegtype=site_meta[2].strip()
    site_canht=site_meta[3].strip()
    site_measht=site_meta[4].strip()
    site_lat=site_meta[7].strip()
    site_lon=site_meta[8].strip()

    SITE_lats.append(float(site_lat))
    SITE_lons.append(float(site_lon))
    SITE_PFT_name.append(site_vegtype)
    SITE_PFT.append(temp_frac.index(max(temp_frac)))
    SITE_markers.append(LC_markers[temp_frac.index(max(temp_frac))])
    
    site_plot_dir = plot_dir+site_name+'/'
    os.system('mkdir -p '+site_plot_dir)

    print 'Site: '+site_name
    print 'Plots saved in: '+site_plot_dir
    
    ######################################
    # Read in-situ site data
    SITE_fname = SITE_data_dir+site_name+'Fluxnet.1.4_flux'+Tres_Stag+'.nc'
    print SITE_fname
    Sinf       = nc.Dataset(SITE_fname,'r')
    ############################################################
    # Extract data and convert to same units (Site data units)
    #'time':nctime.num2date(Sinf.variables['time'][:],units=Sinf.variables['time'].units),       \
    S_data = {'NEE':Sinf.variables['NEE'][:].squeeze()*-1.0,                    \
              'TER':Sinf.variables['TER'][:].squeeze(), \
              'GPP':Sinf.variables['GPP'][:].squeeze(), \
              'LHF':Sinf.variables['Qle'][:].squeeze(),                         \
              'SHF':Sinf.variables['Qh'][:].squeeze(),                          \
              }

    #         'TER':np.ma.masked_less(Sinf.variables['TER'][:].squeeze(),0.0), \
    #         'GPP':np.ma.masked_less(Sinf.variables['GPP'][:].squeeze(),0.0), \

    S_seconds=Sinf.variables['time'][:].astype('float64')
    S_seconds=np.round(S_seconds/1800.)*1800.
    S_time = nctime.num2date(S_seconds, \
                             units=Sinf.variables['time'].units)
    S_data['FQW'] = S_data['LHF']/Lc_H2O
    S_data['WUE'] = S_data['GPP']/S_data['FQW']
    
    S_panda = pd.DataFrame( S_data, index=S_time, )
    
    # Store Annual C uptake in a list for later plots
    SITE_Net_Annual_Evap.append( S_panda['FQW'].mean()*seconds_in_year )
    SITE_Net_Annual_Heat.append( S_panda['SHF'].mean()*seconds_in_year )

    SITE_Net_Annual_Cuptake.append( (S_panda['NEE'].mean()/kgC_to_umolsCO2_factor)\
                                    * seconds_in_year )
    print SITE_Net_Annual_Cuptake[-1]
    SITE_Gross_Annual_Cuptake.append( (S_panda['GPP'].mean()/kgC_to_umolsCO2_factor)\
                                    * seconds_in_year )
    print SITE_Gross_Annual_Cuptake[-1]
    Sinf.close()

    ############################################
    # Loop round JULES sources and extract data to dictionary of lists
    J_pandas=[]

    for Jcnt in range(nJSOURCES):
        iSOURCE = J_SOURCES[Jcnt]
        #print iSOURCE, Jcnt, len(JULES_sources), nJSOURCES
        Jsource = JULES_sources[iSOURCE]
        ######################################
        # Construct filename and open:
        JULES_fname=Jsource[1]['data_dir']+site_name+'.'+Tres+'.nc'
        print JULES_fname
        Jinf = nc.Dataset(JULES_fname,'r')
        J_data = { 'NPP':Jinf.variables[Jsource[1]['npp_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor,              \
                   'GPP':np.ma.masked_less(Jinf.variables[Jsource[1]['gpp_name']][:].squeeze() \
                                  *kgC_to_umolsCO2_factor, 1e-5 ),               \
                   'resp_s':Jinf.variables[Jsource[1]['resp_s_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor,                   \
                   'resp_p':Jinf.variables[Jsource[1]['resp_p_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor,                   \
                   'FQW':Jinf.variables['fqw_gb'][:].squeeze(),                 \
                   'LHF':Jinf.variables['latent_heat'][:].squeeze(),            \
                   'SHF':Jinf.variables['ftl_gb'][:].squeeze(),                 \
                   }
        
        # Quick fix to add diurnal cycle to J-E-F output
        if (Tres=='tstep') & ('-E-' in Jsource[1]['longname']):
            print 'Calculating Diurnal NPP from original NPP and the ratio of '+\
                  ' npp_nuptake_out/npp_nuptake_out' 
            J_data['NPP']= Jinf.variables['npp_gb'][:].squeeze() * \
                           ( Jinf.variables['npp_nuptake_out_gb'][:].squeeze() / \
                             Jinf.variables['npp_nuptake_in_gb'][:].squeeze()  ) \
                                  *   kgC_to_umolsCO2_factor
        
        J_data['NEE']= J_data['NPP']-J_data['resp_s'] 
        J_data['WUE']= J_data['GPP']/J_data['FQW']
        J_seconds=Jinf.variables['time_bounds'][:,0].astype('float64')
        J_seconds=np.round(J_seconds/1800.)*1800.
        J_time=nctime.num2date(J_seconds, units=Jinf.variables['time'].units  )
         
        J_panda = pd.DataFrame(J_data,index=J_time)
        J_pandas.append(J_panda.copy())
        
        Jinf.close()
        
        # Calculate Annual uptakes/evaporation/energ flux
        JULES_Net_Annual_Evap[Jcnt].append( J_panda['FQW'].mean()*seconds_in_year )
        JULES_Net_Annual_Heat[Jcnt].append( J_panda['SHF'].mean()*seconds_in_year )
        
        JULES_Net_Annual_Cuptake[Jcnt].append(J_panda['NEE'].mean()*   \
                                              seconds_in_year/kgC_to_umolsCO2_factor )
        JULES_Gross_Annual_Cuptake[Jcnt].append(J_panda['GPP'].mean()* \
                                                seconds_in_year/kgC_to_umolsCO2_factor )
        
        # Calculate and store the Specific Heat Flux Statistics  
        SHF_meanbias[Jcnt].append((J_panda['SHF']-S_panda['SHF']).mean()) 
        SHF_stddev[Jcnt].append( (J_panda['SHF']-S_panda['SHF']).std() )
        SHF_correlation[Jcnt].append( J_panda['SHF'].corr(S_panda['SHF']))
        
        # Calculate and store the Latent Heat Flux Statistics  
        LHF_meanbias[Jcnt].append((J_panda['LHF']-S_panda['LHF']).mean()) 
        LHF_stddev[Jcnt].append( (J_panda['LHF']-S_panda['LHF']).std() )
        LHF_correlation[Jcnt].append( J_panda['LHF'].corr(S_panda['LHF']))
        
        # Calculate and store the NEE Statistics  
        NEE_meanbias[Jcnt].append((J_panda['NEE']-S_panda['NEE']).mean()) 
        NEE_stddev[Jcnt].append( (J_panda['NEE']-S_panda['NEE']).std() )
        NEE_correlation[Jcnt].append( J_panda['NEE'].corr(S_panda['NEE']))
        
        # Calculate and store the GPP Statistics 
        print (J_panda['GPP']-S_panda['GPP']).mean()
        GPP_meanbias[Jcnt].append((J_panda['GPP']-S_panda['GPP']).mean()) 
        print (J_panda['GPP']-S_panda['GPP']).std() 
        GPP_stddev[Jcnt].append( (J_panda['GPP']-S_panda['GPP']).std() )
        print  J_panda['GPP'].corr(S_panda['GPP'])
        GPP_correlation[Jcnt].append( J_panda['GPP'].corr(S_panda['GPP']))
        
        
    #Create Parameter Pandas
    LHF_panda=pd.concat([S_panda['LHF']]+[J_panda['LHF'] for J_panda in J_pandas],axis=1)
    LHF_panda.columns=(PS_names)
    #LHF_panda=LHF_panda.dropna(0)
    SHF_panda=pd.concat([S_panda['SHF']]+[J_panda['SHF'] for J_panda in J_pandas],axis=1)
    SHF_panda.columns=(PS_names)
    #SHF_panda=SHF_panda.dropna(0)
    NEE_panda=pd.concat([S_panda['NEE']]+[J_panda['NEE'] for J_panda in J_pandas],axis=1)
    NEE_panda.columns=(PS_names)
    #NEE_panda=NEE_panda.dropna(0)
    GPP_panda=pd.concat([S_panda['GPP']]+[J_panda['GPP'] for J_panda in J_pandas],axis=1)
    GPP_panda.columns=(PS_names)
    #GPP_panda=GPP_panda.dropna(0)
    NPP_panda=pd.concat([J_panda['NPP'] for J_panda in J_pandas],axis=1)
    NPP_panda.columns=(PS_names[1:])
    #NPP_panda=NPP_panda.dropna(0)
    RESPS_panda=pd.concat([J_panda['resp_s'] for J_panda in J_pandas],axis=1)
    RESPS_panda.columns=(PS_names[1:])
    #RESPS_panda=RESPS_panda.dropna(0)
    RESPP_panda=pd.concat([J_panda['resp_p'] for J_panda in J_pandas],axis=1)
    RESPP_panda.columns=(PS_names[1:])
    #RESPP_panda=RESPP_panda.dropna(0)
    TER_panda=pd.concat([S_panda['TER']],axis=1)
    TER_panda.columns=(PS_names[0:1])
    

    #########################################################################################
    ##########################################################
    # Plot time-series of SHF, LHF, and NEE
    if (PLOTS[0]=='Y'):

        FIG = plt.figure(figsize=(14,12))
        
        AX = FIG.add_subplot(3,1,1)
        SHF_panda.plot(lw=1.5,ax=AX,color=PS_colours,legend=False)
        AX.text(1.02,0.22,"Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.08,"SHF Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            AX.text(1.02+(Jcnt*0.07),0.16,"%6.2f"%(SHF_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.07),0.02,"%6.2f"%(SHF_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
        AX.set_ylabel('Sensible Heat Flux (W $m^{-2}$)', fontsize=14)
        
        AX = FIG.add_subplot(3,1,2)
        LHF_panda.plot(lw=1.5,ax=AX,color=PS_colours,legend=False)
        AX.text(1.02,0.22,"LHF Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.08,"LHF Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            AX.text(1.02+(Jcnt*0.07),0.16,"%6.2f"%(LHF_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.07),0.02,"%6.2f"%(LHF_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
        AX.set_ylabel('Latent Heat Flux (W $m^{-2}$)', fontsize=14)
        handles,labels=AX.get_legend_handles_labels()
        
        AX = FIG.add_subplot(3,1,3)
        NEE_panda.plot(lw=1.5,ax=AX,color=PS_colours,legend=False)
        GPP_panda.plot(lw=1.5,ax=AX,ls='-.',color=PS_colours,legend=False)
        AX.text(1.02,0.46,"NEE Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.32,"NEE Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.16,"GPP Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.02,"GPP Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            AX.text(1.02+(Jcnt*0.07),0.38,"%6.2f"%(NEE_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.07),0.26,"%6.2f"%(NEE_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.07),0.10,"%6.2f"%(GPP_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.07),-0.04,"%6.2f"%(GPP_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.075+(Jcnt*0.06),0.71,"%6.2f"%(JULES_Net_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.075+(Jcnt*0.06),0.87,"%6.2f"%(JULES_Gross_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )

        AX.set_ylabel('NEE and GPP ($\mu$mol $m^{-2} s^{-1}$)')

        AX.text(1.03,0.95,"Gross C Uptake: ",\
                transform=AX.transAxes, fontsize=14)
        AX.text(1.01,0.87,"%6.2f"%(SITE_Gross_Annual_Cuptake[cnt]),\
                transform=AX.transAxes, fontsize=14, color=SITE_colour )
        AX.text(1.03,0.79,"Net C Uptake:",\
                transform=AX.transAxes, fontsize=14)
        AX.text(1.01,0.71,"%6.2f"%(SITE_Net_Annual_Cuptake[cnt]),\
                transform=AX.transAxes, fontsize=14, color=SITE_colour )
        AX.text(1.03,0.63,"(gC y$^{-1}$ m$^{-2}$)",\
                transform=AX.transAxes, fontsize=14)
        
        FIG.legend( handles,labels, loc=8, ncol=nJSOURCES+1)
        FIG.tight_layout(rect=[0,0.05,0.8,0.96])
        FIG.suptitle(site_name+'flux tower time-series, obs and J-E-F',\
                     fontsize=24)
        FIG.text(0.8,0.92,'Country: '+site_country,fontsize=16)
        FIG.text(0.8,0.895,'PFT: '+site_vegtype,fontsize=16)
        FIG.text(0.8,0.87,'Lat: '+site_lat+'$^o$E',fontsize=16)
        FIG.text(0.8,0.845,'Lon: '+site_lon+'$^o$N',fontsize=16)
        FIG.text(0.8,0.82,'Can. hgt: '+site_canht+' m',fontsize=16)
        FIG.text(0.8,0.795,'Meas. hgt: '+site_measht+' m',fontsize=16)
        
        
        FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_Time-Series.png', \
                    bbox_inches='tight')
    
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
            
             
    ##########################################################
    # Plot scatter plots of LHF and SHF data
        
    if (PLOTS[1]=='Y'):
        FIG = plt.figure(figsize=(18,9))
        
        # Plot the LHF data
        AX = FIG.add_subplot(1,2,1)
                
        for Jcnt in range(nJSOURCES):
            Jname = J_names[Jcnt] 
            AX.scatter( LHF_panda['Site'],LHF_panda[Jname], \
                        c=PS_colours[Jcnt+1], label=Jname,linewidth=0 )
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(LHF_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        
        AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
                fontsize=18, color='k' )
        limits=[min(AX.get_xlim()[0],AX.get_ylim()[0]),\
                max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
        AX.plot(np.array(limits),np.array(limits),c='k')
        AX.set_title('Latent Heat Flux',fontsize=14)
        AX.set_xlabel('Site (W $m^{-2}$)')
        AX.set_xlim(limits)
        AX.set_ylabel('JULES (W $m^{-2}$)')
        AX.set_ylim(limits)
        AX.grid(True)
        AX.legend( bbox_to_anchor=(1.1,-0.1),loc=8,borderaxespad=0.,ncol=nJSOURCES+1 )

        # Plot the SHF data
        AX = FIG.add_subplot(1,2,2)
        for Jcnt in range(nJSOURCES):
            Jname = J_names[Jcnt] 
            AX.scatter( SHF_panda['Site'],LHF_panda[Jname], \
                        c=PS_colours[Jcnt+1], label=Jname,linewidth=0 )
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(SHF_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        
        AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
                fontsize=18, color='k' )
        limits=[min(AX.get_xlim()[0],AX.get_ylim()[0]),\
                max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
        AX.plot(np.array(limits),np.array(limits),c='k')
        AX.set_title('Sensible Heat Flux',fontsize=14)
        AX.set_xlabel('Site (W $m^{-2}$)')
        AX.set_xlim(limits)
        AX.set_ylabel('JULES (W $m^{-2}$)')
        AX.set_ylim(limits)
        AX.grid(True)
        
        FIG.tight_layout(rect=[0,0.02,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, Heat Flux Comparison',\
                     fontsize=24)
        FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_HFscatter.png', \
                    bbox_inches='tight')
        
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
        
    ##########################################################
    # Plot scatter plot of NEE and GPP data
        
    if (PLOTS[2]=='Y'):
        FIG = plt.figure(figsize=(18,9))
        
        # Plot the NEE data
        AX = FIG.add_subplot(1,2,1)
        for Jcnt in range(nJSOURCES):
            Jname = J_names[Jcnt] 
            AX.scatter( NEE_panda['Site'],NEE_panda[Jname], \
                        c=PS_colours[Jcnt+1], label=Jname,linewidth=0 )
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(NEE_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        
        AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
                fontsize=18, color='k' )
        limits=[min(AX.get_xlim()[0],AX.get_ylim()[0]),\
                max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
        AX.plot(np.array(limits),np.array(limits),c='k')
        AX.set_title('Net Ecosystem Exchange',fontsize=14)
        AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_xlim(limits)
        AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylim(limits)
        AX.grid(True)
        AX.legend( bbox_to_anchor=(1.1,-0.1),loc=8,borderaxespad=0.,ncol=nJSOURCES+1 )

        # Plot the GPP data
        AX = FIG.add_subplot(1,2,2)
        for Jcnt in range(nJSOURCES):
            Jname = J_names[Jcnt] 
            AX.scatter( GPP_panda['Site'],GPP_panda[Jname], \
                        c=PS_colours[Jcnt+1], label=Jname ,linewidth=0)
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(GPP_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        
        AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
                fontsize=18, color='k' )
        limits=[min(AX.get_xlim()[0],AX.get_ylim()[0]),\
                max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
        AX.plot(np.array(limits),np.array(limits),c='k')
        AX.set_title('Gross Primary Productivity',fontsize=14)
        AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_xlim(limits)
        AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylim(limits)
        AX.grid(True)
        FIG.tight_layout(rect=[0,0.02,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, Carbon Flux Comparison',\
                     fontsize=24)
        FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_NEEGPPscatter.png', \
                    bbox_inches='tight')
        
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
               
    ##########################################################
    # Plot scatter plot of LHF against NEE
        
    if (PLOTS[3]=='Y'):
        FIG = plt.figure(figsize=(12,12))
            
        AX = FIG.add_subplot(1,1,1)
        J_slope=[]
        J_intercept=[]
        J_rvals=[]
        for Jcnt in range(nJSOURCES):
            AX.scatter( J_pandas[Jcnt]['LHF'],J_pandas[Jcnt]['GPP'],\
                         c=PS_colours[Jcnt+1], linewidth=0 )
            slope, intercept, r, p, stderr = stats.linregress(J_pandas[Jcnt].dropna(0)['LHF'],\
                                                              J_pandas[Jcnt].dropna(0)['GPP'] )
            J_slope.append(slope)
            J_intercept.append(intercept)
            J_rvals.append(r)
        
        AX.scatter(S_panda['LHF'],S_panda['GPP'],c=PS_colours[0],marker='x')
        S_slope, S_intercept, S_rval, p, stderr = stats.linregress(S_panda.dropna(0)['LHF'],\
                                                                   S_panda.dropna(0)['GPP'] )

        xlims = [0,AX.get_xlim()[1]]
        ylims = [0,AX.get_ylim()[1]]
        
        for Jcnt in range(nJSOURCES):
            AX.plot( np.array(xlims), (np.array(xlims)*J_slope[Jcnt])+J_intercept[Jcnt],   \
                     c=PS_colours[Jcnt+1],linewidth=2,                 \
                     label=J_names[Jcnt]+', $y=%4.2fx+%4.2f$, $r^2= %4.2f$' % \
                           (J_slope[Jcnt],J_intercept[Jcnt],J_rvals[Jcnt])  )
        
        AX.plot( np.array(xlims), (np.array(xlims)*S_slope)+S_intercept,   \
                 c=PS_colours[0],linewidth=2,                 \
                 label='Site, $y=%4.2fx+%4.2f$, $r^2= %4.2f$'%(S_slope,S_intercept,S_rval) )
        
        AX.set_xlabel('Latent Heat Flux (W $m^{-2}$)')
        AX.set_xlim(xlims)
        AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylim(ylims)
        AX.legend(loc='upper left')
        
        FIG.tight_layout(rect=[0,0,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, LHF vs GPP',\
                     fontsize=24)
        FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_LHFvGPPscatter.png', \
                    bbox_inches='tight')
        
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
            
    #########################################################################################
    ##########################################################
    # Plot time-series breakdown of CO2 fluxes
    if (PLOTS[4]=='Y'):
        FIG,AXES = plt.subplots(nrows=2,ncols=1,figsize=(18,10))
        
        AX = AXES[0]
        #print GPP_panda
        for Jcnt in range(nJSOURCES):
            Jname = J_names[Jcnt] 
            AX.plot(GPP_panda.index,GPP_panda[Jname],lw=1.5,color='g',label='GPP') 
            AX.plot(NPP_panda.index,NPP_panda[Jname],lw=1.5,color='greenyellow',label='NPP')
            AX.plot(RESPS_panda.index,RESPS_panda[Jname]*-1.,lw=1.5,color='saddlebrown',label='Soil Resp.')
            AX.plot(RESPP_panda.index,RESPP_panda[Jname]*-1.,lw=1.5,color='orange',label='Plant Resp.')

        AX.grid(b=True)
        AX.set_ylabel('CO2 flux ($\mu$mol $m^{-2} s^{-1}$)')
        handles,labels=AX.get_legend_handles_labels()
        
        AX=AXES[1]
        NEE_panda.plot(lw=2,ax=AX,color=PS_colours,legend=False)
        AX.set_ylabel('CO2 flux ($\mu$mol $m^{-2} s^{-1}$)')
        handles2,labels2=AX.get_legend_handles_labels()
        handles=handles+handles2
        labels=labels+labels2
        
        FIG.legend( handles,labels, loc=8, ncol=6)
        FIG.tight_layout(rect=[0,0.05,0.96,0.96])
        FIG.suptitle(site_name+'flux tower time-series, obs and J-E-F',\
                     fontsize=24)
        FIG.savefig(site_plot_dir+site_name+'_CO2_Flux_Breakdown_Time-Series.png', \
                    bbox_inches='tight')
    
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
            
             
    ##########################################################
    #########################################################################################
    # Plot Diurnal Behaviour
    #
    # Full Simulation
    if ((PLOTS[10]=='Y')|(PLOTS[12]=='Y')|(PLOTS[14]=='Y')) & (Tres=='tstep'):
        # Create Diurnal Panda Data Frames for NEE and GPP
        NEE_diurnal_panda=NEE_panda.dropna(0).copy()
        NEE_diurnal_panda['Time']=NEE_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        NEE_diurnal_panda = NEE_diurnal_panda.groupby('Time').describe().unstack()
        NEE_diurnal_panda.index = pd.to_datetime(NEE_diurnal_panda.index.astype(str))

        GPP_diurnal_panda=GPP_panda.dropna(0).copy()
        GPP_diurnal_panda['Time']=GPP_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        GPP_diurnal_panda = GPP_diurnal_panda.groupby('Time').describe().unstack()
        GPP_diurnal_panda.index = pd.to_datetime(GPP_diurnal_panda.index.astype(str))

        NPP_diurnal_panda=NPP_panda.dropna(0).copy()
        NPP_diurnal_panda['Time']=NPP_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        NPP_diurnal_panda = NPP_diurnal_panda.groupby('Time').describe().unstack()
        NPP_diurnal_panda.index = pd.to_datetime(NPP_diurnal_panda.index.astype(str))

        RESPP_diurnal_panda=RESPP_panda.dropna(0).copy()
        RESPP_diurnal_panda['Time']=RESPP_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        RESPP_diurnal_panda = RESPP_diurnal_panda.groupby('Time').describe().unstack()
        RESPP_diurnal_panda.index = pd.to_datetime(RESPP_diurnal_panda.index.astype(str))

        RESPS_diurnal_panda=RESPS_panda.dropna(0).copy()
        RESPS_diurnal_panda['Time']=RESPS_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        RESPS_diurnal_panda = RESPS_diurnal_panda.groupby('Time').describe().unstack()
        RESPS_diurnal_panda.index = pd.to_datetime(RESPS_diurnal_panda.index.astype(str))

        TER_diurnal_panda=TER_panda.dropna(0).copy()
        TER_diurnal_panda['Time']=TER_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        TER_diurnal_panda = TER_diurnal_panda.groupby('Time').describe().unstack()
        TER_diurnal_panda.index = pd.to_datetime(TER_diurnal_panda.index.astype(str))
                
        # Store mean amplitude of GPP and NEE for multi-site stats
        if (PLOTS[12]=='Y'):
            SITE_mean_diurnal_GPP_Amplitude.append(GPP_diurnal_panda['Site']['mean'].max() - \
                                                   GPP_diurnal_panda['Site']['mean'].min()   )
            SITE_mean_diurnal_NEE_Amplitude.append(NEE_diurnal_panda['Site']['mean'].max() - \
                                                   NEE_diurnal_panda['Site']['mean'].min()   )
            
            for Jcnt in range(nJSOURCES):
                JULES_mean_diurnal_GPP_Amplitude[Jcnt].append( \
                                    GPP_diurnal_panda[PS_names[Jcnt+1]]['mean'].max() - \
                                    GPP_diurnal_panda[PS_names[Jcnt+1]]['mean'].min()   )
                JULES_mean_diurnal_NEE_Amplitude[Jcnt].append( \
                                    NEE_diurnal_panda[PS_names[Jcnt+1]]['mean'].max() - \
                                    NEE_diurnal_panda[PS_names[Jcnt+1]]['mean'].min()   )
                    
        
        if (PLOTS[10]=='Y'):
            # Create Figure and Axes for diurnal plots
            FIG, AXES = plt.subplots(nrows=2,ncols=1,figsize=[12,9])
            # Plot GPP in upper plot
            for iPS in range(nPS):
                AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[PS_names[iPS]]['mean'], \
                             c=PS_colours[iPS],ls='-',lw=2,label=PS_names[iPS])
                AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[PS_names[iPS]]['50%'], \
                             c=PS_colours[iPS],ls='-.',lw=1.5)
                AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[PS_names[iPS]]['25%'], \
                             c=PS_colours[iPS],ls=':')
                AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[PS_names[iPS]]['75%'], \
                             c=PS_colours[iPS],ls=':')

            AXES[0].set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)',fontsize=14)
            AXES[0].set_ylim([0,AXES[0].get_ylim()[1]])
            AXES[0].grid(True)

            # Plot GPP in upper plot
            for iPS in range(nPS):
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['mean'], \
                             c=PS_colours[iPS],ls='-',lw=2)
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['50%'], \
                             c=PS_colours[iPS],ls='-.',lw=1.5)
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['25%'], \
                             c=PS_colours[iPS],ls=':')
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['75%'], \
                             c=PS_colours[iPS],ls=':')
        
            AXES[1].set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)',fontsize=14)
            NEE_lim = np.max( np.absolute(AXES[1].get_ylim()) ) 
            AXES[1].set_ylim([-NEE_lim,NEE_lim])
            AXES[1].grid(True)
            
            handles,labels=AXES[0].get_legend_handles_labels()
            FIG.legend(handles,labels,loc=8,ncol=nPS)
            FIG.suptitle(site_name+'('+site_country+') flux tower, GPP and NEE diurnal cycle',\
                         fontsize=24)
            
            FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_GPPNEEdiurnalcycle.png', \
                        bbox_inches='tight')
        
            if iDISPLAY=='Y':
                plt.show()
            else:
                plt.close()

        if (PLOTS[14]=='Y'):
            # Create Figure and Axes for diurnal plots
            FIG, AXES = plt.subplots(nrows=2,ncols=1,figsize=[12,9])
            
            #for Jcnt in range(nJSOURCES):
            Jcnt=0
            Jname = J_names[Jcnt] 
            AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[PS_names[0]]['mean'], \
                         c='b',lw=1.5,label='GPP-SITE')
            AXES[0].plot(TER_diurnal_panda.index, TER_diurnal_panda[PS_names[0]]['mean'], \
                         c='r',lw=1.5,label='TER-SITE')
            AXES[0].plot(GPP_diurnal_panda.index, GPP_diurnal_panda[Jname]['mean'], \
                         c='g',lw=1.5,label='GPP-JULES')
            AXES[0].plot(NPP_diurnal_panda.index, NPP_diurnal_panda[Jname]['mean'], \
                         c='greenyellow',lw=1.5,label='NPP')
            AXES[0].plot(RESPP_diurnal_panda.index, RESPP_diurnal_panda[Jname]['mean'], \
                         c='orange',lw=1.5,label='Plant Resp.')
            AXES[0].plot(RESPS_diurnal_panda.index, RESPS_diurnal_panda[Jname]['mean'], \
                         c='saddlebrown',lw=1.5,label='Soil Resp.')

            AXES[0].legend(loc=8,bbox_to_anchor=(0.5,-0.2),ncol=6)
            AXES[0].set_ylabel('Carbon Flux ($\mu$mol $m^{-2} s^{-2}$)',fontsize=14)
            #AXES[0].set_ylim([0,AXES[0].get_ylim()[1]])
            AXES[0].grid(True)

            # Plot GPP in upper plot
            for iPS in range(2):
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['mean'], \
                             c=PS_colours[iPS],ls='-',lw=2,label=PS_names[iPS])
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['50%'], \
                             c=PS_colours[iPS],ls='-.',lw=1.5)
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['25%'], \
                             c=PS_colours[iPS],ls=':')
                AXES[1].plot(NEE_diurnal_panda.index, NEE_diurnal_panda[PS_names[iPS]]['75%'], \
                             c=PS_colours[iPS],ls=':')
        
            AXES[1].set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)',fontsize=14)
            NEE_lim = np.max( np.absolute(AXES[1].get_ylim()) ) 
            AXES[1].set_ylim([-NEE_lim,NEE_lim])
            AXES[1].grid(True)
            
            FIG.suptitle(site_name+'('+site_country+') flux tower, Carbon Flux Diurnal Cycle, '+Jname,\
                         fontsize=24)
            
            handles,labels=AXES[1].get_legend_handles_labels()
            FIG.legend(handles,labels,loc=8,ncol=nPS)

            FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_GPPNEERESPdiurnalcycle.png', \
                        bbox_inches='tight')
        
            if iDISPLAY=='Y':
                plt.show()
            else:
                plt.close()

    # Seasonal
    if ((PLOTS[11]=='Y')|(PLOTS[13]=='Y')) & (Tres=='tstep'):
        # Create Diurnal Panda Data Frames for NEE and GPP
        NEE_diurnal_panda=NEE_panda.dropna(0).copy()
        NEE_diurnal_panda['Time']=NEE_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        NEE_diurnal_panda['Seas']=NEE_diurnal_panda.index.month/4
        NEE_diurnal_panda = NEE_diurnal_panda.groupby('Seas')

        GPP_diurnal_panda=GPP_panda.dropna(0).copy()
        GPP_diurnal_panda['Time']=GPP_diurnal_panda.index.map(lambda x: x.strftime("%H:%M"))
        GPP_diurnal_panda['Seas']=GPP_diurnal_panda.index.month/4
        GPP_diurnal_panda = GPP_diurnal_panda.groupby('Seas')

        NEE_seasdiur_panda_list = []
        GPP_seasdiur_panda_list = []

        for i in range(nSEAS):
            NEE_seasdiur_panda_list.append( \
                    NEE_diurnal_panda.get_group(i).groupby('Time').describe().unstack()  )
            NEE_seasdiur_panda_list[i].index = pd.to_datetime(NEE_seasdiur_panda_list[i].index.astype(str))

            GPP_seasdiur_panda_list.append( \
                    GPP_diurnal_panda.get_group(i).groupby('Time').describe().unstack()  )
            GPP_seasdiur_panda_list[i].index = pd.to_datetime(GPP_seasdiur_panda_list[i].index.astype(str))

            # Store mean amplitude of GPP and NEE for multi-site stats
            if (PLOTS[13]=='Y'):
                SITE_mean_diurnal_GPP_Amplitude_BySeas[i].append( \
                        GPP_seasdiur_panda_list[i]['Site']['mean'].max() - \
                        GPP_seasdiur_panda_list[i]['Site']['mean'].min()   )
                SITE_mean_diurnal_NEE_Amplitude_BySeas[i].append( \
                        NEE_seasdiur_panda_list[i]['Site']['mean'].max() - \
                        NEE_seasdiur_panda_list[i]['Site']['mean'].min()   )
              
                for Jcnt in range(nJSOURCES):
                    JULES_mean_diurnal_GPP_Amplitude_BySeas[i][Jcnt].append( \
                            GPP_seasdiur_panda_list[i][PS_names[Jcnt+1]]['mean'].max() - \
                            GPP_seasdiur_panda_list[i][PS_names[Jcnt+1]]['mean'].min()   )
                    JULES_mean_diurnal_NEE_Amplitude_BySeas[i][Jcnt].append( \
                            NEE_seasdiur_panda_list[i][PS_names[Jcnt+1]]['mean'].max() - \
                            NEE_seasdiur_panda_list[i][PS_names[Jcnt+1]]['mean'].min()   )

        if (PLOTS[11]=='Y'):
            # Group info for simpler plotting
            panda_lists = [GPP_seasdiur_panda_list,NEE_seasdiur_panda_list]
            var_names = ['GPP','NEE']
            nVARS = len(var_names)
            max_limits = [0 for i in range(nVARS)]
            units='($\mu$mol $m^{-2} s^{-2}$)'
            
            # Create Figure and Axes for diurnal plots
            FIG, AXES = plt.subplots(nrows=2,ncols=4,figsize=[24,7]) 
            # Plot GPP in upper plot
            for iVAR in range(nVARS):
                for iSEAS in range(nSEAS):
                    iAXES = (iVAR*nVARS)+iSEAS+1
                    panda_df = panda_lists[iVAR][iSEAS]
                    
                    for iPS in range(nPS):
                        AXES[iVAR,iSEAS].plot(panda_df.index, panda_df[PS_names[iPS]]['mean'], \
                                              c=PS_colours[iPS],ls='-',lw=2,label=PS_names[iPS])
                        AXES[iVAR,iSEAS].plot(panda_df.index, panda_df[PS_names[iPS]]['50%'], \
                                              c=PS_colours[iPS],ls='-.',lw=1.5)
                        AXES[iVAR,iSEAS].plot(panda_df.index, panda_df[PS_names[iPS]]['25%'], \
                                              c=PS_colours[iPS],ls=':',lw=1.2)
                        AXES[iVAR,iSEAS].plot(panda_df.index, panda_df[PS_names[iPS]]['75%'], \
                                              c=PS_colours[iPS],ls=':',lw=1.2)
                        
                    if iVAR==0:
                        AXES[iVAR,iSEAS].set_title(SEAS_names[iSEAS],fontsize=17)
                    
                    if iSEAS==0:
                        AXES[iVAR,iSEAS].set_ylabel(var_names[iVAR]+' '+units,fontsize=14)
                    
                    xticks=AXES[iVAR,iSEAS].get_xticks()[::2]
                    xticks=np.append(xticks,xticks[0]+1)
                    AXES[iVAR,iSEAS].set_xticks(xticks)
                    AXES[iVAR,iSEAS].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
                    AXES[iVAR,iSEAS].grid(True)
                    
                    # Store max ylimits for both GPP and NEE seperately so all seasons plotted
                    # on the same scale
                    max_limits[iVAR] = max([max_limits[iVAR]]+list(np.absolute(AXES[iVAR,iSEAS].get_ylim())))

            for iSEAS in range(nSEAS):
                AXES[0,iSEAS].set_ylim( [0,max_limits[0]] )
                AXES[1,iSEAS].set_ylim([-max_limits[1],max_limits[1]])

               
            handles,labels=AXES[0,0].get_legend_handles_labels()
            FIG.legend(handles,labels,loc=8,ncol=nPS)
            FIG.suptitle(site_name+'('+site_country+') flux tower, GPP and NEE diurnal cycle by season',\
                         fontsize=24)
            FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_GPPNEEdiurnalcycleBySeason.png', \
                        bbox_inches='tight')
            if iDISPLAY=='Y':
                plt.show()
            else:
                plt.close()

#########################################################################################
GRASS_TREE_PFT = np.array(SITE_PFT)
GRASS_TREE_PFT[GRASS_TREE_PFT==1]=0
unique_PFTs = list(set(GRASS_TREE_PFT))

#unique_PFTs = list(set(SITE_PFT))
# plot JULES vs SITE Cuptake (GROSS and NET)
# Scatter and Taylor
#  2x2 plots, columns: GPP, NEE
#             rows: scatter, taylor
if (PLOTS[5]=='Y'):
    # Create nice big figure:
    FIG=plt.figure(figsize=(18,12))
    # Plot Scatter in upper 2/3 of figure and create lists of arrays for taylor plots
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,0), rowspan=2))
    Taylor_std_list=[]
    Taylor_corr_list=[]
    Taylor_col_list=[]
    Taylor_marker_list=[]
    Taylor_lw_list=[]
    Taylor_s_list=[]
    print 'Gross C uptake stats:'
    print 5*'%15s' % ('Jsource','PFT','mean bias','std dev','correlation')
    for Jcnt in range(nJSOURCES):
        for iPFT in unique_PFTs: 
            #index = np.array(SITE_PFT)==iPFT
            index = GRASS_TREE_PFT==iPFT
            Jname=J_names[Jcnt]
            #print index
            print '%15s%15d%15.3f%15.3f%15.3f' % \
                (Jname, iPFT+1, \
                 np.mean(  np.array(JULES_Gross_Annual_Cuptake[Jcnt])[index]        \
                         - np.array(SITE_Gross_Annual_Cuptake)[index] ),   \
                 np.std(  np.array(JULES_Gross_Annual_Cuptake[Jcnt])[index]        \
                        - np.array(SITE_Gross_Annual_Cuptake)[index] ),   \
                 stats.pearsonr(np.array(JULES_Gross_Annual_Cuptake[Jcnt])[index], \
                                np.array(SITE_Gross_Annual_Cuptake)[index] )[0] \
                 )
                 
            Taylor_std_list.append(np.array(GPP_stddev[Jcnt])[index])
            Taylor_corr_list.append(np.array(GPP_correlation[Jcnt])[index])
            Taylor_col_list.append(PS_colours[Jcnt+1])
            Taylor_marker_list.append(LC_markers[iPFT])
            Taylor_lw_list.append(LC_mews[iPFT])
            Taylor_s_list.append(LC_symsizes[iPFT]/2.)
            AX.scatter( np.array(SITE_Gross_Annual_Cuptake)[index],  \
                        np.array(JULES_Gross_Annual_Cuptake)[Jcnt][index], \
                        c=PS_colours[Jcnt+1],\
                        marker=LC_markers[iPFT], s=LC_symsizes[iPFT],linewidth=LC_mews[iPFT] )
       
        correlation=stats.pearsonr(np.array(SITE_Gross_Annual_Cuptake),      \
                                   np.array(JULES_Gross_Annual_Cuptake[Jcnt]))[0]
        AX.text(0.1+(Jcnt*0.1),0.92,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        mean_bias = np.mean(  np.array(JULES_Gross_Annual_Cuptake)[Jcnt][:]  \
                             -np.array(SITE_Gross_Annual_Cuptake)[:]  )
        AX.text(0.1+(Jcnt*0.1),0.84,"%6.2f"%(mean_bias),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        stddev = np.std(  np.array(SITE_Gross_Annual_Cuptake)[:]  \
                        - np.array(JULES_Gross_Annual_Cuptake)[Jcnt][:])
        AX.text(0.1+(Jcnt*0.1),0.76,"%6.3f"%(stddev),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        
    AX.text(0.05,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.87,"Mean Bias:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.79,"Standard Deviation:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Gross C uptake',fontsize=20)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
    # Plot Taylor series below scatter
    PT.plot_taylor_diagram( Taylor_std_list, Taylor_corr_list, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(325), \
                            COLORS=Taylor_col_list, \
                            MARKERS=Taylor_marker_list, \
                            LINEWIDTHS=Taylor_lw_list, \
                            SYMSIZES=Taylor_s_list, \
    )
    
    # Plot Scatter in upper 2/3 of figure and create lists of arrays for taylor plots
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,1), rowspan=2))
    Taylor_std_list=[]
    Taylor_corr_list=[]
    Taylor_col_list=[]
    Taylor_marker_list=[]
    Taylor_lw_list=[]
    Taylor_s_list=[]
    print 'Net C uptake stats:'
    print 5*'%15s' % ('Jsource','PFT','mean bias','std dev','correlation')
    for Jcnt in range(nJSOURCES):
        for iPFT in unique_PFTs: 
            #index = np.array(SITE_PFT)==iPFT
            index = GRASS_TREE_PFT==iPFT
            Jname=J_names[Jcnt]
            print '%15s%15d%15.3f%15.3f%15.3f' % \
                (Jname, iPFT+1, \
                 np.mean(  np.array(JULES_Net_Annual_Cuptake[Jcnt])[index]        \
                         - np.array(SITE_Net_Annual_Cuptake)[index] ),   \
                 np.std(  np.array(JULES_Net_Annual_Cuptake[Jcnt])[index]        \
                        - np.array(SITE_Net_Annual_Cuptake)[index] ),   \
                 stats.pearsonr(np.array(JULES_Net_Annual_Cuptake[Jcnt])[index], \
                                np.array(SITE_Net_Annual_Cuptake)[index] )[0] \
                 )
                 
            Taylor_std_list.append(np.array(NEE_stddev[Jcnt])[index])
            Taylor_corr_list.append(np.array(NEE_correlation[Jcnt])[index])
            Taylor_col_list.append(PS_colours[Jcnt+1])
            Taylor_marker_list.append(LC_markers[iPFT])
            Taylor_lw_list.append(LC_mews[iPFT])
            Taylor_s_list.append(LC_symsizes[iPFT]/2.)
            AX.scatter( np.array(SITE_Net_Annual_Cuptake)[index],  \
                        np.array(JULES_Net_Annual_Cuptake)[Jcnt][index], \
                        c=PS_colours[Jcnt+1], \
                        marker=LC_markers[iPFT], s=LC_symsizes[iPFT],linewidth=LC_mews[iPFT] )
        
        correlation=stats.pearsonr(np.array(SITE_Net_Annual_Cuptake),      \
                                   np.array(JULES_Net_Annual_Cuptake[Jcnt]))[0]
        AX.text(0.1+(Jcnt*0.1),0.92,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        mean_bias = np.mean(  np.array(JULES_Net_Annual_Cuptake)[Jcnt][:]  \
                             -np.array(SITE_Net_Annual_Cuptake)[:]  )
        AX.text(0.1+(Jcnt*0.1),0.84,"%6.2f"%(mean_bias),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        stddev    = np.std(  np.array(SITE_Net_Annual_Cuptake)[:]  \
                           - np.array(JULES_Net_Annual_Cuptake)[Jcnt][:])
        AX.text(0.1+(Jcnt*0.1),0.76,"%6.3f"%(stddev),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )

    AX.text(0.05,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.87,"Mean Bias:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.79,"Standard Deviation:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Net C uptake',fontsize=20)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
    # Plot Taylor series below scatter
    PT.plot_taylor_diagram( Taylor_std_list, Taylor_corr_list, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(326), \
                            COLORS=Taylor_col_list, \
                            MARKERS=Taylor_marker_list, \
                            LINEWIDTHS=Taylor_lw_list, \
                            SYMSIZES=Taylor_s_list, \
                            )
    
    # Create a custom set of handles and labels - symbols for pfts and lines for jules runs
    handles = [LC_handles[iPFT] for iPFT in unique_PFTs] + PS_handles[1:]
    labels  = [LC_names[iPFT] for iPFT in unique_PFTs] + PS_names[1:]
    #handles,labels=AX.get_legend_handles_labels()
    FIG.legend(handles,labels,loc=8) 
    FIG.tight_layout(rect=[0,0,1,0.94])
    FIG.suptitle('Comparison of Annual Carbon Uptake for all PALS sites',\
                 fontsize=24)
    FIG.savefig(plot_dir+'ALLSITES_Flux_Comparison_CarbonUptakeScatter.png', \
                bbox_inches='tight')
    
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()
#######################################################################################    
# plot JULES vs SITE Latent and Sensible
# Scatter and Taylor
#  2x2 plots, columns: GPP, NEE
#             rows: scatter, taylor
if (PLOTS[6]=='Y'):
    # Create nice big figure:
    FIG=plt.figure(figsize=(18,12))
    # Plot Scatter in upper 2/3 of figure and create lists of arrays for taylor plots
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,0), rowspan=2))
    Taylor_std_list=[]
    Taylor_corr_list=[]
    Taylor_col_list=[]
    Taylor_marker_list=[]
    Taylor_lw_list=[]
    Taylor_s_list=[]
    print 'Latent Heat stats:'
    print 5*'%15s' % ('Jsource','PFT','mean bias','std dev','correlation')
    for Jcnt in range(nJSOURCES):
        for iPFT in unique_PFTs: 
            #index = np.array(SITE_PFT)==iPFT
            index = GRASS_TREE_PFT==iPFT
            Jname=J_names[Jcnt]
            #print index
            print '%15s%15d%15.3f%15.3f%15.3f' % \
                (Jname, iPFT+1, \
                 np.mean(  np.array(JULES_Net_Annual_Evap[Jcnt])[index]        \
                         - np.array(SITE_Net_Annual_Evap)[index] ),   \
                 np.std(  np.array(JULES_Net_Annual_Evap[Jcnt])[index]        \
                        - np.array(SITE_Net_Annual_Evap)[index] ),   \
                 stats.pearsonr(np.array(JULES_Net_Annual_Evap[Jcnt])[index], \
                                np.array(SITE_Net_Annual_Evap)[index] )[0] \
                 )
                 
            Taylor_std_list.append(np.array(LHF_stddev[Jcnt])[index])
            Taylor_corr_list.append(np.array(LHF_correlation[Jcnt])[index])
            Taylor_col_list.append(PS_colours[Jcnt+1])
            Taylor_marker_list.append(LC_markers[iPFT])
            Taylor_lw_list.append(LC_mews[iPFT])
            Taylor_s_list.append(LC_symsizes[iPFT]/2.)
            AX.scatter( np.array(SITE_Net_Annual_Evap)[index],  \
                        np.array(JULES_Net_Annual_Evap)[Jcnt][index], \
                        c=PS_colours[Jcnt+1],\
                        marker=LC_markers[iPFT], s=LC_symsizes[iPFT],linewidth=LC_mews[iPFT] )
       
        correlation=stats.pearsonr(np.array(SITE_Net_Annual_Evap),      \
                                   np.array(JULES_Net_Annual_Evap[Jcnt]))[0]
        AX.text(0.1+(Jcnt*0.1),0.92,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        mean_bias = np.mean(  np.array(JULES_Net_Annual_Evap)[Jcnt][:]  \
                             -np.array(SITE_Net_Annual_Evap)[:]  )
        AX.text(0.1+(Jcnt*0.1),0.84,"%6.2f"%(mean_bias),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        stddev = np.std(  np.array(SITE_Net_Annual_Evap)[:]  \
                        - np.array(JULES_Net_Annual_Evap)[Jcnt][:])
        AX.text(0.1+(Jcnt*0.1),0.76,"%6.3f"%(stddev),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        
    AX.text(0.05,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.87,"Mean Bias:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.79,"Standard Deviation:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Net Evaporation',fontsize=20)
    AX.set_xlabel('Site ($kg (water)$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($kg (water)$)')
    AX.set_ylim(limits)
    AX.grid(True)
    # Plot Taylor series below scatter
    PT.plot_taylor_diagram( Taylor_std_list, Taylor_corr_list, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(325), \
                            COLORS=Taylor_col_list, \
                            MARKERS=Taylor_marker_list, \
                            LINEWIDTHS=Taylor_lw_list, \
                            SYMSIZES=Taylor_s_list, \
    )
    
    # Plot Scatter in upper 2/3 of figure and create lists of arrays for taylor plots
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,1), rowspan=2))
    Taylor_std_list=[]
    Taylor_corr_list=[]
    Taylor_col_list=[]
    Taylor_marker_list=[]
    Taylor_lw_list=[]
    Taylor_s_list=[]
    print 'Net Sensible Heat stats:'
    print 5*'%15s' % ('Jsource','PFT','mean bias','std dev','correlation')
    for Jcnt in range(nJSOURCES):
        for iPFT in unique_PFTs: 
            #index = np.array(SITE_PFT)==iPFT
            index = GRASS_TREE_PFT==iPFT
            Jname=J_names[Jcnt]
            print '%15s%15d%15.3f%15.3f%15.3f' % \
                (Jname, iPFT+1, \
                 np.mean(  np.array(JULES_Net_Annual_Heat[Jcnt])[index]        \
                         - np.array(SITE_Net_Annual_Heat)[index] ),   \
                 np.std(  np.array(JULES_Net_Annual_Heat[Jcnt])[index]        \
                        - np.array(SITE_Net_Annual_Heat)[index] ),   \
                 stats.pearsonr(np.array(JULES_Net_Annual_Heat[Jcnt])[index], \
                                np.array(SITE_Net_Annual_Heat)[index] )[0] \
                 )
                 
            Taylor_std_list.append(np.array(SHF_stddev[Jcnt])[index])
            Taylor_corr_list.append(np.array(SHF_correlation[Jcnt])[index])
            Taylor_col_list.append(PS_colours[Jcnt+1])
            Taylor_marker_list.append(LC_markers[iPFT])
            Taylor_lw_list.append(LC_mews[iPFT])
            Taylor_s_list.append(LC_symsizes[iPFT]/2.)
            AX.scatter( np.array(SITE_Net_Annual_Heat)[index],  \
                        np.array(JULES_Net_Annual_Heat)[Jcnt][index], \
                        c=PS_colours[Jcnt+1], \
                        marker=LC_markers[iPFT], s=LC_symsizes[iPFT],linewidth=LC_mews[iPFT] )
        
        correlation=stats.pearsonr(np.array(SITE_Net_Annual_Heat),      \
                                   np.array(JULES_Net_Annual_Heat[Jcnt]))[0]
        AX.text(0.1+(Jcnt*0.1),0.92,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        mean_bias = np.mean(  np.array(JULES_Net_Annual_Heat)[Jcnt][:]  \
                             -np.array(SITE_Net_Annual_Heat)[:]  )
        AX.text(0.1+(Jcnt*0.1),0.84,"%6.2f"%(mean_bias),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )
        stddev    = np.std(  np.array(SITE_Net_Annual_Heat)[:]  \
                           - np.array(JULES_Net_Annual_Heat)[Jcnt][:])
        AX.text(0.1+(Jcnt*0.1),0.76,"%6.3f"%(stddev),       \
                transform=AX.transAxes, fontsize=15, color=PS_colours[Jcnt+1] )

    AX.text(0.05,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.87,"Mean Bias:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    AX.text(0.05,0.79,"Standard Deviation:",transform=AX.transAxes,\
            fontsize=15, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Net Sensible Heat uptake',fontsize=20)
    AX.set_xlabel('Site ($J$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($J$)')
    AX.set_ylim(limits)
    AX.grid(True)
    # Plot Taylor series below scatter
    PT.plot_taylor_diagram( Taylor_std_list, Taylor_corr_list, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(326), \
                            COLORS=Taylor_col_list, \
                            MARKERS=Taylor_marker_list, \
                            LINEWIDTHS=Taylor_lw_list, \
                            SYMSIZES=Taylor_s_list, \
                            )
    
    # Create a custom set of handles and labels - symbols for pfts and lines for jules runs
    handles = [LC_handles[iPFT] for iPFT in unique_PFTs] + PS_handles[1:]
    labels  = [LC_names[iPFT] for iPFT in unique_PFTs] + PS_names[1:]
    #handles,labels=AX.get_legend_handles_labels()
    FIG.legend(handles,labels,loc=8) 
    FIG.tight_layout(rect=[0,0,1,0.94])
    FIG.suptitle('Comparison of Annual Latent and Sensible Heat Exchange for all PALS sites',\
                 fontsize=24)
    print plot_dir+'ALLSITES_Flux_Comparison_HeatUptakeScatter.png'
    FIG.savefig(plot_dir+'ALLSITES_Flux_Comparison_HeatUptakeScatter.png', \
                bbox_inches='tight')
    
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()
   


#####################################################
# plot JULES vs SITE Diurnal GPP and NEE amplitude
# Scatters
if (PLOTS[12]=='Y'):
    # Create figure:
    FIG=plt.figure(figsize=(12,7))
    # Plot Scatter in upper 2/3 of figure
    AX = FIG.add_subplot(1,2,1)
    pars=[]
    for Jcnt in range(nJSOURCES):
        Jname=J_names[Jcnt]
        xdata=np.array(SITE_mean_diurnal_GPP_Amplitude)
        ydata=np.array(JULES_mean_diurnal_GPP_Amplitude[Jcnt])
        AX.scatter( xdata,ydata, \
                    c=PS_colours[Jcnt+1], label=Jname, linewidth=0 )
        pars.append( np.polyfit(xdata,ydata,1) )
        correlation=stats.pearsonr(xdata,ydata)[0]
        AX.text(0.1+(Jcnt*0.1),0.9,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        AX.text(0.1+(Jcnt*0.1),0.8,"%6.2f"%(pars[Jcnt][0]),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )

    AX.text(0.1,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=18, color='k' )
    AX.text(0.1,0.85,"Slope:",transform=AX.transAxes,\
            fontsize=18, color='k' )

    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    for Jcnt in range(nJSOURCES):
        AX.plot(np.array(limits),(pars[Jcnt][0]*np.array(limits))+pars[Jcnt][1],c=PS_colours[Jcnt+1])
    AX.set_title('Diurnal GPP amplitude',fontsize=14)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
    
    AX = FIG.add_subplot(1,2,2)
    pars=[]
    for Jcnt in range(nJSOURCES):
        Jname=J_names[Jcnt]
        xdata=np.array(SITE_mean_diurnal_NEE_Amplitude)
        ydata=np.array(JULES_mean_diurnal_NEE_Amplitude[Jcnt])
        AX.scatter( xdata,ydata, \
                    c=PS_colours[Jcnt+1], label=Jname, linewidth=0 )
        pars.append( np.polyfit(xdata,ydata,1) )
        correlation=stats.pearsonr(xdata,ydata)[0]
        AX.text(0.1+(Jcnt*0.1),0.9,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
        AX.text(0.1+(Jcnt*0.1),0.8,"%6.2f"%(pars[Jcnt][0]),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )
    
    AX.text(0.1,0.95,"Pearson's R:",transform=AX.transAxes,\
            fontsize=18, color='k' )
    AX.text(0.1,0.85,"Slope:",transform=AX.transAxes,\
            fontsize=18, color='k' )

    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    for Jcnt in range(nJSOURCES):
        AX.plot(np.array(limits),(pars[Jcnt][0]*np.array(limits))+pars[Jcnt][1],c=PS_colours[Jcnt+1])
    AX.set_title('NEE',fontsize=17)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
    
    handles,labels=AX.get_legend_handles_labels()
    FIG.legend(handles,labels,loc=8) 
    FIG.tight_layout(rect=[0,0,1,0.94])
    FIG.suptitle('Comparison of Diurnal GPP and NEE amplitude for all PALS sites',\
                 fontsize=24)
    FIG.savefig(plot_dir+'ALLSITES_Flux_Comparison_GPPNEEdiurnalAmplitude.png', \
                bbox_inches='tight')
    
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()
   
# plot JULES vs SITE Seasonal Diurnal GPP/NEE amplitude
# Scatter
#  2x4 plots, columns: Seasons (DJF,MAM,JJA,SON)
#             rows: GPP, NEE
if (PLOTS[13]=='Y'):
    # Create figure and axes:
    FIG, AXES = plt.subplots(nrows=2,ncols=4,figsize=[18,9]) 
    # Plot Scatter in upper 2/3 of figure
    SITE_Amp_list = [ SITE_mean_diurnal_GPP_Amplitude_BySeas, \
                      SITE_mean_diurnal_NEE_Amplitude_BySeas  ]
    JULES_Amp_list = [ JULES_mean_diurnal_GPP_Amplitude_BySeas, \
                       JULES_mean_diurnal_NEE_Amplitude_BySeas  ]
    var_names = ['GPP','NEE']
    nVARS = len(var_names)
    max_limits = [0 for i in range(nVARS)]
    units='$\mu$mol $m^{-2} s^{-2}$'
    pars = [ [ [ [] for i in range(nJSOURCES) ] for j in range(nSEAS) ] for k in range(nVARS) ]
    for iVAR in range(nVARS):
        for iSEAS in range(nSEAS):
            for Jcnt in range(nJSOURCES):
                Jname=J_names[Jcnt]
                xdata=np.array(SITE_Amp_list[iVAR][iSEAS])
                jdata=np.array(JULES_Amp_list[iVAR][iSEAS][Jcnt])
                AXES[iVAR,iSEAS].scatter( xdata,ydata, \
                                          c=PS_colours[Jcnt+1], label=Jname, linewidth=0 )
                par = np.polyfit(xdata,ydata,1) 
                pars[iVAR][iSEAS][Jcnt]=  par 
                correlation=stats.pearsonr(xdata,ydata)[0]
                AXES[iVAR,iSEAS].text(0.1+(Jcnt*0.1),0.9,"%6.2f"%(correlation),       \
                        transform=AXES[iVAR,iSEAS].transAxes, fontsize=12, color=PS_colours[Jcnt+1] )
                AXES[iVAR,iSEAS].text(0.1+(Jcnt*0.1),0.8,"%6.2f"%(par[0]),       \
                        transform=AXES[iVAR,iSEAS].transAxes, fontsize=12, color=PS_colours[Jcnt+1] )
                
            AXES[iVAR,iSEAS].text(0.1,0.95,"Pearson's R:",transform=AXES[iVAR,iSEAS].transAxes,\
                    fontsize=12, color='k' )
            AXES[iVAR,iSEAS].text(0.1,0.85,"Slope:",transform=AXES[iVAR,iSEAS].transAxes,\
                    fontsize=12, color='k' )
            limits=[min(AXES[iVAR,iSEAS].get_xlim()[0],AXES[iVAR,iSEAS].get_ylim()[0],0),\
                    max(AXES[iVAR,iSEAS].get_xlim()[1],AXES[iVAR,iSEAS].get_ylim()[1]) ]
             
            if iVAR==0:
                AXES[iVAR,iSEAS].set_title(SEAS_names[iSEAS],fontsize=17)
            
            if iVAR==nVARS-1:
                AXES[iVAR,iSEAS].set_xlabel('Site '+var_names[iVAR]+' '+units,fontsize=14)
                    
            if iSEAS==0:
                AXES[iVAR,iSEAS].set_ylabel('JULES '+var_names[iVAR]+' '+units,fontsize=14)

            AXES[iVAR,iSEAS].grid(True)
                    
            # Store max ylimits for both GPP and NEE seperately so all seasons plotted
            # on the same scale
            max_limits[iVAR] = max([max_limits[iVAR]]+limits ) 
    
    for iSEAS in range(nSEAS):
        limits=[0,max_limits[0]]
        AXES[0,iSEAS].set_ylim( limits )
        AXES[0,iSEAS].set_xlim( limits )
        AXES[0,iSEAS].plot(np.array(limits),np.array(limits),c='k')
        for Jcnt in range(nJSOURCES):
            AXES[0,iSEAS].plot(np.array(limits),\
                                ( pars[0][iSEAS][Jcnt][0]*np.array(limits)) + \
                                  pars[0][iSEAS][Jcnt][1], \
                                c=PS_colours[Jcnt+1], linewidth=0)

        limits=[-max_limits[1],max_limits[1]]
        AXES[1,iSEAS].set_ylim( limits )
        AXES[1,iSEAS].set_xlim( limits )
        AXES[1,iSEAS].plot(np.array(limits),np.array(limits),c='k')
        for Jcnt in range(nJSOURCES):
            AXES[1,iSEAS].plot(np.array(limits),\
                                ( pars[1][iSEAS][Jcnt][0]*np.array(limits)) + \
                                  pars[1][iSEAS][Jcnt][1],c=PS_colours[Jcnt+1])

    handles,labels=AX.get_legend_handles_labels()
    FIG.legend(handles,labels,loc=8) 
    FIG.tight_layout(rect=[0,0,1,0.94])
    FIG.suptitle('Comparison of Diurnal GPP and NEE amplitude for all PALS sites',\
                 fontsize=24)
    FIG.savefig(plot_dir+'ALLSITES_Flux_Comparison_GPPNEEdiurnalAmplitude.png', \
                bbox_inches='tight')
    
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


