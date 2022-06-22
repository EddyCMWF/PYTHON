#!/bin/python
#
# Scripts to plot output from the JULES simulations alongside the field data.
# 
#  ./plot_PALS_sites_day.py [INTERACTIVE] [iDISPLAY] [PLOTS] [SITES] [Jsources]
#  ./plot_PALS_sites_day.py Y
#  ./plot_PALS_sites_day.py N N YNNN [0,1,2,3] [0,1]
# 
# 

import netCDF4 as nc
import netcdftime as nctime
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
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
SITE_colour = 'blue'
JULES_colour= 'green'

# T resolution hard coded to daily
Tres = 'day'
Tres_Stag='_daily'

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
        SITES = range(sSITE,eSITE)
        del sSITE
        del eSITE
    else:
        SITES = [int(site) for site in SITES_temp.split(',')]
    
    del SITES_temp
else:
    SITES_temp = sys.argv[4]
    if SITES_temp=='ALL':
        SITES=range(len(meta_data))
    else:
        SITES_temp=filter(None,SITES_temp.replace('[','').replace(']','').split(','))
        SITES=[int(temp) for temp in SITES_temp]

nSITES=len(SITES)
print 'nSITES, SITES = ',nSITES,SITES

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
    
PLOTS=''
if (INTERACTIVE=='Y'):
    PLOTS+=raw_input('Plot time-series? (Y/N) ')  # PLOTS[0]
    PLOTS+=raw_input('Plot LHF and SHF scatter? (Y/N) ')  # PLOTS[1] 
    PLOTS+=raw_input('Plot NEE and GPP scatter? (Y/N) ')  # PLOTS[2]
    PLOTS+=raw_input('Plot LHF vs GPP scatter? (Y/N) ')  # PLOTS[3]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[4]
    PLOTS+=raw_input('Plot scatter and Taylor diagrams of C uptake for all sites? (Y/N) ')  # PLOTS[5]
    PLOTS+='N' # raw_input('(Y/N) ')  # PLOTS[6] 
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[7]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[8]
    PLOTS+='N' #raw_input('(Y/N) ')  # PLOTS[9]
    PLOTS+='N'#raw_input('Plot mean and standard deviation diurnal behviour of GPP and NEE (Y/N) ')  # PLOTS[10]
    PLOTS+='N'#raw_input('Plot seasonal mean and standard deviation diurnal behviour of GPP and NEE (Y/N) ')  # PLOTS[11]
else:
    PLOTS=sys.argv[3]

##########################################################
plot_dir=BASE_plot_dir
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

##########################################################
# If INTERACTIVE, print resubmit command for ease next time:
if (INTERACTIVE=='Y'):
    print 'Resubmit Command: '
    print './plot_PALS_sites_day_pandas.py N '+\
            iDISPLAY+' '+PLOTS+' '+\
            str(SITES).replace(' ','')+' '+\
            str(J_SOURCES).replace(' ','')
    if (plot_dir!=BASE_plot_dir):
        print '-plot_dir '+plot_dir

SITE_lats = []
SITE_lons = []
SITE_PFT  = []

#create lists for storing data to plot for all sites
SITE_Net_Annual_Cuptake = []
SITE_Gross_Annual_Cuptake = []
JULES_Net_Annual_Cuptake = [ [] for i in range(nJSOURCES) ]
JULES_Gross_Annual_Cuptake = [ [] for i in range(nJSOURCES) ]

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

for cnt in range(nSITES):
    iSITE=SITES[cnt]
    site_meta=meta_data[iSITE]
    frac=np.array(frac_list[iSITE].split(),dtype='float')
    site_name=site_meta[0].strip()
    site_country=site_meta[1].strip()
    site_vegtype=site_meta[2].strip()
    site_canht=site_meta[3].strip()
    site_measht=site_meta[4].strip()
    site_lat=site_meta[7].strip()
    site_lon=site_meta[8].strip()

    SITE_lats.append(float(site_lat))
    SITE_lons.append(float(site_lon))
    SITE_PFT.append(site_vegtype)
    
    site_plot_dir = plot_dir+site_name+'/'
    os.system('mkdir -p '+site_plot_dir)

    print 'Site: '+site_name
    print 'Plots saved in: '+site_plot_dir
    
    ######################################
    # Read in-situ site data
    SITE_fname = site_name+'Fluxnet.1.4_flux'+Tres_Stag+'.nc'
    Sinf       = nc.Dataset(SITE_data_dir+SITE_fname,'r')
    ############################################################
    # Extract data and convert to same units (Site data units)
    #'time':nctime.num2date(Sinf.variables['time'][:],units=Sinf.variables['time'].units),       \
    S_data = {'NEE':Sinf.variables['NEE'][:].squeeze(),                         \
              'GPP':np.ma.masked_less(Sinf.variables['GPP'][:].squeeze(),0.01), \
              'LHF':Sinf.variables['Qle'][:].squeeze(),                         \
              'SHF':Sinf.variables['Qh'][:].squeeze(),                          \
              }
    S_time = nctime.num2date(Sinf.variables['time'][:], \
                             units=Sinf.variables['time'].units)
    S_data['FQW'] = S_data['LHF']/Lc_H2O
    S_data['WUE'] = S_data['GPP']/S_data['FQW']
    
    S_panda = pd.DataFrame( S_data, index=S_time, )
    
    # Store Annual C uptake in a list for later plots
    SITE_Net_Annual_Cuptake.append( (S_panda['NEE'].mean()/kgC_to_umolsCO2_factor)\
                                    * seconds_in_year*-1. )
    SITE_Gross_Annual_Cuptake.append( (S_panda['GPP'].mean()/kgC_to_umolsCO2_factor)\
                                    * seconds_in_year )
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
        JULES_fname=site_name+'.'+Tres+'.nc'
        Jinf = nc.Dataset(Jsource[1]['data_dir']+JULES_fname,'r')
        J_data = { 'NPP':Jinf.variables[Jsource[1]['npp_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor,                  \
                   'GPP':np.ma.masked_less(Jinf.variables[Jsource[1]['gpp_name']][:].squeeze() \
                                  *kgC_to_umolsCO2_factor, 0.01 ),               \
                   'resp_s':Jinf.variables[Jsource[1]['resp_s_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor,                   \
                   'FQW':Jinf.variables['fqw_gb'][:].squeeze(),                 \
                   'LHF':Jinf.variables['latent_heat'][:].squeeze(),            \
                   'SHF':Jinf.variables['ftl_gb'][:].squeeze(),                 \
                   }

        J_data['NEE']=(J_data['NPP']-J_data['resp_s'])*-1 
        J_data['WUE']= J_data['GPP']/J_data['FQW']
        
        J_time=nctime.num2date(Jinf.variables['time_bounds'][:,0], \
                               units=Jinf.variables['time'].units  )
        
        J_panda = pd.DataFrame(J_data,index=J_time)
        J_pandas.append(J_panda.copy())
        
        Jinf.close()
        
        JULES_Net_Annual_Cuptake[Jcnt].append(J_panda['NEE'].mean()*   \
                                              seconds_in_year*-1./kgC_to_umolsCO2_factor )
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
        GPP_meanbias[Jcnt].append((J_panda['GPP']-S_panda['GPP']).mean()) 
        GPP_stddev[Jcnt].append( (J_panda['GPP']-S_panda['GPP']).std() )
        GPP_correlation[Jcnt].append( J_panda['GPP'].corr(S_panda['GPP']))
        
        
    #Create Parameter Pandas
    LHF_panda=pd.concat([S_panda['LHF']]+[J_panda['LHF'] for J_panda in J_pandas],axis=1)
    LHF_panda.columns=(PS_names)
    SHF_panda=pd.concat([S_panda['SHF']]+[J_panda['SHF'] for J_panda in J_pandas],axis=1)
    SHF_panda.columns=(PS_names)
    NEE_panda=pd.concat([S_panda['NEE']]+[J_panda['NEE'] for J_panda in J_pandas],axis=1)
    NEE_panda.columns=(PS_names)
    GPP_panda=pd.concat([S_panda['GPP']]+[J_panda['GPP'] for J_panda in J_pandas],axis=1)
    GPP_panda.columns=(PS_names)
        
    ##########################################################
    # Plot time-series of SHF, LHF, and NEE
    if (PLOTS[0]=='Y'):

        FIG = plt.figure(figsize=(13,12))
        
        AX = FIG.add_subplot(3,1,1)
        SHF_panda.plot(lw=1.5,ax=AX,color=PS_colours,legend=False)
        AX.text(1.02,0.22,"Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.08,"SHF Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            AX.text(1.02+(Jcnt*0.1),0.16,"%6.2f"%(SHF_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.1),0.02,"%6.2f"%(SHF_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
        AX.set_ylabel('Sensible Heat Flux (W $m^{-2}$)', fontsize=14)
        
        AX = FIG.add_subplot(3,1,2)
        LHF_panda.plot(lw=1.5,ax=AX,color=PS_colours,legend=False)
        AX.text(1.02,0.22,"LHF Mean Bias: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.02,0.08,"LHF Std. Dev.: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            AX.text(1.02+(Jcnt*0.1),0.16,"%6.2f"%(LHF_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.1),0.02,"%6.2f"%(LHF_stddev[Jcnt][cnt]), \
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
            AX.text(1.02+(Jcnt*0.1),0.38,"%6.2f"%(NEE_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.1),0.26,"%6.2f"%(NEE_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.1),0.10,"%6.2f"%(GPP_meanbias[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.02+(Jcnt*0.1),-0.04,"%6.2f"%(GPP_stddev[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.13+(Jcnt*0.1),0.71,"%6.3f"%(JULES_Net_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )
            AX.text(1.13+(Jcnt*0.1),0.87,"%6.3f"%(JULES_Gross_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=PS_colours[Jcnt+1] )

        AX.set_ylabel('NEE and GPP ($\mu$mol $m^{-2} s^{-1}$)')

        AX.text(1.03,0.95,"Gross C Uptake: ",\
                transform=AX.transAxes, fontsize=14)
        AX.text(1.03,0.87,"%6.3f"%(SITE_Gross_Annual_Cuptake[cnt]),\
                transform=AX.transAxes, fontsize=14, color=SITE_colour )
        AX.text(1.03,0.79,"Net C Uptake:",\
                transform=AX.transAxes, fontsize=14)
        AX.text(1.03,0.71,"%6.3f"%(SITE_Net_Annual_Cuptake[cnt]),\
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
                        c=PS_colours[Jcnt+1], label=Jname )
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
                        c=PS_colours[Jcnt+1], label=Jname )
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
                        c=PS_colours[Jcnt+1], label=Jname )
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
                        c=PS_colours[Jcnt+1], label=Jname )
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
        
        for Jcnt in range(nJSOURCES):
            AX.scatter( J_pandas[Jcnt]['LHF'],J_pandas[Jcnt]['GPP'],\
                      c=PS_colours[Jcnt+1], label=J_names[Jcnt] )
        
        AX.scatter(S_panda['LHF'],S_panda['GPP'],c=PS_colours[0],label='Site',marker='x')

        AX.set_xlabel('Latent Heat Flux (W $m^{-2}$)')
        AX.set_xlim([0,AX.get_xlim()[1]])
        AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylim([0,AX.get_ylim()[1]])
        AX.legend()
        
        FIG.tight_layout(rect=[0,0,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, LHF vs GPP',\
                     fontsize=24)
        FIG.savefig(site_plot_dir+site_name+'_Flux_Comparison_LHFvGPPscatter.png', \
                    bbox_inches='tight')
        
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()


# plot JULES vs SITE Cuptake (GROSS and NET)
# Scatter and Taylor
#  2x2 plots, columns: GPP, NEE
#             rows: scatter, taylor
if PLOTS[5]:
    # Create nice big figure:
    FIG=plt.figure(figsize=(18,12))

    # Plot Scatter in upper 2/3 of figure
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,0), rowspan=2))
    for Jcnt in range(nJSOURCES):
        Jname=J_names[Jcnt]
        AX.scatter( np.array(SITE_Gross_Annual_Cuptake),  \
                    np.array(JULES_Gross_Annual_Cuptake[Jcnt]), \
                    c=PS_colours[Jcnt+1], label=Jname )
        correlation=stats.pearsonr(np.array(SITE_Gross_Annual_Cuptake),      \
                                   np.array(JULES_Gross_Annual_Cuptake[Jcnt]))[0]
        AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )

    AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
            fontsize=18, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Gross C uptake',fontsize=14)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
   
    # Plot Taylor series below scatter
    PT.plot_taylor_diagram( GPP_stddev, GPP_correlation, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(325), \
                            COLORS=PS_colours[1:], \
                            )

    
    AX = FIG.add_subplot(plt.subplot2grid( (3,2), (0,1), rowspan=2))
    for Jcnt in range(nJSOURCES):
        Jname=J_names[Jcnt]
        AX.scatter( np.array(SITE_Net_Annual_Cuptake),  \
                    np.array(JULES_Net_Annual_Cuptake[Jcnt]), \
                    c=PS_colours[Jcnt+1], label=Jname )
        correlation=stats.pearsonr(np.array(SITE_Net_Annual_Cuptake),      \
                                   np.array(JULES_Net_Annual_Cuptake[Jcnt]))[0]
        AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(correlation),       \
                transform=AX.transAxes, fontsize=18, color=PS_colours[Jcnt+1] )

    AX.text(0.6,0.1,"Pearson's R:",transform=AX.transAxes,\
            fontsize=18, color='k' )
    limits=[min(AX.get_xlim()[0],AX.get_ylim()[0],0),\
            max(AX.get_xlim()[1],AX.get_ylim()[1]) ]
    AX.plot(np.array(limits),np.array(limits),c='k')
    AX.set_title('Net C uptake',fontsize=14)
    AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_xlim(limits)
    AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)')
    AX.set_ylim(limits)
    AX.grid(True)
    handles,labels=AX.get_legend_handles_labels()

    PT.plot_taylor_diagram( NEE_stddev, NEE_correlation, \
                            corr_range=[-1,1] , \
                            FIG=FIG,FIG_subplot=(326), \
                            COLORS=PS_colours[1:], \
                            )

    
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
   


    
