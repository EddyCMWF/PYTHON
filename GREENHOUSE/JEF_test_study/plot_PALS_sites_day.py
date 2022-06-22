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
site_metadata_file = PALS_dir+'sites_meta.csv'
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
print 'nJSOURCES, J_SOURCES = ',nJSOURCES,J_SOURCES


    
PLOTS=''
if (INTERACTIVE=='Y'):
    PLOTS+=raw_input('Plot time-series? (Y/N) ')
    PLOTS+=raw_input('Plot LHF and SHF scatter? (Y/N) ')
    PLOTS+=raw_input('Plot NEE and GPP scatter? (Y/N) ')
    PLOTS+=raw_input('Plot LHF vs GPP scatter? (Y/N) ')
    PLOTS+='N' #raw_input('(Y/N) ')
    PLOTS+=raw_input('Plot scatter of C uptake for each site? (Y/N) ')
    PLOTS+='N' #raw_input('Plot LHF/NEE scatter? (Y/N) ')
else:
    PLOTS=sys.argv[3]

##########################################################
plot_dir=BASE_plot_dir
if (INTERACTIVE=='Y'):
    print 'Plot output directory: '+plot_dir
    plot_dir_temp=raw_input('Enter alternative plot output directory or hit return to continue: ')
    if (len(plot_dir_temp)>0):
        plot_dir=plot_dir_temp
# Create plot_dir if doesn't already exist
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

SITE_lats = []
SITE_lons = []
SITE_PFT  = []

#create lists for storing data to plot for all sites
SITE_Net_Annual_Cuptake = []
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
    
    print 'Site: '+site_name
    
    ######################################
    # Read in-situ site data
    SITE_fname = site_name+'Fluxnet.1.4_flux'+Tres_Stag+'.nc'
    Sinf       = nc.Dataset(SITE_data_dir+SITE_fname,'r')
    ############################################################
    # Extract data and convert to same units (Site data units)
    S_data = {'time':nctime.num2date(Sinf.variables['time'][:],                 \
                                     units=Sinf.variables['time'].units),       \
              'NEE':Sinf.variables['NEE'][:].squeeze(),                         \
              'GPP':np.ma.masked_less(Sinf.variables['GPP'][:].squeeze(),1e-6), \
              'LHF':Sinf.variables['Qle'][:].squeeze(),                         \
              'SHF':Sinf.variables['Qh'][:].squeeze(),                          \
              }

    # Store Annual C uptake in a list for later plots
    SITE_Net_Annual_Cuptake.append( (np.mean(S_data['NEE'])/kgC_to_umolsCO2_factor)\
                                    * seconds_in_year*-1. )
    S_data['FQW'] = S_data['LHF']/Lc_H2O
    S_data['WUE'] = S_data['GPP']/S_data['FQW']
    
    Sinf.close()

    ############################################
    # Loop round JULES sources and extract data to dictionary of lists
    J_data = { 'time':[],\
               'NPP':[],\
               'GPP':[],\
               'resp_s':[],\
               'FQW':[],\
               'LHF':[],\
               'SHF':[],\
               'NEE':[],\
               'WUE':[],\
               'GPP_index':[],\
               }

    for Jcnt in range(nJSOURCES):
        iSOURCE = J_SOURCES[Jcnt]
        print iSOURCE, Jcnt, len(JULES_sources), nJSOURCES
        Jsource = JULES_sources[iSOURCE]
        ######################################
        # Construct filename and open:
        JULES_fname=site_name+'.'+Tres+'.nc'
        Jinf = nc.Dataset(Jsource[1]['data_dir']+JULES_fname,'r')
        #temp_time = Jinf.variables['time_bounds'][:,0]
        #J_data['time'].append(nctime.num2date(temp_time, \
        J_data['time'].append(nctime.num2date(Jinf.variables['time_bounds'][:,0], \
                                              units=Jinf.variables['time'].units) )
        J_data['NPP'].append( Jinf.variables[Jsource[1]['npp_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor   )
        #J_data['GPP'].append( Jinf.variables[Jsource[1]['gpp_name']][:].squeeze() \
        #                          *   kgC_to_umolsCO2_factor )
        # MASK out GPP data where GPP is less than 0
        J_data['GPP'].append( np.ma.masked_less( Jinf.variables[Jsource[1]['gpp_name']][:].squeeze() \
                                  *kgC_to_umolsCO2_factor, 1e6 )
        
        J_data['resp_s'].append( Jinf.variables[Jsource[1]['resp_s_name']][:].squeeze() \
                                  *   kgC_to_umolsCO2_factor  )
        J_data['FQW'].append( Jinf.variables['fqw_gb'][:].squeeze() )
        J_data['LHF'].append( Jinf.variables['latent_heat'][:].squeeze() )
        J_data['SHF'].append( Jinf.variables['ftl_gb'][:].squeeze()  )
        Jinf.close()
        
        print 'Site: ',S_data['time'][0], S_data['time'][-1]
        print 'JULES: ',J_data['time'][Jcnt][0], J_data['time'][Jcnt][-1]
        
        # Now store the secondary product
        J_data['NEE'].append( (J_data['NPP'][Jcnt]-J_data['resp_s'][Jcnt])*-1 )
        
        J_data['WUE'].append( J_data['GPP'][Jcnt]/J_data['FQW'][Jcnt] )
        
        #J_data['GPP_index'].append( (J_data['GPP'][Jcnt]>1e-6) & (S_data['GPP']>1e6) )

        JULES_Net_Annual_Cuptake[Jcnt].append( np.mean(J_data['NEE'][Jcnt]) \
                                                * seconds_in_year*-1./kgC_to_umolsCO2_factor )
        JULES_Gross_Annual_Cuptake[Jcnt].append( np.mean(J_data['GPP'][Jcnt]) \
                                                 * seconds_in_year/kgC_to_umolsCO2_factor )
        
        # Calculate and store the Specific Heat Flux Statistics  
        SHF_meanbias[Jcnt].append( TST.meanbias( J_data['SHF'][Jcnt],S_data['SHF'],  \
                                                J_data['time'][Jcnt],S_data['time'] )  )
        SHF_stddev[Jcnt].append( TST.stddev( J_data['SHF'][Jcnt],S_data['SHF'],  \
                                            J_data['time'][Jcnt],S_data['time'] )  )
        SHF_correlation[Jcnt].append( TST.pearsonr( J_data['SHF'][Jcnt],S_data['SHF'],  \
                                                   J_data['time'][Jcnt],S_data['time'] )[0]  )
        
        # Calculate and store the LAtent Heat Flux Stats 
        LHF_meanbias[Jcnt].append( TST.meanbias( J_data['LHF'][Jcnt],S_data['LHF'],  \
                                                J_data['time'][Jcnt],S_data['time'] )  )
        LHF_stddev[Jcnt].append( TST.stddev( J_data['LHF'][Jcnt],S_data['LHF'],  \
                                            J_data['time'][Jcnt],S_data['time'] )  )
        LHF_correlation[Jcnt].append( TST.pearsonr( J_data['LHF'][Jcnt],S_data['LHF'],  \
                                                   J_data['time'][Jcnt],S_data['time'] )[0]  )
    
        # Calculate and store the NEE Stats 
        NEE_meanbias[Jcnt].append( TST.meanbias( J_data['NEE'][Jcnt],S_data['NEE'],  \
                                                J_data['time'][Jcnt],S_data['time'] )  )
        NEE_stddev[Jcnt].append( TST.stddev( J_data['NEE'][Jcnt],S_data['NEE'],  \
                                            J_data['time'][Jcnt],S_data['time'] )  )
        NEE_correlation[Jcnt].append( TST.pearsonr( J_data['NEE'][Jcnt],S_data['NEE'],  \
                                                   J_data['time'][Jcnt],S_data['time'] )[0] )

        # Calculate and store the GPP Stats 
        GPP_meanbias[Jcnt].append( TST.meanbias( J_data['GPP'][Jcnt],S_data['GPP'],  \
                                                 J_data['time'][Jcnt],S_data['time'] )  )
        GPP_stddev[Jcnt].append( TST.stddev( J_data['GPP'][Jcnt],S_data['GPP'],  \
                                             J_data['time'][Jcnt],S_data['time'] )  )
        GPP_correlation[Jcnt].append( TST.pearsonr( J_data['GPP'][Jcnt],S_data['GPP'],  \
                                                    J_data['time'][Jcnt],S_data['time'] )[0] )

    ##########################################################
    # Plot time-series of SHF, LHF, and NEE
    if (PLOTS[0]=='Y'):
        FIG = plt.figure(figsize=(13,12))
    
        AX = FIG.add_subplot(3,1,1)
        AX.plot(S_data['time'],S_data['SHF'],lw=2,c=SITE_colour)
        AX.text(1.03,0.2,"Correlations with site: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.plot( J_data['time'][Jcnt],J_data['SHF'][Jcnt],lw=1.5,\
                      c=Jsource[1]['colour'] )
            AX.text(1.13+(Jcnt*0.1),0.1,"%6.2f"%(SHF_correlation[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour'] )
        AX.set_ylabel('Sensible Heat Flux (W $m^{-2}$)')
        
        AX = FIG.add_subplot(3,1,2)
        AX.plot(S_data['time'],S_data['LHF'],lw=2,c=SITE_colour)
        AX.text(1.03,0.2,"Correlations with site: ",transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.plot( J_data['time'][Jcnt],J_data['LHF'][Jcnt],lw=1.5,\
                      c=Jsource[1]['colour'] )
            AX.text(1.13+(Jcnt*0.1),0.1,"%6.2f"%(LHF_correlation[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour'] )
        AX.set_ylabel('Latent Heat Flux (W $m^{-2}$)')
        
        AX = FIG.add_subplot(3,1,3)
        AX.plot(S_data['time'],S_data['NEE'],lw=2,c=SITE_colour,label='Site')
        AX.plot(S_data['time'],S_data['GPP'],lw=2,c=SITE_colour,ls=':')
        AX.text(1.03,0.5,"Net C Uptake \n (kgC y$^{-1}$ m$^{-2}$):",\
                transform=AX.transAxes, fontsize=14)
        AX.text(1.03,0.42,"%6.3f"%(SITE_Net_Annual_Cuptake[cnt]),\
                transform=AX.transAxes, fontsize=14, color=SITE_colour )
        AX.text(1.03,0.1,"Correlations with site: ",transform=AX.transAxes, fontsize=14)
        AX.text(1.03,0.85,"Gross C Uptake \n(gC y$^{-1}$ m$^{-2}$): ",\
                transform=AX.transAxes, fontsize=14)
        for Jcnt in range(nJSOURCES):
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.plot( J_data['time'][Jcnt],J_data['NEE'][Jcnt],lw=1.5,\
                      c=Jsource[1]['colour'], label=Jsource[1]['longname'] )
            AX.plot( J_data['time'][Jcnt],J_data['GPP'][Jcnt],lw=1.5,\
                    c=Jsource[1]['colour'], ls=':' )
            AX.text(1.13+(Jcnt*0.1),0.0,"%6.2f"%(NEE_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour'] )
            AX.text(1.13+(Jcnt*0.1),0.42,"%6.3f"%(JULES_Net_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour']    )
            AX.text(1.13+(Jcnt*0.1),0.77,"%6.3f"%(JULES_Gross_Annual_Cuptake[Jcnt][cnt]), \
                    transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour']      )
        AX.set_ylabel('NEE and GPP ($\mu$mol $m^{-2} s^{-1}$)')
        
        AX.legend( bbox_to_anchor=(0.5,-0.25),loc=8,borderaxespad=0.,ncol=nJSOURCES+1)
        
        FIG.tight_layout(rect=[0,0.05,0.8,0.96])
        FIG.suptitle(site_name+'flux tower time-series, obs and J-E-F',\
                     fontsize=24)
        FIG.text(0.8,0.92,'Country: '+site_country,fontsize=16)
        FIG.text(0.8,0.89,'PFT: '+site_vegtype,fontsize=16)
        FIG.text(0.8,0.86,'Lat: '+site_lat+'$^o$E',fontsize=16)
        FIG.text(0.8,0.83,'Lon: '+site_lon+'$^o$N',fontsize=16)
        FIG.text(0.8,0.80,'Can. hgt: '+site_canht+' m',fontsize=16)
        FIG.text(0.8,0.77,'Meas. hgt: '+site_measht+' m',fontsize=16)
        
        
        FIG.savefig(plot_dir+site_name+'_Flux_Comparison_Time-Series.png', \
                    bbox_inches='tight')
    
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
            
             
    ##########################################################
    # Plot scatter plot of LHF and SHF data
        
    if (PLOTS[1]=='Y'):
        FIG = plt.figure(figsize=(18,9))
        
        # Plot the LHF data
        AX = FIG.add_subplot(1,2,1)
        for Jcnt in range(nJSOURCES):
            overlapJ,overlapS=TST.overlap_indexes(J_data['time'][Jcnt],S_data['time'])
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.scatter( S_data['LHF'][overlapS],J_data['LHF'][Jcnt][overlapJ],\
                      c=Jsource[1]['colour'], label=Jsource[1]['longname'] )
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(LHF_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=Jsource[1]['colour'] )
        
        AX.set_xlabel('Site Latent Heat Flux (W $m^{-2}$)')
        AX.set_ylabel('JULES Latent Heat Flux (W $m^{-2}$)')
        AX.legend( bbox_to_anchor=(1.0,-0.25),loc=8,borderaxespad=0.,ncol=nJSOURCES+1)

        # Plot the SHF data
        AX = FIG.add_subplot(1,2,2)
        for Jcnt in range(nJSOURCES):
            overlapJ,overlapS=TST.overlap_indexes(J_data['time'][Jcnt],S_data['time'])
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.scatter( S_data['SHF'][overlapS],J_data['SHF'][Jcnt][overlapJ],\
                      c=Jsource[1]['colour'])
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(SHF_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=Jsource[1]['colour'] )
        
        AX.set_xlabel('Site Sensible Heat Flux (W $m^{-2}$)')
        AX.set_ylabel('JULES Sensible Heat Flux (W $m^{-2}$)')

        FIG.tight_layout(rect=[0,0,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, Heat Flux Comparison',\
                     fontsize=24)
        FIG.savefig(plot_dir+site_name+'_Flux_Comparison_HFscatter.png', \
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
            overlapJ,overlapS=TST.overlap_indexes(J_data['time'][Jcnt],S_data['time'])
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.scatter( S_data['NEE'][overlapS],J_data['NEE'][Jcnt][overlapJ],\
                      c=Jsource[1]['colour'], label=Jsource[1]['longname'] )
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(NEE_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=Jsource[1]['colour'] )
        
        AX.set_xlabel('Site Net Ecosystem Exchange ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylabel('JULES Net Ecosystem Exchange ($\mu$mol $m^{-2} s^{-2}$)')
        AX.legend( bbox_to_anchor=(1.0,-0.25),loc=8,borderaxespad=0.,ncol=nJSOURCES+1)
        

        # Plot the GPP data
        AX = FIG.add_subplot(1,2,2)
        for Jcnt in range(nJSOURCES):
            overlapJ,overlapS=TST.overlap_indexes(J_data['time'][Jcnt],S_data['time'])
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.scatter( S_data['GPP'][overlapS],J_data['GPP'][Jcnt][overlapJ],\
                      c=Jsource[1]['colour'])
            AX.text(0.6+(Jcnt*0.1),0.05,"%6.2f"%(GPP_correlation[Jcnt][cnt]),       \
                    transform=AX.transAxes, fontsize=18, color=Jsource[1]['colour'] )
        
        AX.set_xlabel('Site Gross Primary Productivity ($\mu$mol $m^{-2} s^{-2}$)')
        AX.set_ylabel('JULES Gross Primary Productivity ($\mu$mol $m^{-2} s^{-2}$)')

        FIG.tight_layout(rect=[0,0,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, Heat Flux Comparison',\
                     fontsize=24)
        FIG.savefig(plot_dir+site_name+'_Flux_Comparison_NEEGPPscatter.png', \
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
        
        AX.scatter(S_data['LHF'],S_data['GPP'],c=SITE_colour,label='Site',marker='x')

        for Jcnt in range(nJSOURCES):
            Jsource = JULES_sources[J_SOURCES[Jcnt]]
            AX.scatter( J_data['LHF'][Jcnt],J_data['GPP'][Jcnt],\
                      c=Jsource[1]['colour'], label=Jsource[1]['longname'] )
            #AX.text(1.13+(Jcnt*0.1),0.0,"%6.2f"%(NEE_correlation[Jcnt][cnt]),       \
            #        transform=AX.transAxes, fontsize=14, color=Jsource[1]['colour'] )
        
        AX.set_xlabel('Latent Heat Flux (W $m^{-2}$)')
        AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)')
        AX.legend()
        
        FIG.tight_layout(rect=[0,0,1,0.94])
        FIG.suptitle(site_name+'('+site_country+') flux tower, LHF vs GPP',\
                     fontsize=24)
        FIG.savefig(plot_dir+site_name+'_Flux_Comparison_LHFvGPPscatter.png', \
                    bbox_inches='tight')
        
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()
    

# plot JULES vs SITE Cuptake (GROSS and NET) 
#if PLOTS[5]:
    

