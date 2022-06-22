#!/usr/bin/python
#
# Python to plot uk site data against jules output
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# April 2015
#
# Contains
import os, sys
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import csv
import plot_tools as PT
#import pylab as plt
import matplotlib.pyplot as plt
import UKsite_info
import maths_tools as MT
#
JULES_output_DIR='/users/eow/edwcom/UKsitecomp/JULES_output/'
#
UKsites_DIR='/data/grp/fluxdata/uk_sites/data/'
#
OUTPUT_DIR='/users/eow/edwcom/UKsitecomp/plots/'
#
site_names  = [ 'AliceHolt', 'Cardington','EasterBush' ]  #  'Auchencorth','EasterBush', , 'GriffinForrest'

# kgC to umol conversion factor
kgC2umol=(1e9/12.)

for site in ['EastSaltoun']:    # site_names:     ##['EasterBush']:    # site_dict:
    #
    print site
    site_outdir=OUTPUT_DIR+site+'/'
    if not os.path.exists(site_outdir):
        os.makedirs(site_outdir)
    #
    site_dict=UKsite_info.get_site_dict(site)
    sitefname=site_dict['fname']
    start_year=site_dict['years'][0]
    end_year=site_dict['years'][1]
    Jindex=site_dict['Jindex']
    PFT_index=site_dict['PFT_index']
    CF_n=site_dict['CF_n']
    Rd_n=site_dict['Rd_n']
    H_n=site_dict['H_n']
    LE_n=site_dict['LE_n']
    Site_int=site_dict['Site_int']
    CO2_range=site_dict['CO2_range']
    CO2_small_range=site_dict['CO2_close_range']
    CO2_smF=site_dict['CO2_smF']
    #
    # construct fluxsite filename
    if Site_int==1800.:
        fluxsite_fname = sitefname+'_'+str(start_year)[2:]+'-'+str(end_year)[2:]+'_halfhourly.nc'
    elif Site_int==3600.:
        fluxsite_fname = sitefname+'_'+str(start_year)[2:]+'-'+str(end_year)[2:]+'_hourly.nc'
    else:
        print 'Site_int =',Site_int,', Not recognised interval'
    
    # open file
    print 'opening site file: '+UKsites_DIR+fluxsite_fname
    SITE_inf = nc.Dataset(UKsites_DIR+fluxsite_fname,'r')
    #    
    # Read in flux site data
    site_FC=SITE_inf.variables[CF_n][:]
    site_FC_units=SITE_inf.variables[CF_n].units
    site_FC_longname=SITE_inf.variables[CF_n].long_name
    
    #if (Rd_n!=None):
    #    site_Rd=SITE_inf.variables[Rd_n][:]
    #    site_Rd_units=SITE_inf.variables[Rd_n].units
    #    site_Rd_longname=SITE_inf.variables[Rd_n].long_name
    
    #if (LE_n!=None):        
    #    site_LE=SITE_inf.variables[LE_n][:]
    #    site_LE_units=SITE_inf.variables[LE_n].units
    #    site_LE_longname=SITE_inf.variables[LE_n].long_name

    #if (H_n!=None):
    #    site_H=SITE_inf.variables[H_n][:]
    #    site_H_units=SITE_inf.variables[H_n].units
    #    site_H_longname=SITE_inf.variables[H_n].long_name

    site_time= nctime.num2date(SITE_inf.variables['time'][:], \
                               units=SITE_inf.variables['time'].units, \
                               calendar='standard')
    SITE_inf.close()
    
    # No 2013 model data, so chop off 2013 site data
    if (end_year>2012):
        end_year=2012
        end_time_stamp=nctime.num2date(0.,units='seconds since 2013-01-01 00:00:00')
        site_ep = np.where(site_time==end_time_stamp)[0]
        print 'site_ep = ',site_ep
        site_FC=site_FC[:site_ep+1]
        site_time=site_time[:site_ep+1]
    
    # loop round relevent JULES files and read in data
    # initiate numpy arrays to store data
    J_time=np.array([],dtype='float64')
    J_time_num=np.array([],dtype='float64')
    J_swdown=np.array([],dtype='float64')
    J_t1p5m=np.array([],dtype='float64')
    J_smc_tot=np.array([],dtype='float64')
    J_gpp=np.array([],dtype='float64')
    J_npp=np.array([],dtype='float64')
    J_resp_p=np.array([],dtype='float64')
    J_resp_s=np.array([],dtype='float64')
    J_anetc=np.array([],dtype='float64')
    J_rdc=np.array([],dtype='float64')
    
    for year in range(start_year,end_year+1):
        #construct filename
        JULES_fname='JULES_v4.3_TRIFFID_RsQ10.halfhourly.'+str(year)+'.nc'
        
        # openfile
        Jinf = nc.Dataset(JULES_output_DIR+JULES_fname,'r')
        
        # Read in data
        temp_J_time=Jinf.variables['time'][:]
        temp_J_time=np.zeros_like(temp_J_time,dtype='float64')+temp_J_time
        temp_J_time=np.round(temp_J_time/1800.)*1800.
        J_time=np.append(J_time,nctime.num2date(temp_J_time,             \
                                                units=Jinf.variables['time'].units,    \
                                                calendar=Jinf.variables['time'].calendar)  )
        J_time_num=np.append(J_time_num,temp_J_time)
        #
        #J_swdown=np.append(J_swdown,Jinf.variables['sw_down'][:,0,Jindex])
        #J_t1p5m=np.append(J_t1p5m,Jinf.variables['t1p5m'][:,PFT_index,0,Jindex])
        #J_smc_tot=np.append(J_smc_tot,Jinf.variables['smc_tot'][:,0,Jindex])
        J_gpp=np.append(J_gpp,Jinf.variables['gpp'][:,PFT_index,0,Jindex])
        J_npp=np.append(J_npp,Jinf.variables['npp'][:,PFT_index,0,Jindex])
        J_resp_p=np.append(J_resp_p,Jinf.variables['resp_p'][:,PFT_index,0,Jindex])
        J_resp_s=np.append(J_resp_s,np.sum(Jinf.variables['resp_s'][:,:,0,Jindex],axis=1))
        J_anetc=np.append(J_anetc,Jinf.variables['anetc'][:,PFT_index,0,Jindex])
        J_rdc=np.append(J_rdc,Jinf.variables['rdc'][:,PFT_index,0,Jindex])
        #
        Jinf.close()
    #
    
    # Put Site data onto J grid
    # first find start and end points using the time arrays
    J_sp = np.where(J_time==site_time[0])[0]
    J_ep = np.where(J_time==site_time[-1])[0]
    #
    if (Site_int==1800.):
        site_to_Jgrid_index=np.arange(J_sp[0],J_ep[0]+1,1)
    elif (Site_int==3600.):
        # for hourly data only fill every other point
        site_to_Jgrid_index=np.arange(J_sp[0],J_ep[0]+1,2)
    #
    nDays = int(len(J_time)/48.)
    #
    # create diurnal time arrays.
    DIURNAL_time_hours = (np.round(J_time_num/1800.)/2.) % 24.
    DIURNAL_time_plot  = np.append(DIURNAL_time_hours[:47],24.)
    #
    # create seasonal indexes
    month_arr = np.array([dt.month for dt in J_time])
    seasonal_indexes = [ np.where((month_arr==12)|(month_arr==1)|(month_arr==2)), \
                         np.where((month_arr==3)|(month_arr==4)|(month_arr==5)),  \
                         np.where((month_arr==6)|(month_arr==7)|(month_arr==8)),  \
                         np.where((month_arr==10)|(month_arr==11)|(month_arr==9)) ]
    seasonal_names = [ 'DJF','MAM','JJA','SON' ]
    month_DIURNAL_grid   = month_arr.reshape(nDays,48)[:,0]
    
    site_FC_Jgrid = np.zeros_like(J_npp)+site_FC.fill_value
    site_FC_Jgrid[site_to_Jgrid_index]=site_FC
    site_FC_Jgrid=np.ma.masked_equal(site_FC_Jgrid,site_FC.fill_value)
    site_FC_index=np.where(site_FC_Jgrid.mask==False)
    
    J_NEE_homemade = J_gpp-J_resp_p-J_resp_s
    
    J_time_DIURNAL_grid  = J_time.reshape(nDays,48)
    J_time_num_DIURNAL_grid  = J_time_num.reshape(nDays,48)
    J_npp_DIURNAL_grid   = J_npp.reshape(nDays,48)
    J_gpp_DIURNAL_grid   = J_gpp.reshape(nDays,48)
    J_anetc_DIURNAL_grid = J_anetc.reshape(nDays,48)
    J_rdc_DIURNAL_grid   = J_rdc.reshape(nDays,48)
    J_resp_p_DIURNAL_grid   = J_resp_p.reshape(nDays,48)
    J_resp_s_DIURNAL_grid   = J_resp_s.reshape(nDays,48)
    J_NEE_homemade_DIURNAL_grid  = J_NEE_homemade.reshape(nDays,48)
    
    site_FC_DIURNAL_grid = site_FC_Jgrid.reshape(nDays,48)   
    
    site_FC_dailymean        = np.mean(site_FC_DIURNAL_grid,axis=1)
    J_NEE_dailymean          = np.mean(J_NEE_homemade_DIURNAL_grid,axis=1)
    J_npp_dailymean          = np.mean(J_npp_DIURNAL_grid,axis=1)
    J_gpp_dailymean          = np.mean(J_gpp_DIURNAL_grid,axis=1)
    J_resp_p_dailymean       = np.mean(J_resp_p_DIURNAL_grid,axis=1)
    J_resp_s_dailymean       = np.mean(J_resp_s_DIURNAL_grid,axis=1)
    
    # Smoothed Daily means
    site_FC_dailymean_smooth  = MT.FFT_smooth(site_FC_dailymean,1,smooth_factor=CO2_smF*2.)
    J_NEE_dailymean_smooth    = MT.FFT_smooth(J_NEE_dailymean,1,smooth_factor=CO2_smF*2.)
    J_npp_dailymean_smooth    = MT.FFT_smooth(J_npp_dailymean,1,smooth_factor=CO2_smF*2.)
    J_gpp_dailymean_smooth    = MT.FFT_smooth(J_gpp_dailymean,1,smooth_factor=CO2_smF*2.)
    J_resp_p_dailymean_smooth = MT.FFT_smooth(J_resp_p_dailymean,1,smooth_factor=CO2_smF*2.)
    J_resp_s_dailymean_smooth = MT.FFT_smooth(J_resp_s_dailymean,1,smooth_factor=CO2_smF*2.)
    
    
    # Calculate climatologies
    J_date=np.array([ J_time_DIURNAL_grid[i,0].day for i in range(nDays) ] )
    J_month=np.array([ J_time_DIURNAL_grid[i,0].month for i in range(nDays) ] )
    SEASONAL_time_plot = J_time_DIURNAL_grid[:365,0]
    no_leap_index=np.where((J_month!=2)|(J_date!=29))[0]
    nYears=end_year-start_year+1
    site_FC_climatology = np.mean( site_FC_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    J_NEE_climatology = np.mean( J_NEE_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    J_npp_climatology = np.mean( J_npp_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    J_gpp_climatology = np.mean( J_gpp_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    J_resp_p_climatology = np.mean( J_resp_p_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    J_resp_s_climatology = np.mean( J_resp_s_dailymean[no_leap_index].reshape(nYears,365), axis=0)
    
    # Calculate de-seasonalised trend
    site_FC_anomaly = np.array( [ site_FC_dailymean[no_leap_index][i:i+365] -
                                  site_FC_climatology for i in range(0,nYears*365,365) ] ).flatten()
    J_NEE_anomaly   = np.array( [ J_NEE_dailymean[no_leap_index][i:i+365] -
                                  J_NEE_climatology for i in  range(0,nYears*365,365) ] ).flatten()
    J_npp_anomaly   = np.array( [ J_npp_dailymean[no_leap_index][i:i+365] -
                                  J_npp_climatology for i in  range(0,nYears*365,365) ] ).flatten()
    J_gpp_anomaly   = np.array( [ J_gpp_dailymean[no_leap_index][i:i+365] -
                                  J_gpp_climatology for i in  range(0,nYears*365,365) ] ).flatten()
    J_resp_p_anomaly= np.array( [ J_resp_p_dailymean[no_leap_index][i:i+365] -
                                  J_resp_p_climatology for i in  range(0,nYears*365,365) ] ).flatten()
    J_resp_s_anomaly= np.array( [ J_resp_s_dailymean[no_leap_index][i:i+365] -
                                  J_resp_s_climatology for i in  range(0,nYears*365,365) ] ).flatten()
    
    # FFT smoothed de-seasonalised trend
    site_FC_anomaly_smooth  = MT.FFT_smooth(site_FC_anomaly,1,smooth_factor=CO2_smF/2.)
    J_NEE_anomaly_smooth    = MT.FFT_smooth(J_NEE_anomaly,1,smooth_factor=CO2_smF/2.)
    J_npp_anomaly_smooth    = MT.FFT_smooth(J_npp_anomaly,1,smooth_factor=CO2_smF/2.)
    J_gpp_anomaly_smooth    = MT.FFT_smooth(J_gpp_anomaly,1,smooth_factor=CO2_smF/2.)
    J_resp_p_anomaly_smooth = MT.FFT_smooth(J_resp_p_anomaly,1,smooth_factor=CO2_smF/2.)
    J_resp_s_anomaly_smooth = MT.FFT_smooth(J_resp_s_anomaly,1,smooth_factor=CO2_smF/2.)
    
    ##################################################################################################
    #  PLOT SECTION
    ##################################################################################################
    
    # Plot Full Time-Series multi-plot (Daily means)
    fig=plt.figure(figsize=(13,10))
    AX=fig.add_subplot(2,1,1)
    AX.plot(J_time_DIURNAL_grid[:,0],site_FC_dailymean*-1,color='black',lw=1.5,label='NEE')
    AX.set_ylim([-10,20])
    AX.legend(loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('Observations')

    AX=fig.add_subplot(2,1,2)
    AX.plot(J_time_DIURNAL_grid[:,0],J_npp_dailymean*kgC2umol,\
            color='yellow',lw=1.,label='npp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_gpp_dailymean*kgC2umol,\
            color='green',lw=1.,label='gpp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_p_dailymean*kgC2umol*-1,\
            color='darkorange',lw=1.,label='resp_p')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_s_dailymean*kgC2umol*-1,\
            color='brown',lw=1.,label='resp_s')
    AX.plot(J_time_DIURNAL_grid[:,0],J_NEE_dailymean*kgC2umol,\
            color='red',lw=1.5,label='NEE')
    AX.set_ylim(CO2_range)
    AX.legend(ncol=5,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('JULES')
    #plt.show()
    plt.savefig(site_outdir+'NEE_fulltimeseries_multi.png',bbox_inches='tight')
    plt.close()

    # Plot Full Time-Series (Daily means)
    fig=plt.figure(figsize=(15,7))
    AX=fig.add_subplot(1,1,1)
    AX.plot(J_time_DIURNAL_grid[:,0],site_FC_dailymean*-1,color='black',lw=1.,label='NEE - Obs')

    AX.plot(J_time_DIURNAL_grid[:,0],J_npp_dailymean*kgC2umol,\
            color='yellow',lw=0.7,label='npp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_gpp_dailymean*kgC2umol,\
            color='green',lw=0.7,label='gpp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_p_dailymean*kgC2umol*-1,\
            color='darkorange',lw=0.7,label='resp_p')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_s_dailymean*kgC2umol*-1,\
            color='brown',lw=0.7,label='resp_s')
    AX.plot(J_time_DIURNAL_grid[:,0],J_NEE_dailymean*kgC2umol,\
            color='red',lw=1.5,label='NEE')
    AX.set_ylim(CO2_range)
    AX.legend(ncol=6,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('NEE time-series - Observations vs JULES')
    #plt.show()
    plt.savefig(site_outdir+'NEE_fulltimeseries.png',bbox_inches='tight')
    plt.close()
    
    # Plot Smoothed Full Time-Series (Daily means)
    fig=plt.figure(figsize=(15,7))
    AX=fig.add_subplot(1,1,1)
    AX.plot(J_time_DIURNAL_grid[:,0],site_FC_dailymean_smooth*-1,\
            color='black',lw=1.5,label='NEE - Obs')

    AX.plot(J_time_DIURNAL_grid[:,0],J_npp_dailymean_smooth*kgC2umol,\
            color='yellow',lw=0.7,label='npp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_gpp_dailymean_smooth*kgC2umol,\
            color='green',lw=0.7,label='gpp')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_p_dailymean_smooth*kgC2umol*-1,\
            color='darkorange',lw=0.7,label='resp_p')
    AX.plot(J_time_DIURNAL_grid[:,0],J_resp_s_dailymean_smooth*kgC2umol*-1,\
            color='brown',lw=0.7,label='resp_s')
    AX.plot(J_time_DIURNAL_grid[:,0],J_NEE_dailymean_smooth*kgC2umol,\
            color='red',lw=1.5,label='NEE')
    AX.set_ylim(CO2_small_range)
    AX.legend(ncol=6,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('Smoothed NEE time-series - Observations vs JULES')
    #plt.show()
    plt.savefig(site_outdir+'NEE_fulltimeseries_smoothed.png',bbox_inches='tight')
    plt.close()
    
    # Plot climatolgy of daily means
    fig=plt.figure(figsize=(15,7))
    AX=fig.add_subplot(1,1,1)
    AX.plot(SEASONAL_time_plot,site_FC_climatology*-1,color='black',lw=1.5,label='NEE Obs')
    AX.plot(SEASONAL_time_plot,J_npp_climatology*kgC2umol,\
            color='yellow',lw=0.7,label='npp')
    AX.plot(SEASONAL_time_plot,J_gpp_climatology*kgC2umol,\
            color='green',lw=0.7,label='gpp')
    AX.plot(SEASONAL_time_plot,J_resp_p_climatology*kgC2umol*-1,\
            color='darkorange',lw=0.7,label='resp_p')
    AX.plot(SEASONAL_time_plot,J_resp_s_climatology*kgC2umol*-1,\
            color='brown',lw=0.7,label='resp_s')
    AX.plot(SEASONAL_time_plot,J_NEE_climatology*kgC2umol,\
            color='red',lw=1.5,label='NEE')
    AX.set_ylim(CO2_range)
    AX.legend(ncol=6,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('Climatology - Obs vs JULES')
    
    #plt.show()
    plt.savefig(site_outdir+'NEE_seasonality.png',bbox_inches='tight')
    plt.close()
    
    # Plot deseasonalised trend
    fig=plt.figure(figsize=(15,7))
    AX=fig.add_subplot(1,1,1)
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],site_FC_anomaly*-1,\
            color='black',lw=1.,label='NEE - Obs')

    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_npp_anomaly*kgC2umol,\
            color='yellow',lw=1.,label='npp')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_gpp_anomaly*kgC2umol,\
            color='green',lw=1.,label='gpp')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_resp_p_anomaly*kgC2umol*-1,\
            color='darkorange',lw=1.,label='resp_p')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_resp_s_anomaly*kgC2umol*-1,\
            color='brown',lw=1.,label='resp_s')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_NEE_anomaly*kgC2umol,\
            color='red',lw=2.,label='NEE')

    AX.set_ylim(CO2_range)
    AX.legend(ncol=6,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('NEE Anomaly time-series - Observations vs JULES')
    #plt.show()
    plt.savefig(site_outdir+'NEE_anomaly_fulltimeseries.png',bbox_inches='tight')
    plt.close()
        
    # Plot Smoothed deseasonalised trend
    fig=plt.figure(figsize=(15,7))
    AX=fig.add_subplot(1,1,1)
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],site_FC_anomaly_smooth*-1,\
            color='black',lw=1.5,label='NEE - Obs')

    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_npp_anomaly_smooth*kgC2umol,\
            color='yellow',lw=1.,label='npp')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_gpp_anomaly_smooth*kgC2umol,\
            color='green',lw=1.,label='gpp')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_resp_p_anomaly_smooth*kgC2umol*-1,\
            color='darkorange',lw=1.,label='resp_p')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_resp_s_anomaly_smooth*kgC2umol*-1,\
            color='brown',lw=1.,label='resp_s')
    AX.plot(J_time_DIURNAL_grid[no_leap_index,0],J_NEE_anomaly_smooth*kgC2umol,\
            color='red',lw=2.,label='NEE')

    AX.set_ylim([-5,5])
    AX.legend(ncol=6,loc='upper center')
    AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
    AX.set_title('NEE Smoothed Anomaly time-series - Observations vs JULES')
    #plt.show()
    plt.savefig(site_outdir+'NEE_anomaly_fulltimeseries_smoothed.png',bbox_inches='tight')
    plt.close()
    
    
    # Plot Seasonal Diurnal Cycles
    fig=plt.figure(figsize=(15,12))
    for seas in range(len(seasonal_indexes)):
        seas_index=seasonal_indexes[seas][0]
        DIURNAL_time_hours_seas = DIURNAL_time_hours[seas_index]
        J_npp_seas_data=J_npp[seas_index]
        J_gpp_seas_data=J_gpp[seas_index]
        J_anetc_seas_data=J_anetc[seas_index]
        J_rdc_seas_data=J_rdc[seas_index]
        J_resp_p_seas_data=J_resp_p[seas_index]
        J_resp_s_seas_data=J_resp_s[seas_index]
        
        J_NEE_homemade_seas = J_NEE_homemade[seas_index]
        
        site_FC_seas_data=site_FC_Jgrid[seas_index]
        site_FC_seas_index= np.where(site_FC_seas_data.mask==False)
        
        nSeasDays=int(len(seas_index)/48.)
        
        J_npp_seas_DIURNAL_grid   = J_npp_seas_data.reshape(nSeasDays,48)
        J_gpp_seas_DIURNAL_grid   = J_gpp_seas_data.reshape(nSeasDays,48)
        J_anetc_seas_DIURNAL_grid = J_anetc_seas_data.reshape(nSeasDays,48)
        J_rdc_seas_DIURNAL_grid   = J_rdc_seas_data.reshape(nSeasDays,48)
        J_resp_p_seas_DIURNAL_grid   = J_resp_p_seas_data.reshape(nSeasDays,48)
        J_resp_s_seas_DIURNAL_grid   = J_resp_s_seas_data.reshape(nSeasDays,48)
        J_NEE_homemade_seas_DIURNAL_grid  = J_NEE_homemade_seas.reshape(nSeasDays,48)
        
        site_FC_seas_DIURNAL_grid = site_FC_seas_data.reshape(nSeasDays,48)
        
        AX=fig.add_subplot(4,2,(seas*2)+1)
        AX.hist2d(DIURNAL_time_hours_seas[site_FC_seas_index],\
                  site_FC_seas_data[site_FC_seas_index]*-1,\
                  (24,50), cmin=1,cmax=1000, cmap='OrRd', \
                  range=[ [0,24], CO2_range ]  )
        
        AX.plot(DIURNAL_time_plot,np.mean(site_FC_seas_DIURNAL_grid*-1,axis=0),color='red',lw=3)
        AX.set_ylim(CO2_range)
        AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
        AX.set_title(seasonal_names[seas])
        
        AX=fig.add_subplot(4,2,(seas*2)+2)
        AX.hist2d(DIURNAL_time_hours_seas,\
                  (J_NEE_homemade_seas*kgC2umol),\
                  (24,50), cmin=1,cmax=1000, cmap='OrRd', \
                  range=[ [0,24], CO2_range ]  )
        
        
        AX.plot(DIURNAL_time_plot,np.mean(J_npp_seas_DIURNAL_grid*kgC2umol,axis=0),\
                color='yellow',lw=2)
        
        AX.plot(DIURNAL_time_plot,np.mean(J_gpp_seas_DIURNAL_grid*kgC2umol,axis=0),\
                color='green',lw=2)
        #AX.plot(DIURNAL_time_plot,np.mean(J_anetc_seas_DIURNAL_grid*kgC2umol,axis=0),\
        #        color='yellow',lw=2)
        #AX.plot(DIURNAL_time_plot,np.mean(J_rdc_seas_DIURNAL_grid*kgC2umol*-1,axis=0),\
        #        color='black',lw=2)
        AX.plot(DIURNAL_time_plot,np.mean(J_resp_p_seas_DIURNAL_grid*kgC2umol*-1,axis=0),\
                color='orange',lw=2)
        AX.plot(DIURNAL_time_plot,np.mean(J_resp_s_seas_DIURNAL_grid*kgC2umol*-1,axis=0),\
                color='brown',lw=2)
        
        AX.plot(DIURNAL_time_plot,np.mean(J_NEE_homemade_seas_DIURNAL_grid*kgC2umol*1,axis=0),\
                color='red',lw=3)
        
        AX.set_ylim(CO2_range)
        AX.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
        AX.set_title(seasonal_names[seas])
        
    #plt.show()
    plt.savefig(site_outdir+'NEE_DiurnalCycles.png',bbox_inches='tight')
    plt.close()
    
    
    #plt.plot(np.arange(nDays),np.mean(site_FC_DIURNAL_grid,axis=1),color='red')
    #plt.plot(np.arange(nDays),np.mean(J_npp_DIURNAL_grid*1e8,axis=1),color='blue')
    #plt.plot(np.arange(nDays),np.mean(J_gpp_DIURNAL_grid*1e8,axis=1),color='green')
    #plt.show()

#fig=plt.figure()
#ax=fig.add_subplot(1,1,1)
#ax.plot(np.arange(10)) 
#ax.set_ylabel('$\mu$mol $CO_2$ $m^2$ $s^-1$')
#plt.show()

