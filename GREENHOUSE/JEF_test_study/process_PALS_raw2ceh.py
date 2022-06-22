#!/bin/python
#
# Script to process raw PALS data to the commeon ceh format
#        also calculates the GPP based on subtraction of the 
#        mean NEE when SWdown=0 for the encompassing 24 hours
#

import numpy as np
import datetime as dt
import netCDF4 as nc
import netcdftime as nctime
import netcdf_tools as nctools
import os
import maths_tools.DateTimeTools as DTT
import maths_tools.MathsTools as MT
from scipy.optimize import curve_fit
import scipy.ndimage as im

#############################################
# Parameters and Options
met_params=['SWdown','LWdown','PSurf','Qair',\
            'Rainf','Tair','Wind']
flux_params=['NEE','Qg','Qh','Qle','Rnet','GPP','DarkResp']
fill_value=-9999.
nTRES = 4  # number of time resolutions
GPP_CW = 96 # size of GPP calculation window in tsteps
half_GPP_CW = int(GPP_CW/2.) # half of above for simpler indexing
#nat_tres = 0.5   # native time resolution in hours
#increments = [ 0, dt.timedelt

#####################################
# Directories and filenames
in_dir  = '/data/grp/fluxdata/raw/GLOBALFLUXDATA/'
out_dir = '/data/grp/fluxdata/PALS_sites_ECP/'

site_list_file = in_dir+'site_list.txt'

####################################
# read list of sites
site_list = open(site_list_file,'r').readlines()
site_list=filter(None,site_list)
site_list=[site.replace('\n','') for site in site_list]

# Loop round each site and process data
for site in [site_list[3]]:   #  site_list:
    print 'Processing: '+site+'\n'
    in_met_file       = in_dir+site+'Fluxnet.1.4_met.nc'
    in_flux_file      = in_dir+site+'Fluxnet.1.4_flux.nc'
    out_met_file_tag  = out_dir+site+'Fluxnet.1.4_met'
    out_flux_file_tag = out_dir+site+'Fluxnet.1.4_flux'
    # First check that there are met and flux files for the site
    #if not (os.path.isfile(in_met_file)&os.path.isfile(in_flux_file)):
    #    continue
    
    # Read in Met Data 
    # aggregate to different time resolutions
    # and write out to appropriate files
    
    # open in_met_file
    minf=nc.Dataset(in_met_file,'r')
    
    # read in time data
    tunits =  minf.variables['time'].units
    in_datetime = nctime.num2date(minf.variables['time'][:], units=tunits )
    #in_date = np.array([date.date() for date in in_datetime])
    
    # extract start/end time/day/month
    start_time = in_datetime[0]
    end_time   = in_datetime[-1]
    start_day  = dt.datetime(start_time.year,start_time.month,start_time.day)
    end_day    = dt.datetime(end_time.year,end_time.month,end_time.day)
    start_month= dt.datetime(start_time.year,start_time.month,1)
    end_month  = dt.datetime(end_time.year,end_time.month,1)
    
    # construct  arrays of time objects and numbers for the different Tres
    time_nat_obj = in_datetime
    time_nat_num = minf.variables['time'][:]
    time_3hr_obj = np.array(DTT.DTarange(start_time,end_time,tdelta=dt.timedelta(hours=3)))
    time_3hr_num = nctime.date2num(time_3hr_obj,units=tunits)
    time_day_obj = np.array(DTT.DTarange(start_day,end_day))
    time_day_num = nctime.date2num(time_day_obj,units=tunits)
    time_month_obj = np.array(DTT.DTarange_months(start_month,end_month))
    time_month_num = nctime.date2num(time_month_obj,units=tunits)
    # Append time obj/num arrays to list
    time_obj_list = [ time_nat_obj,time_3hr_obj,time_day_obj,time_month_obj ]
    time_num_list = [ time_nat_num,time_3hr_num,time_day_num,time_month_num ]
    time_len_list = [ len(temp_time) for temp_time in time_obj_list ]
    
    # build a list of lists of indexes for averaging in subsequent sections
    # Blank first element for the native resolution where no averaging takes place
    time_indices = [ None ]
    for iTRES in range(1,nTRES):
        # select time object array from list
        temp_time_obj=time_obj_list[iTRES]
        # for each tstep create index of where native time is in between 
        # current and subsequent steps
        temp_time_index = [   (time_nat_obj>=temp_time_obj[iSTEP])  \
                               & (time_nat_obj<temp_time_obj[iSTEP+1]) \
                                 for iSTEP in range(len(temp_time_obj)-1) ]
        # append greater than last time step
        temp_time_index.append( time_nat_obj>=temp_time_obj[-1] )
        time_indices.append( temp_time_index )

    # open out files
    moutfiles = [ \
                  nc.Dataset(out_met_file_tag+'.nc','w'),            \
                  nc.Dataset(out_met_file_tag+'_threehourly.nc','w'),\
                  nc.Dataset(out_met_file_tag+'_daily.nc','w'),      \
                  nc.Dataset(out_met_file_tag+'_monthly.nc','w'),    \
                  ]
    
    # Create Dimensions for outfiles
    for iTRES in range(nTRES):
        moutfiles[iTRES].createDimension('land',1)
        moutfiles[iTRES].createDimension('time',time_len_list[iTRES])
    
    # Copy lat/lon and time variables and other static variables
    for iTRES in range(nTRES):
        # lat/lon
        nctools.copy_variable(minf,moutfiles[iTRES],'latitude',dimensions=('land'))
        nctools.copy_variable(minf,moutfiles[iTRES],'longitude',dimensions=('land'))
        # time
        outvar=moutfiles[iTRES].createVariable('time','float32',('time'))
        outvar.units=tunits
        outvar[:]=time_num_list[iTRES]
        # elevation and reference height
        nctools.copy_variable(minf,moutfiles[iTRES],'elevation',dimensions=('land'))
        nctools.copy_variable(minf,moutfiles[iTRES],'reference_height',dimensions=('land'))

    
    for param in met_params:
        data = minf.variables[param][:].squeeze()
        mv   = minf.variables[param].missing_value
        # Don't use the met data flags cos they are crap
        #flag = minf.variables[param+'_qc'][:].squeeze()
        #mask = (flag==0)|(data==mv)
        #data_masked = np.ma.masked_array(data,mask=mask,fill_value=fill_value)
        #data_masked.data[mask==True] = fill_value
        data_masked=data
        if param=='SWdown':
            SWdown_data = data_masked.copy()
        elif param=='Tair':
            Tair_data = data_masked.copy()


        moutvars = []
        for moutf in moutfiles:
            moutvar=moutf.createVariable(param,'float32',\
                                                  ('time','land'),\
                                                  fill_value=fill_value)
            moutvars.append( moutvar )
            for attr in minf.variables[param].ncattrs():
                moutvar.setncattr(str(attr),minf.variables[param].getncattr(str(attr)))
                
        # For native resolution (moutfiles[0]) copy straight over
        moutvars[0][:]=data_masked
        
        for iTRES in range(1,nTRES):
            time_indexes = time_indices[iTRES]
            moutf = moutfiles[iTRES]
            
            # create list of indexed data
            data_masked_Tindexed = [ data_masked[Tindex] for Tindex in time_indexes ]

            mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ]) 
            moutvars[iTRES][:] = mean_data
        
        for iTRES in range(1,nTRES):
            time_indexes = time_indices[iTRES]
            moutf = moutfiles[iTRES]
            
            # create list of indexed data
            data_masked_Tindexed = [ data_masked[Tindex] for Tindex in time_indexes ]

            #####################################################################
            # Loop round the time indexes and calculate the mean, median, 
            #  stddev 
            Npts_data   = np.array([ len(dat) for dat in data_masked_Tindexed ])
            mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ])
            median_data = np.array([ np.median(dat) for dat in data_masked_Tindexed ])
            stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
            # stderr_data = stddev_data/Npts_data
            
            # Store the mean into the pre created variable 
            moutvars[iTRES][:] = mean_data
            moutvars[iTRES].long_name+=' - Mean'

            # Create variables for remaining parameters
            # Median
            moutvar=moutf.createVariable(param+'.median','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=minf.variables[param].units
            moutvar.long_name=minf.variables[param].long_name+' - Median'
            # StdDev
            moutvar=moutf.createVariable(param+'.stddev','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=minf.variables[param].units
            moutvar.long_name=minf.variables[param].long_name+' - Standard Deviation'
            # Npts
            moutvar=moutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=''
            moutvar.long_name=minf.variables[param].long_name+' - Number of Averaged time-steps'
        
        del data_masked
        del data_masked_Tindexed
        del mean_data

    for iTRES in range(nTRES):
        for attr in minf.ncattrs():
            moutfiles[iTRES].setncattr( str(attr), minf.getncattr(str(attr)) )
        moutfiles[iTRES].close()

    minf.close()
    
    
     
    #####################################################################
    # Now repeat process for fluxes, 
    #   For fluxes include the median, std and npts
    # open in_flux_file
    finf=nc.Dataset(in_flux_file,'r')
    
    # read in time data - REPEATING THIS PROCESS MAY BE OVERLY SENSITVE 
    #                     BUT IT GUARENTEES CORRECT TREATMENT OF THE DATA
    tunits =  finf.variables['time'].units
    in_datetime = nctime.num2date(finf.variables['time'][:], units=tunits )
    
    # extract start/end time/day/month
    start_time = in_datetime[0]
    end_time   = in_datetime[-1]
    start_day  = dt.datetime(start_time.year,start_time.month,start_time.day)
    end_day    = dt.datetime(end_time.year,end_time.month,end_time.day)
    start_month= dt.datetime(start_time.year,start_time.month,1)
    end_month  = dt.datetime(end_time.year,end_time.month,1)
    
    # construct  arrays of time objects and numbers for the different Tres
    time_nat_obj = in_datetime
    time_nat_num = finf.variables['time'][:]
    time_3hr_obj = np.array(DTT.DTarange(start_time,end_time,tdelta=dt.timedelta(hours=3)))
    time_3hr_num = nctime.date2num(time_3hr_obj,units=tunits)
    time_day_obj = np.array(DTT.DTarange(start_day,end_day))
    time_day_num = nctime.date2num(time_day_obj,units=tunits)
    time_month_obj = np.array(DTT.DTarange_months(start_month,end_month))
    time_month_num = nctime.date2num(time_month_obj,units=tunits)
    # Append time obj/num arrays to list
    time_obj_list = [ time_nat_obj,time_3hr_obj,time_day_obj,time_month_obj ]
    time_num_list = [ time_nat_num,time_3hr_num,time_day_num,time_month_num ]
    time_len_list = [ len(temp_time) for temp_time in time_obj_list ]
    
    # build a list of lists of indexes for averaging in subsequent sections
    # Blank first element for the native resolution where no averaging takes place
    time_indices = [ None ]
    for iTRES in range(1,nTRES):
        # select time object array from list
        temp_time_obj=time_obj_list[iTRES]
        # for each tstep create index of where native time is in between 
        # current and subsequent steps
        temp_time_index = [   (time_nat_obj>=temp_time_obj[iSTEP])  \
                               & (time_nat_obj<temp_time_obj[iSTEP+1]) \
                                 for iSTEP in range(len(temp_time_obj)-1) ]
        # append greater than last time step
        temp_time_index.append( time_nat_obj>=temp_time_obj[-1] )
        time_indices.append( temp_time_index )


    
    foutfiles = [ \
                  nc.Dataset(out_flux_file_tag+'.nc','w'),            \
                  nc.Dataset(out_flux_file_tag+'_threehourly.nc','w'),\
                  nc.Dataset(out_flux_file_tag+'_daily.nc','w'),      \
                  nc.Dataset(out_flux_file_tag+'_monthly.nc','w'),    \
                  ]
    
    # Create Dimensions for outfiles
    for iTRES in range(nTRES):
        foutfiles[iTRES].createDimension('land',1)
        foutfiles[iTRES].createDimension('time',time_len_list[iTRES])
    
    # Copy lat/lon and time variables
    for iTRES in range(nTRES):
        # lat/lon
        nctools.copy_variable(finf,foutfiles[iTRES],'latitude',dimensions=('land'))
        nctools.copy_variable(finf,foutfiles[iTRES],'longitude',dimensions=('land'))
        # time
        outvar=foutfiles[iTRES].createVariable('time','float32',('time'))
        outvar.units=tunits
        outvar[:]=time_num_list[iTRES]
        # elevation and reference height
        nctools.copy_variable(finf,foutfiles[iTRES],'elevation',dimensions=('land'))
        nctools.copy_variable(finf,foutfiles[iTRES],'reference_height',dimensions=('land'))
        
    # loop round flux parameters (except Dark Respiration and GPP as we have to calculate these)
    for param in flux_params[:-2]:
        data = finf.variables[param][:].squeeze()
        mv   = finf.variables[param].missing_value
        flag = finf.variables[param+'_qc'][:].squeeze()
        mask = (flag==0)|(data==mv)
        data_masked = np.ma.masked_array(data,mask=mask,fill_value=fill_value)
        data_masked.data[mask==True] = fill_value
        # Store NEE data for GPP and Dark Resp calcualation below
        if param=='NEE':
            NEE_data      = data_masked.copy()
        
        foutvars = []
        for foutf in foutfiles:
            foutvar=foutf.createVariable(param,'float32',\
                                                  ('time','land'),\
                                                  fill_value=fill_value)
            foutvars.append( foutvar )
            for attr in finf.variables[param].ncattrs():
                foutvar.setncattr(str(attr),finf.variables[param].getncattr(str(attr)))
                
        # For native resolution (foutfiles[0]) copy straight over
        foutvars[0][:]=data_masked
        
        for iTRES in range(1,nTRES):
            time_indexes = time_indices[iTRES]
            foutf = foutfiles[iTRES]
            
            # create list of indexed data
            data_masked_Tindexed = [ data_masked[Tindex] for Tindex in time_indexes ]

            #####################################################################
            # Loop round the time indexes and calculate the mean, median, 
            #  stddev 
            Npts_data   = np.array([ len(dat) for dat in data_masked_Tindexed ])
            mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ])
            median_data = np.array([ np.median(dat) for dat in data_masked_Tindexed ])
            stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
            # stderr_data = stddev_data/Npts_data
            
            # Store the mean into the pre created variable 
            foutvars[iTRES][:] = mean_data
            foutvars[iTRES].long_name+=' - Mean'

            # Create variables for remaining parameters
            # Median
            foutvar=foutf.createVariable(param+'.median','float32',('time','land'),\
                                          fill_value=fill_value                     )
            foutvar.units=finf.variables[param].units
            foutvar.long_name=finf.variables[param].long_name+' - Median'
            # StdDev
            foutvar=foutf.createVariable(param+'.stddev','float32',('time','land'),\
                                          fill_value=fill_value                     )
            foutvar.units=finf.variables[param].units
            foutvar.long_name=finf.variables[param].long_name+' - Standard Deviation'
            # Npts
            foutvar=foutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                          fill_value=fill_value                     )
            foutvar.units=''
            foutvar.long_name=finf.variables[param].long_name+' - Number of Averaged time-steps'


    # Calculate Dark respiration and  GPP
    # Dark Respiration is calculated by fitting the Night-time data
    #    to a Q10 function
    # Night-time data is defined as data which is not within 2 timesteps
    #    of registered SWdown
     
    # Daytime mask for calculation of GPP
    DAYTIME_mask= SWdown_data>10
    DAYTIME_mask_copy=DAYTIME_mask.copy()
    # add points within 2 timesteps of daytime to mask
    # copy mask to prevent masking all elements
    SW_smooth_int=3
    for iSMOOTH in range(SW_smooth_int,len(DAYTIME_mask)-SW_smooth_int):
        if any(DAYTIME_mask_copy[iSMOOTH-SW_smooth_int:iSMOOTH+SW_smooth_int]==True):
            DAYTIME_mask[iSMOOTH]=True
    del DAYTIME_mask_copy
    
    # Remove NEE spikes with scipy.ndimage.median_filter
    NEE_median_filter_size=4
    temp = NEE_data.data
    temp[NEE_data.mask==True]=np.nan
    NEE_data_filtered = im.median_filter(temp,NEE_median_filter_size)
    # Reapply mask
    NEE_data_filtered[NEE_data.mask==True]=fill_value
    # Mask out additional elements which became nan during filtering
    NEE_data_filtered[NEE_data_filtered!=NEE_data_filtered]=fill_value
    NEE_data_filtered=np.ma.masked_equal(NEE_data_filtered,fill_value)
    del temp
    
    dark_mask     = DAYTIME_mask|NEE_data_filtered.mask
    dark_NEE_data = np.ma.masked_array(NEE_data_filtered,mask=dark_mask,fill_value=fill_value)
    dark_NEE_data.data[dark_mask==True]=fill_value
    dark_Tair_data= np.ma.masked_array(Tair_data,mask=dark_mask,fill_value=fill_value)
    dark_Tair_data.data[dark_mask==True]=fill_value
    
    # group data into temprature bins
    bin_temp = 2.  
    min_temp = np.floor(dark_Tair_data.min())
    max_temp = np.ceil(dark_Tair_data.max())
    dark_Tair_lims = np.arange(min_temp,max_temp+.1,bin_temp)
    dark_Tair_bins = np.arange(min_temp,max_temp,bin_temp)+(bin_temp/2.)
    dark_NEE_binmean = np.zeros_like(dark_Tair_bins)+fill_value
    nbins = len(dark_Tair_bins)
    
    for iTEMP in range(nbins):
        dark_NEE_binmean[iTEMP] = \
                np.mean( dark_NEE_data[ (dark_Tair_data>=dark_Tair_lims[iTEMP])  & \
                                        (dark_Tair_data<dark_Tair_lims[iTEMP+1]) & \
                                        (dark_NEE_data.mask==False) ] )
        
        
    
    # Optimise the Q10 fit parameter for local conditions by fitting:
    #    DarkResp = Q10^( m(Tair+c) ), 
    # the parameter fit is performed on rolling time frame covering the encompassing 48 hours
    dark_Resp_data = 0.-dark_NEE_data
    Resp_data = np.zeros_like(NEE_data_filtered.data)+fill_value
    
    
    
    Flambda = lambda Tair,Q_10,Resp_10: MT.Q10_func(Tair,Q_10,Resp_10,0.1,273.15)
    popt,pcov=curve_fit(Flambda,\
                        dark_Tair_data,dark_Resp_data,\
                        p0=(2.,np.mean(dark_Resp_data)),maxfev=1500 )

    for iCW in range(half_GPP_CW,len(dark_Resp_data)-half_GPP_CW):
        # fit calc window data using least squares fit routine courtesy of scipy
        popt,pcov = curve_fit( Flambda,                                            \
                                dark_Tair_data[iCW-half_GPP_CW:iCW+half_GPP_CW],   \
                                dark_Resp_data[iCW-half_GPP_CW:iCW+half_GPP_CW],   \
                                p0=(2.,dark_Resp_data[iCW]),maxfev=1500 )
        Resp_data[iCW] = Flambda
    DarkResp  = MT.Q10_func(Tair_data,popt[0],popt[1],popt[2])
    
    
    for iTRES in range(nTRES):
        for attr in finf.ncattrs():
            foutfiles[iTRES].setncattr( str(attr), finf.getncattr(str(attr)) )
        foutfiles[iTRES].close()

    finf.close()

test=im.median_filter(test,5)

i=20000
length=300

plt.plot(time_nat_obj[i:i+length],NEE_data[i:i+length],label='NEE',ls=':')
plt.plot(time_nat_obj[i:i+length],test[i:i+length],label='NEE_filtered')

plt.plot(time_nat_obj[i:i+length],SWdown_data[i:i+length]/100.,label='SW_down')
plt.plot(time_nat_obj[i:i+length],Tair_data[i:i+length]-273.15,label='Tair')
plt.plot(time_nat_obj[i:i+length],np.zeros_like(NEE_data.data[i:i+length]),color='k')

plt.plot(time_nat_obj[i:i+length],dark_NEE_data[i:i+length],label='dark_NEE_filtered',lw=2)
plt.ylim(-5,25)
plt.legend()

plt.show()




plt.scatter(dark_Tair_data,dark_NEE_data)
plt.show()
