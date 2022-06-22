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
import maths_tools.DateTimeTools as DTT
import maths_tools.FluxSiteTools as FST
import matplotlib.pyplot as plt

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
for site in site_list:
    if site in [site_list[itemp] for itemp in [10,17,21,28,39]]:
        continue
    print('Processing: '+site+'\n')
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
    SWdown_data = minf.variables['SWdown'][:].squeeze()
    Tair_data   = minf.variables['Tair'][:].squeeze()
    minf.close()
    
    #####################################################################
    # Follw same process as met data for fluxes, 
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
        # Not all params are available at all sites
        if (param in finf.variables):
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
                Npts_data   = np.array([ len(dat[dat.mask==False]) for dat in data_masked_Tindexed ])
                mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ])
                mean_data[mean_data!=mean_data] = fill_value
                median_data = np.array([ np.median(dat) for dat in data_masked_Tindexed ])
                median_data[median_data!=median_data] = fill_value
                stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
                stddev_data[stddev_data!=stddev_data] = fill_value
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
                foutvar[:]=median_data
                # StdDev
                foutvar=foutf.createVariable(param+'.stddev','float32',('time','land'),\
                                             fill_value=fill_value                     )
                foutvar.units=finf.variables[param].units
                foutvar.long_name=finf.variables[param].long_name+' - Standard Deviation'
                foutvar[:]=stddev_data
                # Npts
                foutvar=foutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                             fill_value=fill_value                     )
                foutvar.units=''
                foutvar.long_name=finf.variables[param].long_name+' - Number of Averaged time-steps'
                foutvar[:]=Npts_data


    # Calculate Respiration (TER) and  GPP
    # Respiration is calculated by fitting the Night-time data (DarkResp)
    #    to a Q10 function
    # Night-time data is defined as data which is not within 2 timesteps
    #    of registered SWdown
     
    #GPP_data, TER_data, GPPFitParams = \
    #                FST.GPP_from_NEE_SW_T(NEE_data,SWdown_data,Tair_data,\
    #                                      spike_filter='gaussian',\
    #                                      SW_smooth_int=3, \
    #                                      ReturnResp=True, \
    #                                      ReturnFitParams=True )  
    #            #                         ReturnNEEfiltered=True )
    #GPP_data, TER_data, GPPFitParams = \
    #                FST.GPP_from_NEE_SW_T(NEE_data,SWdown_data,Tair_data,\
    #                                      spike_filter='stddev_mask', FIT_maxfev=2000, \
    #                                      STDmask_window=2,STDmask_limit=np.std(NEE_data.data),   \
    #                                      ReturnResp=True, \
    #                                      ReturnFitParams=True )  
    GPP_data, TER_data, GPPFitParams = \
                    FST.GPP_from_NEE_SW_T(NEE_data,SWdown_data,Tair_data,\
                                          spike_filter='median', FIT_maxfev=2000, \
                                          ReturnResp=True, \
                                          ReturnFitParams=True )  
    
    plt.plot(time_nat_obj,NEE_data,label='NEE')
    plt.plot(time_nat_obj,GPP_data,label='GPP')
    #plt.plot(time_nat_obj,NEE_filtered,label='NEE-filtered')
    plt.plot(time_nat_obj,TER_data,label='TER')
    plt.plot(time_nat_obj,np.zeros_like(NEE_data.data),color='k')
    plt.legend(ncol=3)
    plt.savefig(out_flux_file_tag+'.png')
    plt.close()


    # now create and store the GPP and TER
    param='GPP'
    foutvars = []
    for foutf in foutfiles:
        foutvar=foutf.createVariable(param,'float32',\
                                     ('time','land'),\
                                     fill_value=fill_value)
        for attr in finf.variables[param].ncattrs():
            foutvar.setncattr(str(attr),finf.variables[param].getncattr(str(attr)))
        
        foutvars.append(foutvar)
    data_masked=GPP_data
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
        Npts_data   = np.array([ len(dat[dat.mask==False]) for dat in data_masked_Tindexed ])
        mean_data   = np.array([ np.mean(dat[dat.mask==False]) for dat in data_masked_Tindexed ])
        mean_data[mean_data!=mean_data] = fill_value
        median_data = np.array([ np.median(dat[dat.mask==False]) for dat in data_masked_Tindexed ])
        median_data[median_data!=median_data] = fill_value
        stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
        stddev_data[stddev_data!=stddev_data] = fill_value
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
        foutvar[:]=median_data
        # StdDev
        foutvar=foutf.createVariable(param+'.stddev','float32',('time','land'),\
                                     fill_value=fill_value                     )
        foutvar.units=finf.variables[param].units
        foutvar.long_name=finf.variables[param].long_name+' - Standard Deviation'
        foutvar[:]=stddev_data
        # Npts
        foutvar=foutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                     fill_value=fill_value                     )
        foutvar.units=''
        foutvar.long_name=finf.variables[param].long_name+' - Number of Averaged time-steps'
        foutvar[:]=Npts_data
        
    
    # now create and store the TER
    param='TER'
    plongname = 'Total Ecosystem Respiration'
    punits    = u'W/m^2'
    foutvars = []
    for foutf in foutfiles:
        foutvar=foutf.createVariable(param,'float32',\
                                     ('time','land'),\
                                     fill_value=fill_value)
        foutvar.units=punits
        foutvar.long_name=plongname 
        
        foutvars.append(foutvar)
    
    data_masked=TER_data
    
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
        Npts_data   = np.array([ len(dat[dat.mask==False]) for dat in data_masked_Tindexed ])
        mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ])
        mean_data[mean_data!=mean_data] = fill_value
        median_data = np.array([ np.median(dat) for dat in data_masked_Tindexed ])
        median_data[median_data!=median_data] = fill_value
        stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
        stddev_data[stddev_data!=stddev_data] = fill_value
        # stderr_data = stddev_data/Npts_data
        
        # Store the mean into the pre created variable 
        foutvars[iTRES][:] = mean_data
        foutvars[iTRES].long_name+=' - Mean'
        
        # Create variables for remaining parameters
        # Median
        foutvar=foutf.createVariable(param+'.median','float32',('time','land'),\
                                     fill_value=fill_value                     )
        foutvar.units=punits
        foutvar.long_name=plongname+' - Median'
        foutvar[:]=median_data
        # StdDev
        foutvar=foutf.createVariable(param+'.stddev','float32',('time','land'),\
                                     fill_value=fill_value                     )
        foutvar.units=punits
        foutvar.long_name=plongname+' - Standard Deviation'
        foutvar[:]=stddev_data
        # Npts
        foutvar=foutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                     fill_value=fill_value                     )
        foutvar.units=''
        foutvar.long_name=plongname+' - Number of Averaged time-steps'
        foutvar[:]=Npts_data
        

    for iTRES in range(nTRES):
        for attr in finf.ncattrs():
            foutfiles[iTRES].setncattr( str(attr), finf.getncattr(str(attr)) )
        foutfiles[iTRES].note='GPP and TER calculated by Edward Comyn-Platt, edwcom@ceh.ac.uk'
        foutfiles[iTRES].close()

    finf.close()

