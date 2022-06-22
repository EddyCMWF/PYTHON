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
        data_masked=data
        
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

            #####################################################################
            # Loop round the time indexes and calculate the mean, median, stddev
            Npts_data   = np.array([ len(dat) for dat in data_masked_Tindexed ])
            mean_data   = np.array([ np.mean(dat) for dat in data_masked_Tindexed ])
            median_data = np.array([ np.median(dat) for dat in data_masked_Tindexed ])
            stddev_data = np.array([ np.std(dat) for dat in data_masked_Tindexed ])
            # Filter out bad data, typically the last tstep
            Npts_data[Npts_data!=Npts_data] = fill_value
            mean_data[mean_data!=mean_data] = fill_value
            median_data[median_data!=median_data] = fill_value
            stddev_data[stddev_data!=stddev_data] = fill_value
            
            # Store the mean into the pre created variable 
            moutvars[iTRES][:] = mean_data
            moutvars[iTRES].long_name+=' - Mean'

            # Create variables for remaining parameters
            # Median
            moutvar=moutf.createVariable(param+'.median','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=minf.variables[param].units
            moutvar.long_name=minf.variables[param].long_name+' - Median'
            moutvar[:]=median_data
            # StdDev
            moutvar=moutf.createVariable(param+'.stddev','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=minf.variables[param].units
            moutvar.long_name=minf.variables[param].long_name+' - Standard Deviation'
            moutvar[:]=stddev_data
            # Npts
            moutvar=moutf.createVariable(param+'.Ntsteps','float32',('time','land'),\
                                          fill_value=fill_value                     )
            moutvar.units=''
            moutvar.long_name=minf.variables[param].long_name+' - Number of Averaged time-steps'
            moutvar[:]=Npts_data
        
        del data_masked
        del data_masked_Tindexed
        del mean_data

    for iTRES in range(nTRES):
        for attr in minf.ncattrs():
            moutfiles[iTRES].setncattr( str(attr), minf.getncattr(str(attr)) )
        moutfiles[iTRES].close()

    minf.close()
    
    