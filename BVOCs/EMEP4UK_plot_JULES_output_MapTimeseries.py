#!/usr/bin/python
#
# Python module to plot the hwsd dat on EMEP grid
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import argparse
import netCDF4 as nc
import netcdftime as nctime
import matplotlib.pyplot as plt
import plot_tools as PT
import EMEP4UK_tools as ET
#
###################################################################################################################
# Define class
###################################################################################################################
class plotemep4uk:
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
        # optional
        parser.add_argument('--data_dir',type=str,help='Data Directory',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/')
        parser.add_argument('--outdir',type=str,help='output directory' ,required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_Output/' )
        parser.add_argument('--run_id',type=str,help='JULES run ID',required=False, \
                             default='EMEP4UK_BC')
        parser.add_argument('--profile_id',type=str,help='JULES profile ID',required=False, \
                             default='hourly')
        parser.add_argument('--file_period',type=str,help='JULES file period',required=False, \
                             default='monthly')
        parser.add_argument('--start_year',type=int,help='Start Year of analysis period' ,required=False, \
                             default=2001)
        parser.add_argument('--end_year',type=int,help='End Year of analysis period' ,required=False, \
                             default=2014)
        parser.add_argument('--Mean_Period',type=str,help='Period over which to mean the data' ,required=False, \
                             default='monthly')
        parser.add_argument('--Map_Type',type=str,help='Which map to produce (Max, Mean, Total)' ,required=False, \
                            default='Max')
        parser.add_argument('--varnames',type=str,help='Variables to plot' ,required=False, \
                            default='isoprene_gb')

        #
        # positional
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.data_dir, args.outdir, \
                args.run_id, args.profile_id, args.file_period, \
                args.start_year, args.end_year, \
                args.Mean_Period, args.Map_Type, \
                args.varnames
#
###################################################################################################################
# Define functions
###################################################################################################################
# Round to given significant figures
def round2SignifFigs(vals,n):
    mags = 10.0**np.floor(np.log10(np.abs(vals)))  # order of mag's
    outvals = np.around(vals/mags,n-1)*mags             # round(val/omag)*omag
    try:
        outvals[np.where(np.isnan(vals))] = 0.0           # where order of mag = 0, set to zero
    except:
        if np.isnan(outvals):
            outvals=0.0
    #
    return outvals
#
###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    pemep=plotemep4uk()
    data_dir, outdir, \
     run_id, profile_id, file_period, \
     start_year, end_year, \
     Mean_Period, Map_Type,  \
     varnames,               \
    = pemep.parse_input()

    #
    #data_dir   = '/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/'
    #run_id     = 'EMEP4UK_BC'
    #profile_id = 'hourly'
    #varnames='isoprene_gb,terpene_gb'
    #outdir     = '/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_Output/032016/'
    #start_year = 2001
    #end_year   = 2001
    #file_period= 'monthly'
    #Mean_Period='monthly'
    #Map_Type='Max'
    
    EMEP4UK_2D_indexfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc'
    EMEP4UK_2D_latlonfile='/users/eow/edwcom/EMEP/EMEP4UK/topidx_emep4uk_grid.nc'
    EMEP4UK_2D_latname='cen_lat'
    EMEP4UK_2D_lonname='cen_lon'
    
    varnames=varnames.split(',')
    nYEARs=end_year-start_year+1

    month_strings = [ '01','02','03','04','05','06','07','08','09','10','11','12' ]
    month_days    = [  31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31  ]
    #
    NLEVELSs   = [ 11 for counter in range(len(varnames))]
    data_ranges  = [ [0.,2], [0.,4],  ]
    conv_factors = [ 3600*24*365*1000., 3600*24*365*1000., ]
    data_units   = [ '$g.m^{-2}.yr^{-1}$','$g.m^{-2}.yr^{-1}$', ]
     
    #
    #
    # UK data
    #
    # Fetch 2D lat lon grids:
    inf=nc.Dataset(EMEP4UK_2D_latlonfile,'r')
    lats_2D=inf.variables[EMEP4UK_2D_latname][:]
    lons_2D=inf.variables[EMEP4UK_2D_lonname][:]
    inf.close()
    #
    lats_1D=None
    lons_1D=None
    index_2D=None
    save_grid_file=True
    data = {}
    if file_period=='monthly':
        for year in range(start_year,end_year+1):
            for mnth in range(12):
                
                month = month_strings[mnth]
                print year, month
                loop_file=data_dir+run_id+'/'+ \
                            run_id+'.'+profile_id+'.'+str(year)+month+'.nc'
                print loop_file
                inf=nc.Dataset(loop_file,'r')
                if (lats_1D==None):
                    lats_1D=inf.variables['latitude'][:].squeeze()
                if (lons_1D==None):
                    lons_1D=inf.variables['longitude'][:].squeeze()
                if index_2D==None:
                    index_2D=np.zeros_like(lats_2D,dtype='int')-999
                    for i,lat,lon in zip(range(len(lats_1D)),lats_1D,lons_1D):
                        index_2D[(lats_2D==lat)&(lons_2D==lon)] = int(i)
                        index_2D=np.ma.masked_equal(index_2D,-999)
                
                if save_grid_file:
                    groutf=nc.Dataset(EMEP4UK_2D_indexfile,'w')
                    groutf.createDimension('i_EMEP',index_2D.shape[0])
                    groutf.createDimension('j_EMEP',index_2D.shape[1])
                    groutvar=groutf.createVariable('lats','float32',('i_EMEP','j_EMEP'))
                    groutvar[:]=lats_2D
                    groutvar=groutf.createVariable('lons','float32',('i_EMEP','j_EMEP'))
                    groutvar[:]=lons_2D
                    groutvar=groutf.createVariable('land_index','int',('i_EMEP','j_EMEP'),fill_value=-1)
                    groutvar[:]=index_2D
                    groutf.note='index to convert JULES output onto EMEP4UK 2D grid'
                    groutf.close()
                    save_gird_file=False
    

                temptime=inf.variables['time'][:]
                if Mean_Period=='monthly':
                    temptime=np.array([np.mean(temptime,axis=0)])
                elif Mean_Period=='daily':
                    days=month_days[mnth]
                    # add extra day to february in leap year
                    if ( (year==2004)|(year==2008)|(year==2012)|(year==2016) ) & (month=='02'):
                        days+=1
                    # average for each day
                    old_dims = list(temptime.shape)
                    new_dims = [days,old_dims[0]/days]
                    temptime=np.mean(temptime.reshape(new_dims),axis=1)

                print temptime.shape 
                if 'time' not in data.viewkeys():
                    print 'a'
                    data['time']= nctime.num2date( np.array(temptime),  \
                                                   units=inf.variables['time'].units, \
                                                   calendar='standard'                ) 
                else:
                    print 'b'
                    data['time']= np.append(data['time'],nctime.num2date( temptime,\
                                                            units=inf.variables['time'].units, \
                                                                          calendar='standard'),axis=0  )  


                for var in varnames:
                    tempdata=inf.variables[var][:].squeeze()
                    if Mean_Period=='monthly':
                        tempdata=np.mean(tempdata,axis=0)
                        tempdata=tempdata.reshape([1]+list(tempdata.shape))
                    elif Mean_Period=='daily':
                        old_dims=list(tempdata.shape)
                        new_dims=[days,old_dims[0]/days]+old_dims[1:]
                        tempdata=np.mean(tempdata.reshape(new_dims),axis=1)

                    if var not in data.viewkeys():
                        data[var]={}
                        data[var]['longname']=inf.variables[var].long_name
                        data[var]['units']=inf.variables[var].units
                        data[var]['data']=tempdata
                    else:
                        data[var]['data']=np.append( data[var]['data'], tempdata , axis=0)
                inf.close()
    
    #data['time']=np.array(data['time'])
    #for var in varnames:
    #    data[var]['data']=np.array(data[var]['data'])
    
    lonrange=[-13,11]
    latrange=[51.5,56.5]
    #
    for ivar in range(len(varnames)):
        var = varnames[ivar]
        print 'plotting maps for '+var
        data_range  = data_ranges[ivar] 
        cbar_title = str(var+'  ('+data_units[ivar]+')')
        cbar_title=cbar_title.replace('_',' ')
        #
        plot_title = str(data[var]['longname'])
        plot_title = plot_title.replace('_',' ')
        #
        print data_range
        #
        # Plot max map
        if Map_Type=='Max':
            plot_data=np.ma.masked_array( np.max(data[var]['data'],axis=0)[index_2D], \
                                               mask=index_2D.mask )    \
                                               *conv_factors[ivar]
        elif Map_Type=='Mean':
            plot_data=np.ma.masked_array( np.mean(data[var]['data'],axis=0)[index_2D], \
                                               mask=index_2D.mask )    \
                                               *conv_factors[ivar]
        elif Map_Type=='Total':
            plot_data=np.ma.masked_array( np.mean(data[var]['data'],axis=0)[index_2D]*nYEARs, \
                                               mask=index_2D.mask )    \
                                               *conv_factors[ivar]


        plotname=outdir+'EMEP4UK_'+str(var)+'_map_'+Map_Type+'_'+Mean_Period+'means_'+\
                    str(start_year)+'_'+str(end_year)+'.png'
        PT.plot_map(plot_data,lons_2D,lats_2D,                                  \
                    DATA_RANGE=data_range,                                      \
                    COLOURS=['beige','greenyellow','darkgreen'],                \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,                        \
                    CBAR_ORIENTATION='vertical',                                \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' '+Map_Type, CBAR_LABEL=cbar_title,  \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='i',                       \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
        
        print 'plot produced and saved to: '+plotname
        
        #
        #
        # Plot time-series of all data
        # Put Data into list
        mean_TS = np.mean(data[var]['data'],axis=1)*conv_factors[ivar]
        min_TS  = np.min(data[var]['data'],axis=1)*conv_factors[ivar]
        max_TS  = np.max(data[var]['data'],axis=1)*conv_factors[ivar]
        std_TS  = np.std(data[var]['data'],axis=1)*conv_factors[ivar]

        PLOT_DATA  = [ mean_TS, mean_TS+std_TS ]  # max_TS
        TIME_DATA  = [ data['time']  for i in range(len(PLOT_DATA)) ]
        legend_names= [ 'mean','mean+std']
        colours    = ['black','black']
        linestyles = ['-','--']
        plotname   = outdir+str(var)+'_timeseries_'+Mean_Period+'means_'+\
                      str(start_year)+'_'+str(end_year)+'.png'
        PT.plot_timeseries(PLOT_DATA, TIME_DATA, \
                           PLOT_TITLE=plot_title, \
                           Y_LABEL=cbar_title,    \
                           FONTSIZES=[15,15,18,15], \
                           COLOURS=colours, LINESTYLES=linestyles, \
                           FILE_PLOT=plotname,HEIGHT=4, \
                           LEGEND='upper right', LEGEND_DATANAMES=legend_names)


        

