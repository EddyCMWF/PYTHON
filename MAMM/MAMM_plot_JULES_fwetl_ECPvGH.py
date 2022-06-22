#!/usr/bin/python
#
# Python module to plot JULES output on the MAMM grid
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
import netcdftime as ncdt
import matplotlib.pyplot as plot
import plot_tools as PT
#
###################################################################################################################
# Define class
###################################################################################################################
class plotMAMM:
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
        # optional
        parser.add_argument('--variable',type=str,help='netCDF variable to compare',required=False, \
                             default='fwetl')
        parser.add_argument('--outdir',type=str,help='output directory' ,required=False, \
                             default='/users/eow/edwcom/MAMM/plots/' )
        parser.add_argument('--met_data',type=str,help='select which met data', required=False, \
                             default='cruncep' )
        #
        # positional
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.variable, args.outdir, args.met_data
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
    pmamm=plotMAMM()
    variable,outdir,met_data =pmamm.parse_input()
    #
    #variable='fwetl'
    #variable='fch4_wetl'
    #outdir='/users/eow/edwcom/MAMM/plots/'
    #
    #met_data='cruncep'
    if (met_data=='cruncep'):
        GH_dir='/prj/ALANIS/jules_data/jules_output_MAMM/MAMM_Scandinavia_job001_NCEP_CRU_v4_wd_cirrus/z_EXTRACT/'
        #
        infile_1='/users/eow/edwcom/MAMM/JULES_output/JULESv4.2_MPI_MAMM_DougDiag.monthly.nc'
        if (variable=='fwetl'):
            infile_2=GH_dir+'JULES_m046_monthly_1980_2012_wetland_fraction_scan_20140729.nc'
        elif (variable=='fch4_wetl'):
            infile_2=GH_dir+'JULES_m046_monthly_1980_2012_wetland_emissions_scan_20140729.nc'
            #
        outdir=outdir+'cruncep/'
    elif (met_data=='wfdei'):
        GH_dir='/prj/ALANIS/jules_data/jules_output_MAMM/MAMM_Scandinavia_job001_WFD_EI_cirrus/z_EXTRACT/'
        #
        infile_1='/users/eow/edwcom/MAMM/JULES_output/MPI_MAMM_v42_DougDiag_WFDEI.monthly.nc'
        if (variable=='fwetl'):
            infile_2=GH_dir+'JULES_m046_monthly_1980_2012_wetland_fraction_scan_20140324.nc'
        elif (variable=='fch4_wetl'):
            infile_2=GH_dir+'JULES_m046_monthly_1980_2012_wetland_emissions_scan_20140324.nc'
            #
        outdir=outdir+'wfdei/'
    #
    #os.mkdir(outdir)
    #
    #
    # ECP - JASMIN Data
    #
    inf       = nc.Dataset(infile_1)
    lats1     = inf.variables['latitude'][:].squeeze()
    lons1     = inf.variables['longitude'][:].squeeze()
    time1     = ncdt.num2date(inf.variables['time'][:],          \
                              units=inf.variables['time'].units, \
                              calendar='standard')
    data1_1D  = inf.variables[variable][:].squeeze()
    # get meta data
    long_name = inf.variables[variable].long_name
    units     = inf.variables[variable].units
    inf.close()
    #
    # GH  - cirrus data
    #
    inf      = nc.Dataset(infile_2)
    lats2    = inf.variables['latitude'][:]
    lons2    = inf.variables['longitude'][:]
    time2     = ncdt.num2date(inf.variables['time'][:],          \
                              units=inf.variables['time'].units, \
                              calendar='standard')
    data2_2D = inf.variables[variable][:]
    inf.close()
    #
    # GH data on 2D grid, reconstruct 2D lat lon grids then map ECP data onto this
    #
    lons_2D, lats_2D = np.meshgrid(lons2,lats2)
    #
    index = [ np.argmin( ((lats_2D-lat)**2) + ((lons_2D-lon)**2) ) \
              for lat,lon in zip(list(lats1),list(lons1)) ]
    #
    data1_2D=np.zeros(data2_2D.shape)+data2_2D.fill_value
    for cnt in range((data1_2D.shape)[0]):
        data1_2D[cnt,:,:].flat[index]=data1_1D[cnt,:]
    data1_2D = np.ma.masked_equal(data1_2D,data2_2D.fill_value)
    
    if (variable=='fwetl'):
        data_range=[0,0.3]
        diff_range=[-0.01,0.01]
        max_diff_range=[-0.3,0.3]
        std_range=[0.,0.05]
    elif (variable=='fch4_wetl'):
        data_range=[0,0.1]
        diff_range=[-0.1,0.1]
        max_diff_range=[-0.1,0.1]
        std_range=[0.,0.1]
    else:
        data_range=[0,np.max(data1)]
    
    units=units.replace('_',' ')
    units=units.replace('^','')
    
    
    # Differences
    diff_data = data1_2D-data2_2D
    
    # data names for plots
    data1_name='v4.2-J'
    data2_name='v3.4.1-c'
    diff_name ='v4.2J-minus-v3.4.1c'
    
    plot_maps=True
    if (plot_maps==True):
        # plot mean of data 1 map
        plotname=outdir+str(variable)+data1_name+'_map_mean.png'
        PLOT_data=np.mean(data1_2D,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=data_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+data1_name+' mean', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )
        #
        # plot max of data 1 map
        plotname  = outdir+str(variable)+data1_name+'_map_max.png'
        PLOT_data = np.max(data1_2D,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=data_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+data1_name+' max', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )
        #
        #
        # plot mean of data 2 map
        plotname=outdir+str(variable)+data2_name+'_map_mean.png'
        PLOT_data=np.mean(data2_2D,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=data_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+data2_name+' mean', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )
        #
        # plot max of data 2 map
        plotname  = outdir+str(variable)+data2_name+'_map_max.png'
        PLOT_data = np.max(data2_2D,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=data_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+data2_name+' max', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )
        #
        #
        # plot mean difference between 2 datasets
        plotname  = outdir+str(variable)+diff_name+'_map_mean.png'
        PLOT_data = np.mean(diff_data,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=diff_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+diff_name+' mean', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )
        #
        # plot standard deviation of difference between 2 datasets
        plotname  = outdir+str(variable)+diff_name+'_map_stddev.png'
        PLOT_data = np.std(diff_data,axis=0)
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=std_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+' '+diff_name+' std', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )    
        #
        # plot standard deviation of difference between 2 datasets
        plotname  = outdir+str(variable)+diff_name+'_map_max_diff.png'
        max_data  = np.max(diff_data,axis=0)
        min_data  = np.min(diff_data,axis=0)
        
        PLOT_data = max_data
        PLOT_data[np.where((0.-min_data)>max_data)]=min_data[np.where((0.-min_data)>max_data)]
        
        PT.plot_map(PLOT_data,lons_2D,lats_2D,                      \
                    DATA_RANGE=max_diff_range,                          \
                    NLEVELS=11, MPL_CBAR='RdYlBu_r',                \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',   \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=long_name+diff_name+' max diff', CBAR_LABEL=units, \
                    iDISPLAY='N', FILE_PLOT=plotname,               \
                    LATDEL=5., LONDEL=5., RESOLUTION='h',           \
                    PROJECTION='cyl' )    

    # Plot time-series of all data
    # Put Data into list
    data1_mean_TS = np.mean(data1_1D,axis=1)
    data1_std_TS  = np.std(data1_1D,axis=1)
    data1_max_TS  = np.max(data1_1D,axis=1)

    data2_1D = data2_2D.reshape(data2_2D.shape[0],data2_2D.shape[1]*data2_2D.shape[2])
    data2_mean_TS = np.mean(data2_1D,axis=1)
    data2_std_TS  = np.std(data2_1D,axis=1)
    data2_max_TS  = np.max(data2_1D,axis=1)
   
    PLOT_DATA  = [ data1_mean_TS, data1_max_TS, data1_mean_TS+data1_std_TS, data1_mean_TS-data1_std_TS, \
                   data2_mean_TS, data2_max_TS, data2_mean_TS+data2_std_TS, data2_mean_TS-data2_std_TS  ]
    TIME_DATA  = [ time1, time1, time1, time1, \
                   time2, time2, time2, time2  ]
    colours    = ['red','red','red','red','black','black','black','black']
    linestyles = ['-','-',':',':','-','-',':',':']
    plotname   = outdir+''+variable+'_timeseries.png'
    PT.plot_timeseries(PLOT_DATA, TIME_DATA, \
                       PLOT_TITLE=long_name, \
                       FONTSIZES=[15,15,18,15], \
                       COLOURS=colours, LINESTYLES=linestyles, \
                       FILE_PLOT=plotname,HEIGHT=4)
    
    diff_data_1D = diff_data.reshape(diff_data.shape[0],diff_data.shape[1]*diff_data.shape[2])
    diff_mean_TS = np.mean(diff_data_1D,axis=1)
    diff_std_TS  = np.std(diff_data_1D,axis=1)
    diff_max_TS  = np.max(diff_data_1D,axis=1)
    diff_min_TS  = np.min(diff_data_1D,axis=1)
    
    PLOT_DATA  = [ diff_mean_TS.data, diff_max_TS.data, diff_min_TS.data, 
                   (diff_mean_TS+diff_std_TS).data, (diff_mean_TS-diff_std_TS).data ]
    TIME_DATA  = [ time2, time2, time2, time2, time2 ]    
    print len(time1)
    print len(time2)
    for dat in PLOT_DATA:
        print len(dat)
    colours    = ['red','black','black','black','black']
    linestyles = ['-','-','-',':',':']
    plotname   = outdir+''+variable+'_'+diff_name+'_timeseries.png'
    PT.plot_timeseries(PLOT_DATA, TIME_DATA, \
                       PLOT_TITLE=long_name+' '+diff_name, \
                       FONTSIZES=[15,15,18,15], \
                       COLOURS=colours, LINESTYLES=linestyles, \
                       FILE_PLOT=plotname,HEIGHT=4)
    
    

