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
import netcdftime as ncdt
import matplotlib.pyplot as plot
import plot_tools as PT
#
###################################################################################################################
# Define class
###################################################################################################################
class plotsoilmaps:
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
        # optional
        parser.add_argument('--infileUK',type=str,help='Input file UK',required=False, \
                             default='/prj/ALANIS/deposition/EMEP4UK_Base_emep_4.3_2001_day.nc')
        parser.add_argument('--infileEU',type=str,help='Input file EU',required=False, \
                             default='/prj/ALANIS/deposition/EMEP4UK_Base_EUROPE_emep_4.3_2001_day.nc')
        parser.add_argument('--outdir',type=str,help='output directory' ,required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/plots/Base_Output/' )
        #
        # positional
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.infileUK, args.infileEU, args.outdir
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
    psm=plotsoilmaps()
    infileUK,infileEU,outdir =psm.parse_input()
    #
    #infileUK = '/prj/ALANIS/deposition/EMEP4UK_Base_emep_4.3_2001_day.nc'
    #outdir   = '/users/eow/edwcom/EMEP/EMEP4UK/plots/Base_Output/'
    #
    varnames=['WDEP_OXN','WDEP_RDN','WDEP_SO2','WDEP_HNO3',\
              'Emis_mgm2_BioNatC5H8','Emis_mgm2_BioNatNO','Emis_mgm2_voc','Emis_mgm2_BioNatAPINENE',\
              'DDEP_SOX_m2Grid','DDEP_OXN_m2Grid','DDEP_RDN_m2Grid']
    NLEVELSs   = [ 11 for counter in range(len(varnames))]
    cbars      = ['RdYlBu_r' for counter in range(len(varnames))]
    #data_ranges = [  [0.,1000.], \
    #                 [0.,500.], \
    #                 [0.0,0.05], \
    #                 [0.0,0.05], \
    #                 [-10.0,10.0], \
    #                 [-10.0,10.0], \
    #                 [90e3,105e3], \
    #                 [0.,0.01], \
    #                 [270.,310.]  ]
    #
    #
    # Open infile1 and extract data
    print 'Reading '+infileUK+':'
    inf1=nc.Dataset(infileUK,'r')
    time1=ncdt.num2date(inf1.variables['time'][:],          \
                        units=inf1.variables['time'].units, \
                        calendar='standard')
    lats1=inf1.variables['lat'][:]
    lons1=inf1.variables['lon'][:]
    data1={}
    # extract date_index from data 
    for var in varnames:
        print 'Reading '+var
        data1[var]={}
        data1[var]['longname']=inf1.variables[var].long_name
        data1[var]['units']=inf1.variables[var].units
        #
        tempdata=inf1.variables[var][:]
        tempdims=inf1.variables[var].shape
        # Calculate maps of mean, stddev, max and min of data for the year (first axis)
        data1[var]['map_mean'] = np.mean(tempdata,axis=0)
        data1[var]['map_std']  = np.std(tempdata,axis=0)
        data1[var]['map_max']  = np.amax(tempdata,axis=0)
        data1[var]['map_min']  = np.amin(tempdata,axis=0)
        #
        # Calculate time series of mean, stddev, min and max of data
        # Requires flattening the spatial dimensions
        data1[var]['TS_mean']  = np.mean(tempdata.reshape(tempdims[0],tempdims[1]*tempdims[2]),axis=1)
        data1[var]['TS_std']   = np.std(tempdata.reshape(tempdims[0],tempdims[1]*tempdims[2]),axis=1)
        data1[var]['TS_max']   = np.amax(tempdata.reshape(tempdims[0],tempdims[1]*tempdims[2]),axis=1)
        data1[var]['TS_min']   = np.amin(tempdata.reshape(tempdims[0],tempdims[1]*tempdims[2]),axis=1)
        #
        data1[var]['NLEVELS']  = NLEVELSs[varnames.index(var)]
        data1[var]['cbar']     = cbars[varnames.index(var)]
        
    inf1.close()
    #
    #
    lonrange=[-13,10.5]
    latrange=[51.5,57.]
    #
    #
    for var in varnames:
        print 'plotting maps for '+var
        data_range  = [ 0 ,round2SignifFigs(np.amax(data1[var]['map_mean']),2) ]  
        cbar_title = str(var+'  ('+data1[var]["units"]+')')
        cbar_title=cbar_title.replace('_',' ')
        
        plot_title = str(data1[var]['longname'])
        plot_title = plot_title.replace('_',' ')
        
        print data_range
        
        # Plot mean map
        plotname=outdir+''+str(var)+'_map_mean.png'
        PT.plot_map(data1[var]['map_mean'],lons1,lats1,                         \
                    DATA_RANGE=data_range,                                      \
                    NLEVELS=data1[var]['NLEVELS'], MPL_CBAR=data1[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' mean', CBAR_LABEL=cbar_title,       \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='h',                       \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)

        # Plot min map
        plotname=outdir+''+str(var)+'_map_min.png'
        PT.plot_map(data1[var]['map_min'],lons1,lats1,                         \
                    DATA_RANGE=data_range,                                      \
                    NLEVELS=data1[var]['NLEVELS'], MPL_CBAR=data1[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' min', CBAR_LABEL=cbar_title,       \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='h',                       \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)

        # Plot max map
        plotname=outdir+''+str(var)+'_map_max.png'
        PT.plot_map(data1[var]['map_max'],lons1,lats1,                         \
                    DATA_RANGE=data_range,                                      \
                    NLEVELS=data1[var]['NLEVELS'], MPL_CBAR=data1[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' max', CBAR_LABEL=cbar_title,       \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='h',                       \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
        
        
        # Plot time-series of all data
        # Mean - solid Red Line
        # Mean +/- std - dotted red lines
        # Max/Min - Solid black lines.
        #
        # Put Data into list
        PLOT_DATA  =[ data1[var]['TS_mean'] ] #,                      \
#                      data1[var]['TS_mean']-data1[var]['TS_std'], \
 #                     data1[var]['TS_mean']+data1[var]['TS_std'], \
  #                    data1[var]['TS_min'],                       \
   #                   data1[var]['TS_max']                        ]
        colours    = ['red']#,'red','red','black','black']
        linestyles = ['-']#,'--','--','-','-']
        plotname   = outdir+''+str(var)+'_timeseries.png'
        PT.plot_timeseries(PLOT_DATA, time1, \
                           PLOT_TITLE=plot_title, \
                           FONTSIZES=[15,15,18,15], \
                           COLOURS=colours, LINESTYLES=linestyles, \
                           FILE_PLOT=plotname,HEIGHT=4)

        

