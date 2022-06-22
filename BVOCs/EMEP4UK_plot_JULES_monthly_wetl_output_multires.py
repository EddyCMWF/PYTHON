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
import EMEP4UK_tools as ET
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
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_MPI_DougDiag.monthly_wetl.nc')
        parser.add_argument('--infileEU',type=str,help='Input file EU',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_EUROPE_MPI_DougDiag.monthly_wetl.nc')
        parser.add_argument('--outdir',type=str,help='output directory' ,required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_Output/' )
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
    infileUK = '/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_MPI_DougDiag.monthly_wetl.nc'
    infileEU = '/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_EUROPE_MPI_DougDiag.monthly_wetl.nc'
    outdir   = '/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_Output/'
    EMEP4UK_2D_indexfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc'
    EMEP4UK_2D_latlonfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Landuse.nc'
    EMEP4UK_2D_latname='lat'
    EMEP4UK_2D_lonname='lon'
    EMEP4UK_EUROPE_2D_indexfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_JULES_output_index.nc'
    EMEP4UK_EUROPE_2D_latlonfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_Landuse.nc'
    EMEP4UK_EUROPE_2D_latname='lat'
    EMEP4UK_EUROPE_2D_lonname='lon'
    #
    varnames=['fwetl','fch4_wetl','fwetl_unf','fwetl_all']
    #
    NLEVELSs   = [ 11 for counter in range(len(varnames))]
    cbars      = ['YlGnBu' , 'YlGn' , 'YlGnBu' , 'YlGnBu' ]
    data_ranges_1= [ [0.,0.3], [0.,0.1], [0.,0.3], [0.,0.3] ]
    data_ranges_2= [ [0.,0.3], [0.,0.2], [0.,0.3], [0.,0.3] ]
    #
    #
    # UK data
    #
    # Fetch 2D lat lon grids:
    inf=nc.Dataset(EMEP4UK_2D_latlonfile,'r')
    lats1_2D=inf.variables[EMEP4UK_2D_latname][:]
    lons1_2D=inf.variables[EMEP4UK_2D_lonname][:]
    inf.close()
    #
    # Open infile1 and extract data
    print 'Reading '+infileUK+':'
    inf1=nc.Dataset(infileUK,'r')
    time1=ncdt.num2date(inf1.variables['time'][:],          \
                        units=inf1.variables['time'].units, \
                        calendar='standard')
    lats1=inf1.variables['latitude'][:]
    lons1=inf1.variables['longitude'][:]
    #
    # Following section created the index file for initial run
    #stuff = ET.conv_EMEP4UK_2D_raw(lats1,lons1,lats1,                                              \
    #                            grid_file=EMEP4UK_2D_latlonfile,                                   \
    #                            grinf_latname=EMEP4UK_2D_latname,grinf_lonname=EMEP4UK_2D_lonname, \
    #                            missing_value=-9999.0,outfilename=EMEP4UK_2D_indexfile             )
    #                    
    #
    data1={}
    # extract date_index from data 
    for var in varnames:
        data1[var]={}
        data1[var]['longname']=inf1.variables[var].long_name
        data1[var]['units']=inf1.variables[var].units
        #
        tempdata=inf1.variables[var][:].squeeze()   # Remove redundant spatial dimension
        tempdims=tempdata.shape
        #
        # Calculate maps of mean, stddev, max and min of data for the year (first axis)
        data1[var]['map_mean'] = ET.conv_EMEP4UK_2D(np.mean(tempdata,axis=0))
        data1[var]['map_std']  = ET.conv_EMEP4UK_2D(np.std(tempdata,axis=0))
        data1[var]['map_max']  = ET.conv_EMEP4UK_2D(np.amax(tempdata,axis=0))
        data1[var]['map_min']  = ET.conv_EMEP4UK_2D(np.amin(tempdata,axis=0))
        #
        # Calculate time series of mean, stddev, min and max of data
        # Requires flattening the spatial dimensions
        data1[var]['TS_mean']  = np.mean(tempdata,axis=1)
        data1[var]['TS_std']   = np.std(tempdata,axis=1)
        data1[var]['TS_max']   = np.amax(tempdata,axis=1)
        data1[var]['TS_min']   = np.amin(tempdata,axis=1)
        #
        data1[var]['NLEVELS']  = NLEVELSs[varnames.index(var)]
        data1[var]['cbar']     = cbars[varnames.index(var)]
        data1[var]['drange']   = data_ranges_1[varnames.index(var)]
    #    
    inf1.close()
    #
    # EUROPE data
    #
    # Fetch 2D lat lon grids:
    inf=nc.Dataset(EMEP4UK_EUROPE_2D_latlonfile,'r')
    lats2_2D=inf.variables[EMEP4UK_EUROPE_2D_latname][:]
    lons2_2D=inf.variables[EMEP4UK_EUROPE_2D_lonname][:]
    inf.close()
    #
    # Open infile2 and extract data
    print 'Reading '+infileEU+':'
    inf2=nc.Dataset(infileEU,'r')
    time2=ncdt.num2date(inf2.variables['time'][:],          \
                        units=inf2.variables['time'].units, \
                        calendar='standard')
    lats2=inf2.variables['latitude'][:]
    lons2=inf2.variables['longitude'][:]
    #
    # Following section created the index file for initial run
    # stuff = ET.conv_EMEP4UK_2D_raw(lats2,lons2,lats2,                                              \
    #                               grid_file=EMEP4UK_EUROPE_2D_latlonfile,                                   \
    #                               grinf_latname=EMEP4UK_EUROPE_2D_latname,grinf_lonname=EMEP4UK_EUROPE_2D_lonname, \
    #                               missing_value=-9999.0,outfilename=EMEP4UK_EUROPE_2D_indexfile             )
    #
    data2={}
    # extract date_index from data 
    for var in varnames:
        data2[var]={}
        data2[var]['longname']=inf2.variables[var].long_name
        data2[var]['units']=inf2.variables[var].units
        #
        tempdata=inf2.variables[var][:].squeeze()   # Remove redundant spatial dimension
        tempdims=tempdata.shape
        #
        # Calculate maps of mean, stddev, max and min of data for the year (first axis)
        data2[var]['map_mean'] = ET.conv_EMEP4UK_2D(np.mean(tempdata,axis=0),index_file=EMEP4UK_EUROPE_2D_indexfile)
        data2[var]['map_std']  = ET.conv_EMEP4UK_2D(np.std(tempdata,axis=0),index_file=EMEP4UK_EUROPE_2D_indexfile)
        data2[var]['map_max']  = ET.conv_EMEP4UK_2D(np.amax(tempdata,axis=0),index_file=EMEP4UK_EUROPE_2D_indexfile)
        data2[var]['map_min']  = ET.conv_EMEP4UK_2D(np.amin(tempdata,axis=0),index_file=EMEP4UK_EUROPE_2D_indexfile)
        #
        # Calculate time series of mean, stddev, min and max of data
        # Requires flattening the spatial dimensions
        data2[var]['TS_mean']  = np.mean(tempdata,axis=1)
        data2[var]['TS_std']   = np.std(tempdata,axis=1)
        data2[var]['TS_max']   = np.amax(tempdata,axis=1)
        data2[var]['TS_min']   = np.amin(tempdata,axis=1)
        #
        data2[var]['NLEVELS']  = NLEVELSs[varnames.index(var)]
        data2[var]['cbar']     = cbars[varnames.index(var)]
        data2[var]['drange']   = data_ranges_2[varnames.index(var)]
    #    
    inf2.close()
    #
    latrange = [40.,38.]
    lonrange = [-30.,60.]
    CU_lonrange=[-15,12]
    CU_latrange=[50.,56.]
    #
    #
    for var in varnames:
        print 'plotting maps for '+var
        data_range  = [ 0,round2SignifFigs(np.amax(data1[var]['map_mean']),2) ]  
        cbar_title = str(var+'  ('+data1[var]["units"]+')')
        cbar_title=cbar_title.replace('_',' ')
        #
        plot_title = str(data1[var]['longname'])
        plot_title = plot_title.replace('_',' ')
        #
        print data_range
        #
        # Plot Europe Data
        # Plot mean map
        plotname=outdir+'EMEP4UK_EUROPE_'+str(var)+'_map_mean.png'
        PT.plot_map(data2[var]['map_mean'],lons2_2D,lats2_2D,                   \
                    DATA_RANGE=data1[var]['drange'],                            \
                    NLEVELS=data2[var]['NLEVELS'], MPL_CBAR=data2[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' mean', CBAR_LABEL=cbar_title,       \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=15., LONDEL=15., RESOLUTION='h',                     \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange   )
        #
        # Plot max map
        plotname=outdir+'EMEP4UK_EUROPE_'+str(var)+'_map_max.png'
        PT.plot_map(data2[var]['map_max'],lons2_2D,lats2_2D,                    \
                    DATA_RANGE=data2[var]['drange'],                            \
                    NLEVELS=data2[var]['NLEVELS'], MPL_CBAR=data2[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' max', CBAR_LABEL=cbar_title,        \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=15., LONDEL=15., RESOLUTION='h',                     \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange   )
        #
        #
        # Plot mean map
        plotname=outdir+'EMEP4UK_'+str(var)+'_map_mean.png'
        PT.plot_map(data2[var]['map_mean'],lons2_2D,lats2_2D,                   \
                    DATA2=data1[var]['map_mean'],LONS2=lons1_2D,LATS2=lats1_2D, \
                    DATA_RANGE=data1[var]['drange'],                            \
                    NLEVELS=data1[var]['NLEVELS'], MPL_CBAR=data1[var]['cbar'], \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',               \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title+' mean', CBAR_LABEL=cbar_title,       \
                    iDISPLAY='N', FILE_PLOT=plotname,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='h',                       \
                    PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)
        #
        # Plot max map
        plotname=outdir+'EMEP4UK_'+str(var)+'_map_max.png'
        PT.plot_map(data2[var]['map_max'],lons2_2D,lats2_2D,                      \
                    DATA2=data1[var]['map_max'],LONS2=lons1_2D,LATS2=lats1_2D,    \
                    DATA_RANGE=data1[var]['drange'],                                        \
                    NLEVELS=data1[var]['NLEVELS'], MPL_CBAR=data1[var]['cbar'],   \
                    MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',                 \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                   \
                    PLOT_TITLE=plot_title+' max', CBAR_LABEL=cbar_title,          \
                    iDISPLAY='N', FILE_PLOT=plotname,                             \
                    LATDEL=2., LONDEL=2., RESOLUTION='h',                         \
                    PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange)
        #
        #
        # Plot time-series of all data
        # Put Data into list
        print time1.shape
        print data1[var]['TS_mean'].shape
        print data1[var]['TS_max'].shape
        print time2.shape
        print data2[var]['TS_mean'].shape
        print data2[var]['TS_max'].shape
        PLOT_DATA  = [ data1[var]['TS_mean'], data1[var]['TS_max'],  \
                       data2[var]['TS_mean'], data2[var]['TS_max']   ]
        TIME_DATA  = [ time1, time1, \
                       time2, time2  ]
        colours    = ['red','red','black','black']
        linestyles = ['-','--','-','--']
        plotname   = outdir+''+str(var)+'_timeseries.png'
        PT.plot_timeseries(PLOT_DATA, TIME_DATA, \
                           PLOT_TITLE=plot_title, \
                           FONTSIZES=[15,15,18,15], \
                           COLOURS=colours, LINESTYLES=linestyles, \
                           FILE_PLOT=plotname,HEIGHT=4)


        

