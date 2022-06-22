#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To remap data based on SCRIP output
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
import argparse
import time
#
import remap_tools as RT
#
##################################################################################
#  DEFINE CLASS
##################################################################################
class remap:
    ##############################################################################
    #  parse input
    ##############################################################################
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Remap data according to SCRIP output')
        #
        # Optional
        #
        parser.add_argument('--method',type=str,help='Method of Remapping', \
                                       required=False,default='Conserv')
        parser.add_argument('--order',type=int,help='order of remapping (1,2 or 3)', \
                                       required=False,default=1)
        parser.add_argument('--STDDEV',type=bool,help='Flag to output weighted standard deviation', \
                                       required=False,default=False)
        parser.add_argument('--src_data_name',type=str,help='Name of data parameter in the source data file', \
                                              required=False,default='mu')
        parser.add_argument('--dst_data_name',type=str,help='Name of data parameter for the destination data file', \
                                              required=False,default='mu')
        parser.add_argument('--dst_data_longname',type=str,help='Long Name of data parameter for the destination data file', \
                                                   required=False,default='Mapping Unit')
        parser.add_argument('--out_mv',type=float,help='Missing data value for output', \
                                        required=False,default=-9999.0)
        parser.add_argument('--out_title',type=str,help='Title for output', \
                                          required=False,default='Remapped data using SCRIP output.')
        parser.add_argument('--setUW',type=float,help='Value to set weights less than zero to', \
                                          required=False,default=None)
        # positional
	parser.add_argument('remapfile',help='Remap File')
        parser.add_argument('src_datafile',help='Source Data File')
	parser.add_argument('outfile',help='Output File')
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.remapfile,        \
               args.src_datafile,     \
               args.outfile,          \
               args.method,           \
               args.order,            \
               args.STDDEV,           \
               args.setUW,            \
               args.src_data_name,    \
               args.dst_data_name,    \
               args.dst_data_longname,\
               args.out_mv,           \
               args.out_title
        #
       


################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################

if __name__=='__main__':
    
    # Call the class
    remap=remap()
    
    remapfile, src_datafile, outfile, \
        method, order, STDDEV, setUW, \
        src_data_name, dst_data_name, dst_data_longname, \
        out_mv, out_title, \
        = remap.parse_input()
    
    # Read remap file parameters
    inf_remap     = nc.Dataset(remapfile,'r')
    src_grid_dims = inf_remap.variables['src_grid_dims'][:]
    dst_grid_dims = inf_remap.variables['dst_grid_dims'][:]
    src_mask      = inf_remap.variables['src_grid_imask'][:]
    dst_mask      = inf_remap.variables['dst_grid_imask'][:]
    src_add       = inf_remap.variables['src_address'][:] - 1   # minus one because python index from zero
    dst_add       = inf_remap.variables['dst_address'][:] - 1   # minus one because python index from zero
    map_wts       = inf_remap.variables['remap_matrix'][:]
    dst_cen_lats  = inf_remap.variables['dst_grid_center_lat'][:]
    dst_cen_lons  = inf_remap.variables['dst_grid_center_lon'][:]
    dst_cen_lons[dst_cen_lons>180.]=dst_cen_lons[dst_cen_lons>180.]-360.  # set to -180 to 180 range
    #= inf_remap.variables[''][:]
    
    src_grid_size = len(inf_remap.dimensions['src_grid_size'])
    dst_grid_size = len(inf_remap.dimensions['dst_grid_size'])
    num_links     = len(inf_remap.dimensions['num_links'])
    num_wgts      = len(inf_remap.dimensions['num_wgts'])
    
    inf_remap.close()
    
    #Read source data
    inf_srcdata = nc.Dataset(src_datafile,'r')
    src_array   = inf_srcdata.variables[src_data_name][:].flatten()
    src_units   = inf_srcdata.variables[src_data_name].units
    inf_srcdata.close()
    
    if ( len(src_array) != src_grid_size ):
        print 'source arrays do not match!'
    
    # if setUW keyword set set tl 0 wgts to zero
    if (setUW):
        map_wts[map_wts<0.]=setUW
    
    if (method=='Conserv'):
        print time.asctime()
        remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add, \
                                     dst_grid_size=dst_grid_size, \
                                     order=order, STDDEV=STDDEV, out_mv=out_mv )
    if (STDDEV):
        dst_array=remapped_data[0]
        dst_std_array=remapped_data[1]
    else:
        dst_array=remapped_data
    # insert alternative remapping methods here as elif statements.
    # The methods could/should be modulised
    #
    #
    #############################################
    #  REMOVED and done in subsequent program.
    #   Apply destination mask
    #  dst_array[dst_mask.reshape(dst_grid_dims[::-1])==0.]=-9999.
    #############################################
    
    dst_cen_lats.shape=dst_grid_dims[::-1]   # Reversed to account for Fortran silly
    dst_cen_lons.shape=dst_grid_dims[::-1]   # Reversed to account for Fortran silly
    dst_array.shape=dst_grid_dims[::-1]      # Reversed to account for Fortran silly
    if (STDDEV):
        dst_std_array.shape=dst_grid_dims[::-1]   # Reversed to account for Fortran silly
    
    # Write out data to netcdf file
    print "opening file: "+outfile
    outf=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
    
    # Create Dimensions
    outf.createDimension('x',dst_grid_dims[0])
    outf.createDimension('y',dst_grid_dims[1])
    # outf.createDimension('t', unknown)     #Include a time dimension at some point?
    
    # create variables
    # Lon
    lon=outf.createVariable('lon','float32',('y','x'))
    lon.units='degrees'
    lon.long_name='Longitude'
    lon.valid_min=np.amin(dst_cen_lons)
    lon.valid_max=np.amax(dst_cen_lons)
    lon[:]=dst_cen_lons
    # Lat
    lat=outf.createVariable('lat','float32',('y','x'))
    lat.units='degrees'
    lat.long_name='Latitude'
    lat.valid_min=np.amin(dst_cen_lats)
    lat.valid_max=np.amax(dst_cen_lats)
    lat[:]=dst_cen_lats
    
    outdata=outf.createVariable(dst_data_name,'float32',('y','x'))
    outdata.units=src_units
    outdata.long_name=dst_data_longname
    outdata.missing_value=out_mv
    outdata[:]=dst_array
    
    if (STDDEV):
        outdata_std=outf.createVariable(dst_data_name+'_std','float32',('y','x'))
        outdata_std.units=src_units
        outdata_std.long_name=dst_data_longname+' Standard Deviation'
        outdata_std.missing_value=out_mv
        outdata_std[:]=dst_std_array
    
    outf.title=out_title
    outf.method=method
    outf.order=order
    outf.close()

    








