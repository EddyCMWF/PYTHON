#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To remap start-dump data based on SCRIP output
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
               args.outfile
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
        = remap.parse_input()
    
    #remapfile='/users/eow/edwcom/EMEP/wfd2emep/input_data/remap_files/rmp_wfdei_TO_emep4uk_EUROPE.nc'
    #src_datafile='/users/eow/edwcom/EMEP/wfd2emep/input_data/WFD_EI_global_2D.dump.20000101.0.nc'
    #outfile='/users/eow/edwcom/EMEP/wfd2emep/wfdei_emep4ukEUROPE_startdump.nc'
    
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
    
    # Index of surface type conversions, i.e. closest JULES surf type for each EMEP surf type:
    SurfTypeConvIndex = [2,1,2,1,3,3,5,3,3,5,8,8,8,7,9,6]
    
    
    # read in each start dump variable and regrid according to number of dims for variable
    inf_srcdata = nc.Dataset(src_datafile,'r')
    
    # Read in dimensions
    lon_dim    = len(inf_srcdata.dimensions['lon'])
    lat_dim    = len(inf_srcdata.dimensions['lat'])
    tile_dim   = len(inf_srcdata.dimensions['tile'])
    scpool_dim = len(inf_srcdata.dimensions['scpool'])
    soil_dim   = len(inf_srcdata.dimensions['soil'])
    snow_dim   = len(inf_srcdata.dimensions['snow'])
    
    
    # Open outfile to write to
    outf=nc.Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    
    # Create Dimensions
    outf.createDimension('west_east',dst_grid_dims[0])
    outf.createDimension('south_north',dst_grid_dims[1])
    outf.createDimension('tile',len(SurfTypeConvIndex))   #number of EMEP surf types
    outf.createDimension('scpool',scpool_dim)
    outf.createDimension('z',soil_dim)
    outf.createDimension('snow',snow_dim)
    
    # loop through params with dimensions: [tile,lat,lon]
    tile_params=['canopy','snow_tile','tstar_tile','rho_snow','snow_depth','snow_grnd','nsnow']
    for src_data_name in tile_params:
        print src_data_name
        src_array_in= inf_srcdata.variables[src_data_name][:]
        # put each tile array into seperatre list element
        src_list = [src_array_in[i,:,:].flatten() for i in range(tile_dim)]
        del src_array_in
        # open list to store destination data
        dst_data_list=[]
        # loop through each input tile
        print 'looping through each JULES tile'
        for src_array in src_list:
            if ( len(src_array) != src_grid_size ):
                print 'source arrays do not match!'
            remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add,\
                                            dst_grid_size=dst_grid_size  )
            # reshape data to original 2D dims and append to dst list
            dst_data_list.append(remapped_data.reshape(dst_grid_dims[::-1]))
            del remapped_data
        del src_array
        del src_list
        # Re-index tiles onto EMEP tiles (-1 to account for index from zero)
        dst_out_list=[dst_data_list[index-1] for index in SurfTypeConvIndex]
        del dst_data_list
        outdata=outf.createVariable(src_data_name,'float32',('tile','south_north','west_east'))
        outdata[:]=np.array(dst_out_list,dtype='float32')
        del dst_out_list    
    
    # Soil Carbon data has dimensions: [scpool,lat,lon]
    # scpool_dim=1, so just geo remapping is required
    src_data_name='cs'
    print src_data_name
    src_array= inf_srcdata.variables[src_data_name][:].flatten()
    if ( len(src_array) != src_grid_size ):
        print 'source arrays do not match!'
    remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add,\
                                    dst_grid_size=dst_grid_size  )
    # reshape data to original 2D dims 
    dst_data=remapped_data.reshape(dst_grid_dims[::-1])
    del remapped_data
    del src_array
    outdata=outf.createVariable(src_data_name,'float32',('scpool','south_north','west_east'))
    outdata[:]=dst_data
    del dst_data
    
    # Loop through params with dimensions: [lat,lon]
    geo_params=['gs','sthzw','zw']
    # just geo remapping is required
    for src_data_name in geo_params:
        print src_data_name
        src_array= inf_srcdata.variables[src_data_name][:].flatten()
        if ( len(src_array) != src_grid_size ):
            print 'source arrays do not match!'
        remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add,\
                                        dst_grid_size=dst_grid_size  )
        # reshape data to original 2D dims 
        dst_data=remapped_data.reshape(dst_grid_dims[::-1])
        del remapped_data
        del src_array
        outdata=outf.createVariable(src_data_name,'float32',('south_north','west_east'))
        outdata[:]=dst_data
        del dst_data
 
    # loop through params with dimensions: [soil,lat,lon]
    soil_params=['sthuf','t_soil']
    for src_data_name in soil_params:
        print src_data_name
        src_array_in= inf_srcdata.variables[src_data_name][:]
        # put each tile array into seperatre list element
        src_list = [src_array_in[i,:,:].flatten() for i in range(soil_dim)]
        del src_array_in
        # open list to store destination data
        dst_data_list=[]
        # loop through each input tile
        print 'looping through each Soil layer'
        for src_array in src_list:
            if ( len(src_array) != src_grid_size ):
                print 'source arrays do not match!'
            remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add,\
                                            dst_grid_size=dst_grid_size  )
            # reshape data to original 2D dims and append to dst list
            dst_data_list.append(remapped_data.reshape(dst_grid_dims[::-1]))
            del remapped_data
        del src_array
        del src_list
        # Re-index tiles onto EMEP tiles (-1 to account for index from zero)
        outdata=outf.createVariable(src_data_name,'float32',('z','south_north','west_east'))
        outdata[:]=np.array(dst_data_list,dtype='float32')
        del dst_data_list
    
    # loop through params with dimensions: [snow,tile,,lat,lon]
    snow_tile_params=['snow_ds','snow_ice','snow_liq','tsnow']
    for src_data_name in snow_tile_params:
        print src_data_name
        dst_data_list=[]
        for i in range(snow_dim):
            dst_data_list_temp=[]
            for j in range(tile_dim):
                src_array = inf_srcdata.variables[src_data_name][i,j,:,:].flatten()
                if ( len(src_array) != src_grid_size ):
                    print 'source arrays do not match!'
                remapped_data=RT.conserv_remap( map_wts,src_array,src_add,dst_add,\
                                                dst_grid_size=dst_grid_size  )
                # append for each JULES surf type
                dst_data_list_temp.append( remapped_data.reshape(dst_grid_dims[::-1]) )
                del remapped_data
                del src_array
            # We now reorder the data to match the EMEP surf types,
            # then append for each snow layer
            dst_data_list.append([dst_data_list_temp[index-1] for index in SurfTypeConvIndex])
            del dst_data_list_temp
        outdata=outf.createVariable(src_data_name,'float32',('snow','tile','south_north','west_east'))
        outdata[:]=np.array(dst_data_list,dtype='float32')
        del dst_data_list    

    inf_srcdata.close()
    
    outf.note='Surface type conversions from original JULES index = [2,1,2,1,3,3,5,3,3,5,8,8,8,7,9,6]'
    outf.close()

    








