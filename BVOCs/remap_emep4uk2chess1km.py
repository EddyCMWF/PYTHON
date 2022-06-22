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
import pickle
import os
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
        # Positional
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
       


################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################

if __name__=='__main__':
    
    # Call the class
    #remap=remap()
    # 
    #remapfile, src_datafile, outfile= remap.parse_input()
    
    remapfile='/users/eow/edwcom/EMEP/chess2emep/remap_files/rmp_emep4uk5km_TO_chess1km_incsea.nc'
    src_datafile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Base_emep_4.3_2001_day.nc'
    outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Base_emep_4.3_2001_day_ChessGrid2.nc'
    
    out_fv = -1e10
    
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
    src_grid_size = len(inf_remap.dimensions['src_grid_size'])
    dst_grid_size = len(inf_remap.dimensions['dst_grid_size'])
    num_links     = len(inf_remap.dimensions['num_links'])
    num_wgts      = len(inf_remap.dimensions['num_wgts'])
    inf_remap.close()
    
    # Array of unique destination points:
    uni_dst_add=set( dst_add )
    # Create list of indexes for each dst point
    picklefile=remapfile[:-2]+'pickleNP'
    if not (os.path.isfile(picklefile)):
        uber_index = [ np.where(dst_add==idx) for idx in uni_dst_add ]
        with open(picklefile, 'w') as savefile:
            pickle.dump(uber_index,savefile)
    else:
        with open(picklefile, 'r') as savefile:
            uber_index=pickle.load(savefile)
    
    # Open source file to read from
    inf_srcdata = nc.Dataset(src_datafile,'r')
    indim_i     = len(inf_srcdata.dimensions['i'])
    indim_j     = len(inf_srcdata.dimensions['j'])
    indim_t     = len(inf_srcdata.dimensions['time'])
    in_time_data= inf_srcdata.variables['time'][:]
    
    # Open destination file to write to
    outf = nc.Dataset(outfile,'w')
    # Create Dimensions
    outf.createDimension('i',dst_grid_dims[0])
    outf.createDimension('j',dst_grid_dims[1])
    outf.createDimension('time',indim_t)
    
    outvar=outf.createVariable('x','int64',('i',))
    outvar.longname="Eastings"
    outvar.units="metres"
    outvar[:]= (np.arange(dst_grid_dims[0],dtype='int64')*1000.)+500.
    
    outvar=outf.createVariable('y','int64',('j',))
    outvar.longname="Northings"
    outvar.units="metres"
    outvar[:]= (np.arange(dst_grid_dims[1],dtype='int64')*1000.)+500.
    
    outvar=outf.createVariable('time','int64',('time',))
    outvar.longname="Time at middle of period"
    outvar.units="days since 1900-1-1 0:0:0"
    outvar[:]= in_time_data
    
    outvar=outf.createVariable('cen_lat','float32',('i','j'))
    outvar.longname="Centre latitude of grid cell"
    outvar.units="Degrees North"
    outvar[:]=dst_cen_lats.reshape(dst_grid_dims[::-1])
    
    outvar=outf.createVariable('cen_lon','float32',('i','j'))
    outvar.longname="Centre longitude of grid cell"
    outvar.units="Degrees East"
    outvar[:]=dst_cen_lons.reshape(dst_grid_dims[::-1])
    
    k=0
    for var in inf_srcdata.variables:
        k+=1
        if (k<=10):
            continue
        print var
        invar   = inf_srcdata.variables[var]
        # read in data from src and flatten spatial dimensions
        indata  = invar[:].reshape(indim_t,indim_i*indim_j)
        
        # Calculate outdata by looping through uber_index within loop through days
        # This is remarkably quick given list looping
        outdata = np.array( [                                                      \
                     np.array([ np.sum( indata[daynum,src_add[idx]]*map_wts[idx,0])    \
                               / np.sum( map_wts[idx,0] ) for idx in uber_index ]) \
                                                   for daynum in range(indim_t) ]  )

        outdata_full=np.arange(indim_t*dst_grid_size).reshape(indim_t,dst_grid_size)+out_fv
        outdata_full[:,list(uni_dst_add)]=outdata
        del outdata
        del indata
        if np.any(outdata_full==out_fv):
            outdata_full=np.ma.masked_equal(outdata_full,out_fv)
        
        outvar  =  outf.createVariable(var,'float64',('time','i','j'),fill_value=out_fv)
        outvar.long_name = invar.long_name
        outvar.units     = invar.units
        outvar[:]        = outdata_full.reshape(indim_t,dst_grid_dims[1],dst_grid_dims[0])
     

    
    outf.startdate=invar.current_date_first
    outf.enddate=invar.current_date_last
    outf.period='Daily'
    outf.projection='Chess 1km Grid'
    outf.model='EMEP_MSC-W'
    outf.conventions='CF-1.0'
    outf.about='EMEP4UK 5km data reprojected on the Chess 1km grid using 1st order conervative remapping'
    outf.author_of_run='Unimod Group'
    outf.author_of_reprojection='Edward Comyn-Platt, edwcom@ceh.ac.uk'
    outf.close()

    
