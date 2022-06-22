#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: extract required WRF data for JULES drive data and output as monthly files
#            convert the snow and rain fall to correct units i.e. kg m-2 s-1, not mm
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
import argparse
import time, datetime
import netcdftime as nctime
import glob
#
##################################################################################
#  DEFINE CLASS
##################################################################################
class input:
    ##############################################################################
    #  parse input
    ##############################################################################
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Remap data according to SCRIP output')
        #
        # Optional
        parser.add_argument('--grid',type=str,help='EMEP grid, "d01" (Europe) or "d03" (UK)',required=False,default="d01")
        parser.add_argument('--outfile_tag',type=str,help='tag for outpuf filenames',required=False,default="wrfout")
        #
        # positional
        parser.add_argument('wrf_dir',type=str,help='directory containing the raw wrf data')
        parser.add_argument('start',type=str,help='start month for conversion, "YYYYMM"')
        parser.add_argument('end',type=str,help='end month for conversion, "YYYYMM"')
        parser.add_argument('out_dir',type=str,help='Directory for monthly output files')
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.wrf_dir,  \
               args.start,    \
               args.end,      \
               args.out_dir,  \
               args.grid,     \
               args.outfile_tag
        #

if __name__=='__main__':
    
    input=input()
    
    wrf_dir,start,end,out_dir,grid_name,outfile_tag=input.parse_input()
    
    #wrf_dir='/users/eow/edwcom/EMEP/EMEP4UK/WRF_DATA'
    #start='200107'
    #end='200107'
    #out_dir='/users/eow/edwcom/EMEP/EMEP4UK/JULES_drive/'
    #grid_name="d03"
    #outfile_tag="wrfout"
    
    # add / to directories if necessary
    if (wrf_dir[-1]!='/'):
        wrf_dir=wrf_dir+'/'
    
    if (out_dir[-1]!='/'):
        out_dir=out_dir+'/'
    
    # Seperate start date to month and year
    start_month= start[4:]
    start_year = start[:4]
    
    # Seperate end date to month and year
    end_month= end[4:]
    end_year = end[:4]
    
    # Extra Params, None at present
    XY_params= ['XLAT','XLONG','XLAND']
    T_param  = 'Times'
    output_time_units='days since 2000-01-01 00:00:00'
    
    params   = ['SWDOWN','GLW','RAINNC',  'SNOWNC', 'U10','V10','PSFC','Q2', 'T2' ]
    conv_fact= [ 1.0,     1.0, 1./10800., 1./10800., 1.0,  1.0,  1.0,   1.0,  1.0 ]
    conv_offs= [ 0.0,     0.0,  0.0,       0.0,      0.0,  0.0,  0.0,   0.0,  0.0 ]
    new_units= ['W m-2', 'W m-2','kg m-2 s-1','kg m-2 s-1','m s-1','m s-1','Pa','kg kg-1', 'K' ]
    
    ####  ASSUMING THAT TIME IS FIRST DIMENSION IN WRF FILE
    #### name of Time dimension
    ####T_dimname = 'Time'
    
    #loop round years
    for year in range(int(start_year),int(end_year)+1):
        # create list of months to loop round for eah year
        # if only one year start and end months can vary
        if (start_year==end_year):
            month_range=range(int(start_month),int(end_month)+1)
        # Start month for first year can vary
        elif year==int(start_year):
            month_range=range(int(start_month),12)
        # End month for first year can vary
        elif year==int(end_year):
            month_range=range(1,end_month+1)
        # all other years start at 1 and end at 12
        else:
            month_range=range(1,12)
                        
        for month in month_range:
            # add zero to month string if necessary
            if (month<=9):
                month_string='0'+str(month)
            else:
                month_string=str(month)
            
            # create file tag to search for files
            input_file_tag='wrfout_'+grid_name+'_'+str(year)+'-'+month_string+'*.nc'
            # list files to concacanate
            input_files=glob.glob(wrf_dir+input_file_tag)
            
            # flag to check if we are looking at first file in order to extract metadata
            first_file_check=True
            
            # List to store data for all relevant params from all files for month
            DATA = [ [] for param in params]
            DIMS = []
            TIMES= []
            
            # Extract Data from file and store in list
            for infile in input_files:
                #print infile
                inf=nc.Dataset(infile,'r')
                for param_num in range(len(params)):
                    DATA[param_num].append(inf.variables[params[param_num]][:])
                    if (first_file_check):
                        DIMS.append(inf.variables[params[param_num]].dimensions)
                
                # Read in Times
                for Tstr in inf.variables[T_param][:]:
                    TIMES.append( datetime.datetime.strptime(''.join(Tstr),'%Y-%m-%d_%H:%M:%S') )
                                
                inf.close()
                first_file_check=False
            
            # concatanate along Time dimension, i.e. first dimension
            # intitiate list
            CONC_DATA=[ [] for param in params]
            for param_num in range(len(params)):
                old_dims = DATA[param_num][0].shape
                # calculate new dimensions, Time dimension should be number of file x tsteps per file
                new_dims = [ len(input_files)*old_dims[0],old_dims[1],old_dims[2] ]
                # Reshape to new dims
                CONC_DATA[param_num]=np.array(DATA[param_num]).reshape(new_dims)
            
            del DATA
            
            # Apply conversaion factors and offsets
            CONVERTED_DATA = [ (np.array(DAT)*factor)+offset for DAT,factor,offset in zip(CONC_DATA,conv_fact,conv_offs) ]
            del CONC_DATA
            
            # output to netCDF
            OUTFILE=out_dir+outfile_tag+'_'+grid_name+'_'+str(year)+'-'+month_string+'.nc'
            
            # open output file to write
            outf=nc.Dataset(OUTFILE,'w')
            
            # open first input file for metadata
            inf=nc.Dataset(input_files[0],'r')
            
            for DIM,new_dim in zip(DIMS[0],new_dims):
                outf.createDimension(DIM,new_dim)
            
            # Write Time data
            outvar=outf.createVariable(T_param,'float32',(DIMS[0][0]))
            outvar.units=output_time_units
            outvar[:]=nctime.date2num(TIMES,units=output_time_units)
            
            # loop through XY params from first input file and store in outfile
            # Data is replicated for each time step so we just take the first time step
            for XYp in XY_params:
                # Create variable with dimensions 1 and 2
                outvar=outf.createVariable(XYp, 'float64', (DIMS[0][1:]) )
                #copy attributes directly from in file
                for attr in inf.variables[XYp].ncattrs():
                   outvar.setncattr(attr,inf.variables[XYp].getncattr(attr)) 
                # copy first time step of data directly from file
                outvar[:]=inf.variables[XYp][0,:,:]
                
            # loop through each parameter and store data in outfile
            for param_num in range(len(params)):
                param=params[param_num]
                outvar=outf.createVariable(param,'float64',(DIMS[param_num]))
                # copy attributes directly from in file apart from units which we set to our new units
                for attr in inf.variables[param].ncattrs():
                    if (attr=='units'):
                        outvar.setncattr(attr,new_units[param_num])
                    else:
                        outvar.setncattr(attr,inf.variables[param].getncattr(attr))
                outvar[:]=CONVERTED_DATA[param_num]
                        
            outf.title='WRF parameters required for JULES driving data'
            outf.author='E. Comyn-Platt'
            outf.history='Created: '+time.strftime('%x')
            outf.close()
            
