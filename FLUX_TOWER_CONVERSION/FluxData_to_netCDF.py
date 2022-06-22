#!/bin/env python
######################################################################################
# Program: FluxData_to_netCDF.py
# Author: Edward Comyn-Platt (edwcom@ceh.ac.uk)
#
# Purpose: Convert .csv Flux Tower data into netCDF4 format.
#            csv file and must have an accomponying etadata file which maps the 
#            contents of the csv to 
# 
# Input files: csv data file -  Must conform to naming protocol '[dir]/[site].csv'
#                                or be provided explicitely with --datafile
#              csv metadata file - Must conform to naming protocol '[dir]/[site]_meta.csv'
#                                   or be provided explicitely with --metafile
#                                   file containing the metadata for output and the
#                                   mapping of datafile variable names to output 
#                                   variable names, and information on collating 
#                                   several variables into a single variable, e.g.
#                                   soil temperature profiles.
#                       Format of metafile:
#                           first line is header and ignored
#                           GLOBAL requirements - n_output_profiles
#                                               - time_name
#                                               - time_format
#                                               - out_time_units
#                                               - fill_value
#                           GLOBAL optional notes for netCDF output can be provided by preceding 
#                             note with "globalnote_"
#                           
#                          Each output profile must be enclose between to /
#                           PROFILE requirements - name
#                                                - nvars
#                           PROFILE optional notes for netCDF output can be provided by preceding 
#                             note with "profilenote_"
#                            DIMENSIONS - only spatial dimensions need to be defined (e.g. x and z)
#                                       - must contain atleast one spatial dimension
#                                       - identified by "# Dimension"
#                                       - each dimension needs name, units, length and values
#                            Variables - Variable section identified by "# output_name" 
#                                      - each line for a variable
#                                      - Must Contain:
#                                          - output_name
#                                          - local_name
#                                          - output_units
#                                          - spatial_dimensions
#                                          - [dim_name] (position for each dimension)
#                                      - Optionals:
#                                          - add_offset 
#                                          - scale_factor
#                                          - missing_value
#                           Variable optional notes for netCDF output can be provided by preceding 
#                             note with "varnote_"
#                                  
# Required Parsed Input: site - name of the site to be converted.
#
# Optional Parsed Input: --dir - Directory which contains the .csv file and 
#                                  accompanying metadata file. 
#                                 Default is current directory (./).
#                        --outdir - Directory to store the output netcdf file.
#                                 Default is current directory (./)
#                        --datafile - Full path and name of input data file if it doesn't
#                                      conform to standard naming protocol, i.e.:
#                                       [dir]/[site].csv
#                        --metafile - Full path and name of input metadata file if it doesn't
#                                      conform to naming protocol, i.e.:
#                                       [dir]/[site]_meta.csv
#
# Output: a netcdf file for each profile specified in the csv metadata file with the
#           naming protocol:
#              [outdir]/[site]_[profile_name].nc
# 
# Examples:
#   ./FluxData_to_netCDF.py [SITE]  (converts SITE data in current direcotry) 
#   ./FluxData_to_netCDF.py [SITE] --dir [DIR] (converts SITE data in DIR)
#   ./FluxData_to_netCDF.py [SITE] --dir [DIR] --outdir [OUTDIR] (converts SITE data in DIR and saves to OUTDIR)
#   ./FluxData_to_netCDF.py [SITE] --datafile [DATAFILE] --metafile [METAFILE] --outdir [OUTDIR] 
#                           (converts SITE data in DATAFILE using METAFILE and stores to OUTDIR)
#   
#####################################################################################

import csv
import netCDF4 as nc
import numpy as np
import argparse
import datetime as dt
import netcdftime as nctime

#################################################################
# Define Class
class FluxDataConvertor:
   
    ################################################################################
    # Parse input
    ################################################################################
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Convert Flux Tower csv data to netCDF')
        #
        # optional
        parser.add_argument('--dir',help='Directory containing the data and meta data file'+
                            'and where to write output file to.',default='./',required=False)
        parser.add_argument('--outdir',help='Alternative path for output netCDF file',\
                            default='None',required=False)
        parser.add_argument('--datafile',help='Full path and filename of data csv file',\
                            default='None',required=False)
        parser.add_argument('--metafile',help='Full path and filename of  metadata csv file',\
                            default='None',required=False)
       
        # positional
        parser.add_argument('site',help='Site name, must correspond to file naming convention'+\
                                         ' if filenames not set specifically')
               
        # Parse the arguments
        args=parser.parse_args()
        
        # add '/' to dir if neccessary
        if args.dir[-1:] != '/':
            args.dir+='/'
        
        # copy dir to outdir
        if args.outdir=='None':
            args.outdir=args.dir
        # add '/' to outdir if neccessary
        elif args.outdir[-1:] != '/':
            args.outdir+='/'
        
        # Create filename where neccessary
        if args.datafile=='None':
            args.datafile=args.dir+args.site+'.csv'
        if args.metafile=='None':
            args.metafile=args.dir+args.site+'_meta.csv'
        
        # return site name and filenames
        return args.site, \
            args.datafile, \
            args.metafile, \
            args.outdir

    
    ########################################################################################
    # Read metadata file
    ########################################################################################
    def find_line_start(self,string,lines):
        return np.where([string==line[:len(string)] for line in lines ])[0]
    
    def read_metadata(self,metafile):
        metafile_lines=open(metafile,'r').readlines()
        
        # Top line is header
        metafile_header=metafile_lines.pop(0)
        
        # get the number of output profiles
        temp_index    = self.find_line_start('n_output_profile',metafile_lines)
        n_out_profiles= int(metafile_lines.pop(temp_index).split(',')[1])
        ## get the latitude and longitude in the same way
        #temp_index = self.find_line_start('latitude',metafile_lines)
        #latitude   = float(metafile_lines.pop(temp_index).split(',')[1])
        #temp_index = self.find_line_start('longitude',metafile_lines)
        #longitude  = float(metafile_lines.pop(temp_index).split(',')[1])
                
        # time_name and format in the same way   
        temp_index = self.find_line_start('time_name',metafile_lines)
        time_name  = metafile_lines.pop(temp_index).split(',')[1]
        temp_index = self.find_line_start('time_format',metafile_lines)
        time_format= metafile_lines.pop(temp_index).split(',')[1]
        temp_index = self.find_line_start('out_time_units',metafile_lines)
        out_time_units= metafile_lines.pop(temp_index).split(',')[1]
        
        # do the same for the GLOBAL fill_value
        temp_index = self.find_line_start('fill_value',metafile_lines)
        fill_value  = metafile_lines.pop(temp_index).split(',')[1]

        # append this into a dictionary of global data
        GLOBAL_DATA = { 'n_out_profiles':n_out_profiles, \
                        #'latitude':latitude,       \
                        #'longitude':longitude,     \
                        'time_name':time_name,     \
                        'time_format':time_format, \
                        'out_time_units':out_time_units, \
                        'fill_value':fill_value, \
                        }
        
        # Deal with remaining header lines
        line=' '
        global_notes={}
        while line[:1]!='/':
            line=metafile_lines.pop(0)
            # if line begins with 'outputnote_' add to 
            if line[:11]=='globalnote_':
                split=line.split(',')
                global_notes[split[0][11:]]=split[1]
        
        GLOBAL_DATA['global_notes']=global_notes
        
        # Now read in profiles
        PROFILES=[]
        for iprof in range(n_out_profiles):
            profile={}
            line=metafile_lines.pop(0)
            while line[:1]!='#':
                split=line.split(',')
                profile[split[0]]=split[1]
                line=metafile_lines.pop(0)

            # read in spatial dimensions
            line=line.replace('#','').replace(' ','').replace('\r','').replace('\n','')
            dim_meta=line.split(',')
            dims=[]
            line=metafile_lines.pop(0).replace('\r','').replace('\n','')
            while line[:1]!='#':
                split=line.split(',')
                dim={}
                for imeta in range(len(dim_meta)):
                    if dim_meta[imeta]=='values':
                        dim[dim_meta[imeta]]=split[imeta:imeta+int(dim['length'])]
                    elif dim_meta[imeta]!='':
                        dim[dim_meta[imeta]]=split[imeta]
                dims.append(dim)
                line = metafile_lines.pop(0).replace('\r','').replace('\n','')

            profile['dims']=dims

            # line that begins with # is header line and marks the start of the variable metadata
            line=line.replace('#','').replace(' ','').replace('\r','').replace('\n','')
            var_meta=line.split(',')
            
            vars=[]
            line = metafile_lines.pop(0).replace('\r','').replace('\n','')
            while line[:1]!='/':
                split=line.split(',')
                var={}
                for imeta in range(len(var_meta)):
                    if var_meta[imeta]!='':
                        var[var_meta[imeta]]=split[imeta]
                vars.append(var)
                line = metafile_lines.pop(0).replace('\r','').replace('\n','')
            
            profile['vars']=vars
            
            PROFILES.append(profile)
    
        return GLOBAL_DATA, PROFILES
        
#################################################################
# Main Function
#
if __name__=='__main__':
    
    # call the class
    FC=FluxDataConvertor()
    
    site,datafile,metafile,outdir=FC.parse_input()
    #site,datafile,metafile,outdir='Berambadi_Example','./Berambadi_Example.csv','./Berambadi_Example_meta.csv','./'
    
    # Read meta data file
    GLOBAL_DATA,PROFILES=FC.read_metadata(metafile)
    
    fill_value=float(GLOBAL_DATA['fill_value'])
    
    # open and read csv file and data
    data_inf=open(datafile,'r')
    full_data=list(csv.reader(data_inf,delimiter=','))
    full_data_headers=full_data.pop(0)
    full_data_units=full_data.pop(0)
    data_inf.close()

    # read in the time data as a ddatetime object
    Tindex    = full_data_headers.index(GLOBAL_DATA['time_name'])
    time_data = np.array([dt.datetime.strptime(line[Tindex],GLOBAL_DATA['time_format']) \
                            for line in full_data] )
    time_length=len(time_data)
    # Create time array from out_base_time
    out_time_units=GLOBAL_DATA['out_time_units'].replace('START',str(time_data[0]))
    out_time_data = nctime.date2num(time_data,units=out_time_units)
    
    for profile in PROFILES:
        # Open netCDF file to write to
        outf=nc.Dataset(outdir+site+'_'+profile['name']+'.nc','w')
        
        # create time dimension and time variable
        outf.createDimension('time',time_length)
        outvar=outf.createVariable('time','float32',('time'))
        outvar.units=out_time_units
        outvar[:]=out_time_data
                
        # create spatial dimensions and spatial dimension variables
        for dim in profile['dims']:
            outf.createDimension(dim['Dimension'],int(dim['length']))
            outvar=outf.createVariable(dim['Dimension'],'float32',(dim['Dimension']))
            outvar.units=dim['units']
            outvar[:]=np.array([ float(dimval) for dimval in dim['values'] ])
                
        # Get data for output profile, collating where necessary
        # out_data dictionary for storing data in
        out_data_dict={}
        outvar_dict={}
        for var_meta in profile['vars']:
            # get/create offset
            if 'add_offset' in var_meta:
                if var_meta['add_offset']!='':
                    offset = float(var_meta['add_offset'])
                else:
                    offset=0.
            else:
                offset=0.

            # get/create scale factor
            if 'scale_factor' in var_meta:
                if var_meta['scale_factor']!='':
                    scale_factor = float(var_meta['scale_factor'])
                else:
                    scale_factor=1.
            else:
                scale_factor=1.
            print(var_meta['output_name'])
            out_name=var_meta['output_name']
            out_dim_names=['time']+var_meta['spatial_dimensions'].split('-')
            out_dims= [ len(outf.dimensions[out_dim_name]) for out_dim_name in out_dim_names ]
            
            # create data array and output variable in dictionaries if not already there
            if out_name not in out_data_dict:
                out_data_dict[out_name]=np.zeros(out_dims)
                outvar_dict[out_name]=outf.createVariable(out_name,'float32',\
                                                          tuple(out_dim_names),\
                                                          fill_value=fill_value)
                outvar_dict[out_name].units=var_meta['output_units']
                for meta in var_meta:
                    if (str(meta)[:8]=='varnote_')&(var_meta[meta]!='') :
                        outvar_dict[out_name].setncattr(meta[8:],var_meta[meta])
                
            # column number of data
            if var_meta['local_name']!='MISSING':
                Vindex= full_data_headers.index(var_meta['local_name'])
                # loop over each point and put into appropriate place in dictionary
                for itstep in range(time_length):
                    point = tuple( [itstep]+[ int(var_meta[dimname])-1 \
                                              for dimname in var_meta['spatial_dimensions'].split('-') ] )
                    #print(time_data[itstep],full_data[itstep][Vindex])
                    if full_data[itstep][Vindex]=='':
                       out_data_dict[out_name][point]=GLOBAL_DATA['fill_value']
                    elif float(full_data[itstep][Vindex])!=GLOBAL_DATA['fill_value']:
                        out_data_dict[out_name][point]=(float(full_data[itstep][Vindex])+offset)*scale_factor
                    else:
                        out_data_dict[out_name][point]=GLOBAL_DATA['fill_value']

            else: 
                # get/create fill_value
                if 'missing_value' in var_meta:
                    if var_meta['missing_value']!='':
                        missing_value = float(var_meta['missing_value'])
                    else:
                        missing_value=GLOBAL_DATA['fill_value']
                else:
                    missing_value=GLOBAL_DATA['fill_value']
                out_data_dict[out_name][:]=missing_value
                
        for out_name in out_data_dict:
            outvar_dict[out_name][:]=out_data_dict[out_name][:]
       
        for prof_meta in profile:
            if prof_meta[:12]=='profilenote_':
                outf.setncattr(prof_meta[12:],profile[prof_meta])

        outf.history=str(dt.date.today())+': netCDF file created using Edward Comyn-Platt (edwcom@ceh.ac.uk) proccessing script FluxData_to_netCDF.py'
        outf.close()
        
    


