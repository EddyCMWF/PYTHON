#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import getpass
import datetime as dt
from scipy import stats
from osgeo import gdal

################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class extract:
        #
	def parse_input(self):
                #
		parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
                #                
                #
		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('emepfile',help='EMEP file')
                #
		# Parse the arguments
		args=parser.parse_args()
                #
		return args.infile, \
                        args.outfile, \
                        args.emepfile
        #
        def minmax_latlon(self,emepfile):
                #
                inf=nc.Dataset(emepfile,'r')
                lats=inf.variables['grid_corner_lat'][:]
                lons=inf.variables['grid_corner_lon'][:]
                latmin= np.floor(np.amin(lats))-1.
                latmax= np.ceil(np.amax(lats))+1.
                lonmin= np.floor(np.amin(lons))-1.
                lonmax= np.ceil(np.amax(lons))+1.
                inf.close()
                #
                return latmin, latmax, lonmin, lonmax
        
                
################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################


if __name__=='__main__':
        
	# Call the class
	e=extract()
        #e.parse_input()
	infile, outfile, emepfile = e.parse_input()  #\
        #'/prj/wetlands_africa/topo/ti_marthews/topidx.tif', \
        #'/users/eow/edwcom/EMEP/topo2emep/x7y0/input_data/topidx.nc', \
        #'/users/eow/edwcom/EMEP/topo2emep/x7y0/input_data/grid_files/EMEP4UK_EUROPE_gridfile.nc', \
        #emepfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_gridfile_x1y1.nc'
        
        in_mv,out_mv=-9999.0, -9999.0
        latmin,latmax,lonmin,lonmax=e.minmax_latlon(emepfile)
        print 'lats: ',latmin,latmax
        print 'lons: ', lonmin,lonmax
        #
        orig_res = 15./3600.  # 15 arc-seconds
        #     
        lat=np.arange(latmax-(orig_res/2.),latmin,0-orig_res)
        lon=np.arange(lonmin+(orig_res/2.),lonmax,orig_res)
        #
        x_start=int(((lonmin+180.)/orig_res))
        x_end  =int(((lonmax+180.)/orig_res))
        
        y_start=int(((176.-(latmax+90.))/orig_res))
        y_end  =int(((176.-(latmin+90.))/orig_res))
        print 'x: ', x_start, x_end
        print 'y: ', y_start, y_end
        
	# Open input file
        print "Read input file: %s"%infile
        #
        inf=gdal.Open(infile, gdal.GA_ReadOnly)
        #
        # image size = 86400x34189.
        indata=inf.ReadAsArray( xoff=x_start,yoff=y_start, xsize=(x_end-x_start), ysize=(y_end-y_start))
        del inf
        indata[indata<-9999.]=-9999.0
        indata=np.ma.masked_equal(indata,in_mv)
        
        # Output file
    	print "Opening output file: %s"%outfile
        
	    outf=nc.Dataset(outfile,'w')
        
        outf.createDimension('lon',len(lon))
        outf.createVariable('lon','double',('lon',))
        outf.variables['lon'][:]=lon[:]
        outf.variables['lon'].setncattr('long_name','longitude')
        outf.variables['lon'].setncattr('units','degrees east')
        
        outf.createDimension('lat',len(lat))
        outf.createVariable('lat','double',('lat',))
        outf.variables['lat'][:]=lat[:]
        outf.variables['lat'].setncattr('long_name','latitude')
        outf.variables['lat'].setncattr('units','degrees north')
        
        outf.createDimension('time',1)
        outf.createVariable('time','int32',('time',))
        outf.variables['time'][:]=0
        outf.variables['time'].setncattr('long_name','Seconds since 1900-01-01')
        outf.variables['time'].setncattr('units','seconds')
        
        outvar=outf.createVariable('mu','double',('time','lat','lon'),fill_value=out_mv)
        outvar[0,:]=indata[:]
        outvar.long_name='mapping unit'
        outvar.units='none'
        
        user=getpass.getuser()
        now=dt.datetime.now().strftime('%d-%m-%Y')
        outf.title='Mapping unit extracted from %s on %s by %s'%(infile,now,user)
        
        outf.close()








	

