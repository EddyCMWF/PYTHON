#!/usr/bin/env python

import netCDF4 as nc
import sys
import numpy as np
import argparse

import EMEP4UK_tools as ET

import getpass
import datetime as dt

################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class extract:
        
	def parse_input(self):
                
		parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
                # positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('LSMfile',help='Land/Sea Mask  file')
                
		# Parse the arguments
		args=parser.parse_args()
                
		return args.infile, \
			args.outfile, \
                        arg.LSMfile

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

        #infile, outfile, latmin, latmax, lonmin, lonmax, in_mv, out_mv = e.parse_input()
        infile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Base_EUROPE_emep_4.3_2001_day.nc'
        outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_gridfile.nc'
        WRFfile='/users/eow/edwcom/EMEP/EMEP4UK/wrfout_d01_2001-12-30_00:00:00'
        
        inf=nc.Dataset(infile,'r')
        
        glon = inf.variables['lon'][:].flatten()
        glat = inf.variables['lat'][:].flatten()
        
        i = inf.variables['i_EMEP'][:]
        j = inf.variables['j_EMEP'][:]
        
        nx=len(i)
        ny=len(j)
        
        inf.close()
        
        gi,gj = np.meshgrid(i,j)
        gi    = gi.flatten()
        gj    = gj.flatten()
        
        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
        
        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
	
        WRFin = nc.Dataset(WRFfile,'r')
        gmask  = WRFin.variables['XLAND'][0,:,:].flatten()
        gmask[gmask==2] = 0
        
	# Open output file
        print "opening file: "+outfile
        outf=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')

	# create dimensions
        grid_size = len(glon)
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)

        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=[nx,ny]
        
        gctr_lat=outf.createVariable('i','float32',('grid_size',))
        gctr_lat.units='unitless'
        gctr_lat[:]=gi-np.amin(gi)

        gctr_lat=outf.createVariable('j','float32',('grid_size',))
        gctr_lat.units='unitless'
        gctr_lat[:]=gj-np.amin(gj)
        
        gctr_lat=outf.createVariable('grid_EMEP_i_coord','float32',('grid_size',))
        gctr_lat.units='unitless'
        gctr_lat[:]=gi

        gctr_lat=outf.createVariable('grid_EMEP_j_coord','float32',('grid_size',))
        gctr_lat.units='unitless'
        gctr_lat[:]=gj
        
        gctr_lat=outf.createVariable('grid_center_lat','float32',('grid_size',))
        gctr_lat.units='degrees'
        gctr_lat[:]=glat
        
        gctr_lon=outf.createVariable('grid_center_lon','float32',('grid_size',))
        gctr_lon.units='degrees'
        gctr_lon[:]=glon
        
        gim=outf.createVariable('grid_imask','i',('grid_size',))
        gim.units='unitless'
        gim[:]=gmask.astype('int')
        
        gcrn_lat=outf.createVariable('grid_corner_lat','float32',('grid_size','grid_corners'))
        gcrn_lat.units='degrees'
        gcrn_lat[:]=glatcorn
        
        gcrn_lon=outf.createVariable('grid_corner_lon','float32',('grid_size','grid_corners'))
        gcrn_lon.units='degrees'
        gcrn_lon[:]=gloncorn

        outf.title='EMEP4UK - Europe 50km grid'
        outf.conventions='SCRIP'
        outf.close()

        
        # Output to .csv file
        outf=open(outfile[:-3]+'.csv','w')
        outf.write('#%7s,%8s,%5s,%15s,%15s,%15s,%15s,%15s,%15s,%15s,%15s,%15s,%15s,\n' \
                   % ('i_EMEP','j_EMEP','Land','lon_cen','lat_cen',\
                      'lon_cor_a','lat_cor_a','lon_cor_b','lat_cor_b',\
                      'lon_cor_c','lat_cor_c','lon_cor_d','lat_cor_d') )
        for point in np.arange(grid_size):
                outf.write( '%8.2f,%8.2f,%5i,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,%15.8f,\n' \
                            % ( gi[point], gj[point], gmask[point], \
                                glon[point], glat[point], \
                                gloncorn[point,0], glatcorn[point,0], \
                                gloncorn[point,1], glatcorn[point,1], \
                                gloncorn[point,2], glatcorn[point,2], \
                                gloncorn[point,3], glatcorn[point,3]  )   )
        
        outf.close()
                           
