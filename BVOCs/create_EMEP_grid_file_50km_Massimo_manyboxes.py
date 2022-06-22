#!/usr/bin/env python

import netCDF4 as nc
import sys
import numpy as np
import argparse

import EMEP4UK_tools as ET

import getpass
import datetime as dt

import matplotlib.pyplot as plt

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
#if (True):
	# Call the class
	#e=extract()
        #infile, outfile_base, latmin, latmax, lonmin, lonmax, in_mv, out_mv = e.parse_input()
        infile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Base_EUROPE_emep_4.3_2001_day.nc'
        outfile_base='/users/eow/edwcom/EMEP/EMEP4UK/sections/EMEP4UK_EUROPE_gridfile'
        WRFfile='/users/eow/edwcom/EMEP/EMEP4UK/wrfout_d01_2001-12-30_00:00:00'
        
        inf=nc.Dataset(infile,'r')
        
        glon_full = inf.variables['lon'][:]
        glat_full = inf.variables['lat'][:]
        
        i = inf.variables['i_EMEP'][:]
        j = inf.variables['j_EMEP'][:]
        
        inf.close()
        
        nx=len(i)
        ny=len(j)
        
        nxBoxes = 6.
        nyBoxes = 5.
        
        xinterval = np.ceil(nx/nxBoxes)
        yinterval = np.ceil(ny/nyBoxes)
        
        
        gi_full,gj_full = np.meshgrid(i,j)
        
        WRFin = nc.Dataset(WRFfile,'r')
        gmask_full = WRFin.variables['XLAND'][0,:,:]
        gmask_full[gmask_full==2] = 0
        plotmask  = gmask_full
        plotmask[glat_full>85]=-1
        WRFin.close()
        
        plt.figure(figsize=(15,15))
        listfile = open(outfile_base+'_list.txt','w')
        ###############################################################################################
        k=1
        for ybox in range(int(nyBoxes)):
           for xbox in range(int(nxBoxes)):
                # if statement to remove section not covered in topidx file, i.e. >85 lat
                #if (xbox==1) and (ybox==3):
                #        xstart   = np.floor((xinterval*xbox) + (xinterval/4.))
                #else:
                
                xstart = np.floor(xinterval * xbox)
                ystart = np.floor(yinterval * ybox)
                yend   = min(ystart+yinterval,ny)   #  for both x and y dimensions  
                
                # if statement to remove section not covered in topidx file, i.e. >85 lat
                #if (xbox==0) and (ybox==3):
                #        xend   = np.floor(min((xinterval * xbox)+(3.*xinterval/4.),nx))   # Cap end to length of dimension   
                #else:
                
                xend   = np.floor(min((xinterval * xbox)+xinterval,nx))   # Cap end to length of dimension
                
                glon      = glon_full[ystart:yend,xstart:xend]
                dims      = glon.shape
                glon      = glon.flatten()
                grid_size = len(glon)
                
                glat      = glat_full[ystart:yend,xstart:xend].flatten()
                gi        = gi_full[ystart:yend,xstart:xend].flatten()
                gj        = gj_full[ystart:yend,xstart:xend].flatten()
                gmask     = gmask_full[ystart:yend,xstart:xend].flatten()
                plt.subplot(nyBoxes,nxBoxes,k)
                plt.imshow(plotmask[ystart:yend,xstart:xend],vmin=-1,vmax=1)
                plt.title('x'+str(xbox)+' y'+str(ybox))
                k=k+1
                if (min(gmask)>=0) and (max(gmask)==1):
                        listfile.write('x'+str(xbox)+'y'+str(ybox)+'\n')
                        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
                        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
                        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
                        # Open output file
                        outfile = outfile_base+'_x'+str(int(xbox))+'y'+str(int(ybox))+'.nc'
                        print "opening file: "+outfile
                        outf=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
                        # create dimensions
                        outf.createDimension('grid_size',grid_size)
                        outf.createDimension('grid_corners',4)
                        outf.createDimension('grid_rank',2)
                        # create variables
                        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
                        gdims[:]=dims[::-1]
                        gctr_lat=outf.createVariable('i','float64',('grid_size',))
                        gctr_lat.units='unitless'
                        gctr_lat[:]=gi-np.amin(gi)
                        gctr_lat=outf.createVariable('j','float64',('grid_size',))
                        gctr_lat.units='unitless'
                        gctr_lat[:]=gj-np.amin(gj)
                        gctr_lat=outf.createVariable('grid_EMEP_i_coord','float64',('grid_size',))
                        gctr_lat.units='unitless'
                        gctr_lat[:]=gi
                        gctr_lat=outf.createVariable('grid_EMEP_j_coord','float64',('grid_size',))
                        gctr_lat.units='unitless'
                        gctr_lat[:]=gj
                        gctr_lat=outf.createVariable('grid_center_lat','float64',('grid_size',))
                        gctr_lat.units='degrees'
                        gctr_lat[:]=glat
                        gctr_lon=outf.createVariable('grid_center_lon','float64',('grid_size',))
                        gctr_lon.units='degrees'
                        gctr_lon[:]=glon
                        gim=outf.createVariable('grid_imask','i',('grid_size',))
                        gim.units='unitless'
                        gim[:]=gmask.astype('int')
                        gcrn_lat=outf.createVariable('grid_corner_lat','float64',('grid_size','grid_corners'))
                        gcrn_lat.units='degrees'
                        gcrn_lat[:]=glatcorn
                        gcrn_lon=outf.createVariable('grid_corner_lon','float64',('grid_size','grid_corners'))
                        gcrn_lon.units='degrees'
                        gcrn_lon[:]=gloncorn
                        outf.title='EMEP4UK-EUROPE 50km grid sub-section'
                        outf.conventions='SCRIP'
                        outf.close()
        ###############################################################################################
        listfile.close()
        plt.savefig('/users/eow/edwcom/test_plots/EMEP_50km_boxes.png')
