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
        outfile='/users/eow/edwcom/EMEP/EMEP4UK_EUROPE_gridfile'
        WRFfile='/users/eow/edwcom/EMEP/EMEP4UK/wrfout_d01_2001-12-30_00:00:00'
        
        inf=nc.Dataset(infile,'r')
        
        glon_full = inf.variables['lon'][:]
        glat_full = inf.variables['lat'][:]
        
        i = inf.variables['i_EMEP'][:]
        j = inf.variables['j_EMEP'][:]
        
        nx=len(i)
        ny=len(j)
        
        xmidpoint=80
        ymidpoint=60
        
        inf.close()
        
        gi_full,gj_full = np.meshgrid(i,j)
        
        WRFin = nc.Dataset(WRFfile,'r')
        gmask_full = WRFin.variables['XLAND'][0,:,:]
        gmask_full[gmask_full==2] = 0
        
        
        ###############################################################################################
        # Box 1
        glon      = glon_full[:xmidpoint,:ymidpoint]
        dims      = glon.shape
        glon      = glon.flatten()
        glat      = glat_full[:xmidpoint,:ymidpoint].flatten()
        gi        = gi_full[:xmidpoint,:ymidpoint].flatten()
        gj        = gj_full[:xmidpoint,:ymidpoint].flatten()
        gmask     = gmask_full[:xmidpoint,:ymidpoint].flatten()
        grid_size = len(glon)
        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
	
	# Open output file BOX 1
        print "opening file: "+outfile+'_box1.nc'
        outf=nc.Dataset(outfile+'_box1.nc','w',format='NETCDF3_CLASSIC')
	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)
        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=dims
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
        outf.title='EMEP4UK-EUROPE 50km grid'
        outf.conventions='SCRIP'
        outf.close()
        ###############################################################################################


        ###############################################################################################
        # Box 2
        glon      = glon_full[xmidpoint:,:ymidpoint]
        dims      = glon.shape
        glon      = glon.flatten()
        glat      = glat_full[xmidpoint:,:ymidpoint].flatten()
        gi        = gi_full[xmidpoint:,:ymidpoint].flatten()
        gj        = gj_full[xmidpoint:,:ymidpoint].flatten()
        gmask     = gmask_full[xmidpoint:,:ymidpoint].flatten()
        grid_size = len(glon)
        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
	
	# Open output file BOX 2
        print "opening file: "+outfile+'_box2.nc'
        outf=nc.Dataset(outfile+'_box2.nc','w',format='NETCDF3_CLASSIC')
	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)
        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=dims
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
        outf.title='EMEP4UK-EUROPE 50km grid'
        outf.conventions='SCRIP'
        outf.close()
        ###############################################################################################


        
        ###############################################################################################
        # Box 3
        glon      = glon_full[:xmidpoint,ymidpoint:]
        dims      = glon.shape
        glon      = glon.flatten()
        glat      = glat_full[:xmidpoint,ymidpoint:].flatten()
        gi        = gi_full[:xmidpoint,ymidpoint:].flatten()
        gj        = gj_full[:xmidpoint,ymidpoint:].flatten()
        gmask     = gmask_full[:xmidpoint,ymidpoint:].flatten()
        grid_size = len(glon)
        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
	
	# Open output file BOX 3
        print "opening file: "+outfile+'_box3.nc'
        outf=nc.Dataset(outfile+'_box3.nc','w',format='NETCDF3_CLASSIC')
	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)
        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=dims
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
        outf.title='EMEP4UK-EUROPE 50km grid'
        outf.conventions='SCRIP'
        outf.close()
        ###############################################################################################


        
        ###############################################################################################
        # Box 4
        glon      = glon_full[xmidpoint:,ymidpoint:]
        dims      = glon.shape
        glon      = glon.flatten()
        glat      = glat_full[xmidpoint:,ymidpoint:].flatten()
        gi        = gi_full[xmidpoint:,ymidpoint:].flatten()
        gj        = gj_full[xmidpoint:,ymidpoint:].flatten()
        gmask     = gmask_full[xmidpoint:,ymidpoint:].flatten()
        grid_size = len(glon)
        cor_i = np.array([ gi-0.5 , gi+0.5 , gi+0.5 , gi-0.5 ]).transpose(1,0)
        cor_j = np.array([ gj-0.5 , gj-0.5 , gj+0.5 , gj+0.5 ]).transpose(1,0)
        gloncorn,glatcorn = ET.xy_to_lonlat(cor_i,cor_j)
	
	# Open output file BOX 4
        print "opening file: "+outfile+'_box4.nc'
        outf=nc.Dataset(outfile+'_box4.nc','w',format='NETCDF3_CLASSIC')
	# create dimensions
        outf.createDimension('grid_size',grid_size)
        outf.createDimension('grid_corners',4)
        outf.createDimension('grid_rank',2)
        # create variables
        gdims=outf.createVariable('grid_dims','i',('grid_rank',))
        gdims[:]=dims
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
        outf.title='EMEP4UK-EUROPE 50km grid'
        outf.conventions='SCRIP'
        outf.close()
        ###############################################################################################


