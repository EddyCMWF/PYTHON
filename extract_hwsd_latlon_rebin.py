#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import getpass
import datetime as dt
from scipy import stats

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
                
		# optional
		parser.add_argument('--inmissing',type=float,help='Missing value in input file',required=False,default=None)
		parser.add_argument('--outmissing',type=float,help='Missing value in output file',required=False,default=None)
                parser.add_argument('--lat_binsize',type=float,help='Binsize of latitudinal aggregation',required=False,default=None)
		parser.add_argument('--lon_binsize',type=float,help='Binsize of longitudinal aggregation',required=False,default=None)
                
		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('latmin',help='Minimum latitude',type=float)
		parser.add_argument('latmax',help='Maximum latitude',type=float)
		parser.add_argument('lonmin',help='Minimum longitude',type=float)
		parser.add_argument('lonmax',help='Maximum longitude',type=float)
                
		# Parse the arguments
		args=parser.parse_args()
                
		return args.infile, \
				args.outfile, \
				args.latmin, \
				args.latmax, \
				args.lonmin, \
				args.lonmax, \
				args.inmissing, \
				args.outmissing, \
                                args.lat_binsize, \
                                args.lon_binsize

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
        
	infile, outfile, latmin, latmax, lonmin, lonmax, in_mv, out_mv, lat_binsize, lon_binsize = e.parse_input()
        # '/users/eow/edwcom/data/HWSD/hwsd.bil', '/users/eow/edwcom/data/HWSD/test.nc', \
        # 20., 90., -180., 180., 0.0, -9999.0, 5, 5

	# Open input file
	print "Read input file: %s"%infile
        
	# Get header information
	hfile=infile[:-4]+'.hdr'
	print "Header in file: %s"%hfile
	hdr={}
	for lin in open(hfile,'r').readlines():
		name,value = lin.strip().split()
		if name=='BYTEORDER' or name=='LAYOUT':
			hdr[name]=value
		else:
			hdr[name]=int(value)
        
	# The grid definition
	gfile=infile[:-4]+'.blw'
	print "Grid definition in file: %s"%hfile
        
	dlon,rot1,rot2,dlat,lon0,lat0=[float(f.strip()) for f in open(gfile,'r').readlines() ]
        
	# Set up grid
	# These will be the CENTRES of the grid squares
	lon=np.arange(lon0,lon0+(hdr['NCOLS']*dlon),dlon)
	lonrange=[lon[0]-(dlon/2),lon[-1]+(dlon/2)]
        
	lat=np.arange(lat0,lat0+(hdr['NROWS']*dlat),dlat)
	latrange=[lat[-1]+(dlat/2),lat[0]-(dlat/2)]
        
	# finally read the input data
	indata=np.fromfile(infile,dtype='int16').reshape([hdr['NROWS'],hdr['NCOLS']]).astype(float)
	# mask where data=missing value
	if in_mv!=None:
		indata=np.ma.masked_equal(indata,in_mv)
        
        # Cut down to what we need
	xmin=np.argmin(np.abs(lon-lonmin))
	if lon[xmin]<lonmin:
		xmin+=1
        
	xmax=np.argmin(np.abs(lon-lonmax))
	if lon[xmax]>lonmax:
		xmax-=1
        	
	ymin=np.argmin(np.abs(lat-latmax))
	if lat[ymin]>latmax:
		ymin-=1
        
	ymax=np.argmin(np.abs(lat-latmin))
	if lat[ymax]<latmin:
		ymax-=1
        
        ###################################################################
        # ECP MOD
        # If rebinning ensure dimensions are multiples of bins, 
        # extend nort/east where neccessary
        if (lat_binsize!=None)&(lon_binsize!=None):       
                xmax=((np.ceil(float(xmax-xmin)/lon_binsize))*lon_binsize)+xmin
                ymax=((np.ceil(float(ymax-ymin)/lat_binsize))*lat_binsize)+ymin
        ###################################################################

        outdata=indata[ymin:ymax,xmin:xmax]
	outlat=lat[ymin:ymax]
	outlon=lon[xmin:xmax]
        del indata
        del lat
        del lon
	#outdataindx=np.ma.masked_equal(np.ones_like(outdata)*in_mv,in_mv)
	#j=1
	#for i in np.unique(outdata):
		#outdataindx[np.where(outdata==i)]=j
		#j+=1

	#plt.imshow(outdata.data)
	#plt.colorbar()
	#plt.grid(True)
	#plt.savefig('/users/eow/edwcom/test_plots/HWSD_fullres.png')

        ###################################################################
        # ECP Modification:
        # regrid to new latitude/longitude grid if required
        if (lat_binsize!=None)&(lon_binsize!=None):
                nlats = len(outlat)
                nlons = len(outlon)
                tempdata=outdata.reshape([nlats/lat_binsize,lat_binsize,\
                                          nlons/lon_binsize,lon_binsize])
                tempdata=tempdata.transpose(0,2,1,3)
                tempdata=tempdata.reshape( (nlats/lat_binsize), (nlons/lon_binsize), (lat_binsize*lon_binsize) )
                mode=stats.mode(tempdata.data,axis=2)[0]
                del outdata
                del tempdata
                outdata=np.ma.masked_equal(mode,0).reshape( (nlats/lat_binsize), (nlons/lon_binsize) )
                del mode
                outlat=np.array(np.median(outlat.reshape(nlats/lat_binsize,lat_binsize),axis=1))
                outlon=np.array(np.median(outlon.reshape(nlons/lon_binsize,lon_binsize),axis=1))
                
         #       plt.imshow(outdata.data)
         #       plt.colorbar()
         #       plt.grid(True)
         #       plt.savefig('/users/eow/edwcom/test_plots/HWSD_rebinned.png')
        
        ###################################################################
        
	# Output file
	print "Opening output file: %s"%outfile

	outf=nc.Dataset(outfile,'w')

	outf.createDimension('lon',len(outlon))
	outf.createVariable('lon','double',('lon',))
	outf.variables['lon'][:]=outlon[:]
	outf.variables['lon'].setncattr('long_name','longitude')
	outf.variables['lon'].setncattr('units','degrees east')

	outf.createDimension('lat',len(outlat))
	outf.createVariable('lat','double',('lat',))
	outf.variables['lat'][:]=outlat[:]
	outf.variables['lat'].setncattr('long_name','latitude')
	outf.variables['lat'].setncattr('units','degrees north')

	outf.createDimension('time',1)
	outf.createVariable('time','int32',('time',))
	outf.variables['time'][:]=0
	outf.variables['time'].setncattr('long_name','Seconds since 1900-01-01')
	outf.variables['time'].setncattr('units','seconds')

	outvar=outf.createVariable('mu','double',('time','lat','lon'),fill_value=out_mv)
	outvar[0,:]=outdata[:]
	outvar.long_name='mapping unit'
	outvar.units='none'

	user=getpass.getuser()
	now=dt.datetime.now().strftime('%d-%m-%Y')
	outf.title='Mapping unit extracted from %s on %s by %s'%(infile,now,user)

	outf.close()








	

