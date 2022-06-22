#!/usr/bin/env python

################################################################################
#
# Script to put the CHESS elev file onto netCDF
#
# ELR 16/05/2014
#
################################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse

################################################################################
################################################################################
#
# Define class
#
################################################################################
################################################################################
class elev:

################################################################################
# Parse input
################################################################################

	def parse_input(self):

		parser=argparse.ArgumentParser(description='Create elev file')

		# optional 
		parser.add_argument('-m','--outmissing',help='Missing value for output',required=False,default=-1e20,type=float)

		# positional
		parser.add_argument('infile',help='Input file')
		parser.add_argument('maskfile',help='Mask file')
		parser.add_argument('maskvar',help='Mask var')
		parser.add_argument('outfile',help='Output file')
		parser.add_argument('nxin',help='number of x points in input elev file',type=int)
		parser.add_argument('nyin',help='number of y points in input elev file',type=int)

		# Parse the arguments
		args=parser.parse_args()

		return args.infile, \
				args.maskfile, \
				args.maskvar, \
				args.outfile, \
				args.outmissing, \
				args.nxin, \
				args.nyin

################################################################################
################################################################################
# 
# Main function 
#
################################################################################
################################################################################


if __name__=='__main__':

	# Call the class
	e=elev()

	infile, maskfile, maskvar, outfile, outmv, nxin, nyin = e.parse_input()

	mf=nc.Dataset(maskfile,'r')
	maskdata=mf.variables[maskvar][:]
	maskindx=np.where(maskdata>0)
	eastings=mf.variables['x'][:]
	northings=mf.variables['y'][:]

	nxout=len(eastings)
	nyout=len(northings)

	inf=np.fromfile(infile,dtype='float32')
	inf=inf.byteswap().reshape((nyin,nxin))
	inf=inf[::-1,:]

	elev=np.ma.masked_equal(np.ones_like(inf)*outmv,outmv)
	elev[np.where(maskdata>0)]=inf[np.where(maskdata>0)]

	outf=nc.Dataset(outfile,'w')

	outf.createDimension('x',nxout)
	outf.createDimension('y',nyout)

	outf.createVariable('elev','f',('y','x'),fill_value=outmv)


	outf.variables['elev'][:]=elev[:nyout,:nxout]
	outf.variables['elev'].units='Elevation (m)'

	outf.title='Elevation data copied from %s'%infile
	outf.close()



