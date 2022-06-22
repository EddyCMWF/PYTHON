#!/usr/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse

class pressure:

	def __init__(self):

		self.lapse_rate = -6.5e-3 # K m-1

		self.mol_mass_air = 0.0289644 # kg mol-1

		self.r_air = 8.31432 # N m mol-1 K-1

		self.g = 9.81 # m s-2

################################################################################
# Parse input
################################################################################

	def parse_input(self):

		parser=argparse.ArgumentParser(description='Create pressure at new altitude')

		# optional
		parser.add_argument('--inpsurfvar',help='Input surface pressure variable', required=False, default='psurf')
		parser.add_argument('--intairvar',help='Input variable', required=False, default='tair')
		parser.add_argument('--inaltvar',help='Input altitude variable', required=False, default='z')
		parser.add_argument('--newaltvar',help='Input new altitude variable', required=False, default='z')
		parser.add_argument('--outmv',help='Missing value in output', required=False, default='-1e20',type=float)

		# positional
		parser.add_argument('inpsurf',help='Input surface pressure file')
		parser.add_argument('intair',help='Input air temperature file')
		parser.add_argument('inalt',help='Input altitude file')
		parser.add_argument('newalt',help='Input new altitude file')
		parser.add_argument('outpsurf',help='Output surface pressure file')
		parser.add_argument('outtair',help='Output air temperature file')

		# Parse the arguments
		args=parser.parse_args()

		return args.inpsurf, \
				args.intair, \
				args.inalt, \
				args.newalt, \
				args.outpsurf, \
				args.outtair, \
				args.inpsurfvar, \
				args.intairvar, \
				args.inaltvar, \
				args.newaltvar, \
				args.outmv


################################################################################
################################################################################
#
# Start the main routine
#
################################################################################
################################################################################

if __name__=='__main__':

	p=pressure()

	inpsurf, intair, inalt, newalt, outpsurf, outtair, inpsurfvar, intairvar, inaltvar, newaltvar, outmv = p.parse_input()

	inps=nc.Dataset(inpsurf,'r')
	psurf_1 = inps.variables[inpsurfvar][:]
	nt,ny,nx=psurf_1.shape

	inta=nc.Dataset(intair,'r')
	tair_1 = inta.variables[intairvar][:]

	inz=nc.Dataset(inalt,'r')
	alt_1 = inz.variables[inaltvar][0,:]

	newz=nc.Dataset(newalt,'r')
	alt_2 = newz.variables[newaltvar][:]

	z = alt_2 - alt_1

	tair_2 = tair_1 + p.lapse_rate*z

	psurf_2 = np.ma.masked_where(psurf_1.mask,np.ones_like(psurf_1)*outmv)
	
	indx=np.where(~psurf_1.mask)
	psurf_2[indx] = psurf_1[indx] * ( (tair_1[indx]/tair_2[indx])**(p.g*p.mol_mass_air/(p.r_air*p.lapse_rate)) )

	op = nc.Dataset(outpsurf,'w')

	op.createDimension('x',nx)
	op.createDimension('y',ny)
	op.createDimension('time',nt)

	op.createVariable('time',inps.variables['time'].dtype,inps.variables['time'].dimensions)
	op.variables['time'][:]=inps.variables['time'][:]
	for attr in inps.variables['time'].ncattrs():
		op.variables['time'].setncattr(attr,inps.variables['time'].getncattr(attr))

	#op.createVariable('x',float,'x')
	#op.variables['x'][:]=np.arange(0,1000*nx,1000)
	#op.variables['x'].long_name='eastings'

	#op.createVariable('y',float,'y')
	#op.variables['y'][:]=np.arange(0,1000*ny,1000)
	#op.variables['y'].long_name='northings'

	for var in ['lat','lon']:
		op.createVariable(var,'f',('y','x'))
		op.variables[var][:]=inps.variables[var][:] 
		for attr in inps.variables[var].ncattrs():
			op.variables[var].setncattr(attr,inps.variables[var].getncattr(attr))

	op.createVariable('psurf','f',('time','y','x'),fill_value=outmv)

	op.variables['psurf'][:]=psurf_2[:]

	op.variables['psurf'].units='Pa'
	op.variables['psurf'].long_name='Surface air pressure'

	op.note='Air pressure lapsed to CHESS altitudes'

	op.close()

	ot = nc.Dataset(outtair,'w')

	ot.createDimension('x',nx)
	ot.createDimension('y',ny)
	ot.createDimension('time',nt)

	ot.createVariable('time',inps.variables['time'].dtype,inps.variables['time'].dimensions)
	ot.variables['time'][:]=inps.variables['time'][:]
	for attr in inps.variables['time'].ncattrs():
		ot.variables['time'].setncattr(attr,inps.variables['time'].getncattr(attr))

	#ot.createVariable('x',float,'x')
	#ot.variables['x'][:]=np.arange(0,1000*nx,1000)
	#ot.variables['x'].long_name='eastings'

	#ot.createVariable('y',float,'y')
	#ot.variables['y'][:]=np.arange(0,1000*ny,1000)
	#ot.variables['y'].long_name='northings'

	for var in ['lat','lon']:
		ot.createVariable(var,'f',('y','x'))
		ot.variables[var][:]=inps.variables[var][:] 
		for attr in inps.variables[var].ncattrs():
			ot.variables[var].setncattr(attr,inps.variables[var].getncattr(attr))



	ot.createVariable('tair','f',('time','y','x'),fill_value=outmv)

	ot.variables['tair'][:]=tair_2[:]

	ot.variables['tair'].units='Pa'
	ot.variables['tair'].long_name='Air temperature'

	ot.note='Air pressure lapsed to CHESS altitudes'

	ot.close()


