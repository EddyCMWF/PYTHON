#!/bin/python

import netCDF4 as nc
import numpy as np

file='/users/eow/edwcom/CHESS/CHESS_constant_Ndeposition.nc'

nx=656
ny=1057
ntime=1

outf=nc.Dataset(file,'w')

outf.createDimension('x',nx)
outf.createDimension('y',ny)
outf.createDimension('time',ntime)

outvar=outf.createVariable('deposition_N','float32',('time','y','x'),fill_value=-99999.)
outvar.units='kg m-2 s-1'
outvar[:] = np.ones([ntime,ny,nx])*3.169e-11

outf.close()

