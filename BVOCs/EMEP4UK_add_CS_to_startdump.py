#!/usr/bin/env python
# routine to remap WFD-data onto 2D grid
#
#
import netCDF4 as nc
#import matplotlib.pyplot as plt
import sys
import numpy as np
#import getpass
#import datetime as dt

infile=sys.argv[1]
outfile=sys.argv[2]
CSfile=sys.argv[3]

infile='/users/eow/edwcom/EMEP/EMEP4UK/JULES_startdump/EMEP4UK_VG.dump.spin3.20010101.0.nc'
outfile='/users/eow/edwcom/EMEP/EMEP4UK/JULES_startdump/EMEP4UK_VG.dump.HWSD-CS.nc'

CSfile='/users/eow/edwcom/EMEP/hwsd2emep/EMEP4UK_HWSD_SC.nc'


inf=nc.Dataset(infile,'r')
outf=nc.Dataset(outfile,'w')

CSf=nc.Dataset(CSfile,'r')

inlats=inf.variables['latitude'][:]
inlons=inf.variables['longitude'][:]

CSlats=CSf.variables['lat'][:].squeeze().flatten()
CSlons=CSf.variables['lon'][:].squeeze().flatten()
CSdata2D=CSf.variables['cs'][:].squeeze().flatten()

#CSlons=np.round(CSlons,2)
#CSlats=np.round(CSlats,2)

CSlons[CSlons>180]-=360.

# loop round each
CSdata1D=np.zeros_like(inlats)

for i,lat,lon in zip(range(len(inlats)),inlats,inlons):
    index=np.argmin( (CSlons-lon)**2 + (CSlats-lat)**2 )
    #index=np.where( (np.round(CSlons,3)==np.round(lon,3)) & \
    #                (np.round(CSlats,3)==np.round(lat,3)) )
    #index=np.where((CSlons==np.round(lon,2))&(CSlats==np.round(lat,2)))
    CSdata1D[i]=CSdata2D[index]


for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    #for ncattr in inf.variables[str(var)].ncattrs():
    #    outvar.setncattr(str(ncattr),inf.variables[str(var)].getncattr(str(ncattr)))
    #
    if not (str(var) in ['cs']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        print 'here'
        outvar[:]=CSdata1D

outf.close()
inf.close()















