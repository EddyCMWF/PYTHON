#!/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib
import numpy as np
#from scipy import stats
import ipdb

SWAMPS_dir='/prj/CLIFFTOP/SWAMPS/'
SWAMPS_infile  = SWAMPS_dir+'fw_swamps-glwd_2000-2012.nc'
SWAMPS_outfile = SWAMPS_dir+'fw_swamps-glwd_2000-2012_climatology.nc'


# Open in and out files:
SWAMPSinf  = nc.Dataset(SWAMPS_infile,'r')
SWAMPSoutf = nc.Dataset(SWAMPS_outfile,'w')

for dim in SWAMPSinf.dimensions:
    indim=SWAMPSinf.dimensions[dim]
    if str(dim) not in ['time']:
        SWAMPSoutf.createDimension(str(dim),len(indim))
    else:
        SWAMPSoutf.createDimension(str(dim),12)

for var in ['lat', 'lon', 'time']:
    invar=SWAMPSinf.variables[var]
    outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
    for att in invar.ncattrs():
        outvar.setncattr(str(att),invar.getncattr(att))
    if 'time' not in invar.dimensions:
        outvar[:] = invar[:]
    else:
        outvar[:] = invar[:12]

#ipdb.set_trace()
nlats = len(SWAMPSinf.dimensions['lat'])
nlons = len(SWAMPSinf.dimensions['lon'])
nmonths = 12
nyears = len(SWAMPSinf.dimensions['time'])/12

invar=SWAMPSinf.variables['Fw']

in_FW_data = invar[:].reshape(nyears,nmonths,nlats,nlons) #13,12,360,720)

var = 'Fw_mean'
outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
for att in invar.ncattrs():
    outvar.setncattr(str(att),invar.getncattr(att))
outvar.setncattr('long_name',"Fraction inundated, mean climatology")
outdata = np.mean(in_FW_data,axis=0) 
outvar[:] = outdata

var = 'Fw_median'
outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
for att in invar.ncattrs():
    outvar.setncattr(str(att),invar.getncattr(att))
outvar.setncattr('long_name',"Fraction inundated, median climatology")
outdata = np.median(in_FW_data,axis=0) 
outvar[:] = outdata

var = 'Fw_max'
outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
for att in invar.ncattrs():
    outvar.setncattr(str(att),invar.getncattr(att))
outvar.setncattr('long_name',"Fraction inundated, max climatology")
outdata = np.max(in_FW_data,axis=0) 
outvar[:] = outdata

var = 'Fw_min'
outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
for att in invar.ncattrs():
    outvar.setncattr(str(att),invar.getncattr(att))
outvar.setncattr('long_name',"Fraction inundated, min climatology")
outdata = np.min(in_FW_data,axis=0) 
outvar[:] = outdata

var = 'Fw_mask'
outvar=SWAMPSoutf.createVariable(str(var),invar.dtype,invar.dimensions)
for att in invar.ncattrs():
    outvar.setncattr(str(att),invar.getncattr(att))
outvar.setncattr('long_name',"Inundated mask, min Fw > 0.01")
outdata[outdata>0.01] = 1
outvar[:] = outdata


SWAMPSinf.close()
SWAMPSoutf.close()

