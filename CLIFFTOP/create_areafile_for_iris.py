#!/bin/env python2.7

import netCDF4 as nc
import numpy as np

jules_file='/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/CEN_BCC_MOD_bcc-csm1-1/vn4.8_imogen_CEN_BCC_MOD_bcc-csm1-1_2deg.Annual_carbon.1860.nc'
area_file='/group_workspaces/jasmin2/clifftop/COMMON_DATA/ANCILS/Area_in_iris_format.nc'

inf=nc.Dataset(jules_file,'r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]

# calculate the area of the gid
rearth=6371000.

area=rearth*rearth*np.radians(3.75)*(np.sin(np.radians(lats+2.5))-np.sin(np.radians(lats)))

outf=nc.Dataset(area_file,'w')
outf.createDimension('y',len(inf.dimensions['y']))
outf.createDimension('x',len(inf.dimensions['x']))
outf.createDimension('time',1)
outf.createDimension('nt',2)

for var in ['time_bounds','time','latitude','longitude']:
    outvar=outf.createVariable(var,'float32',inf.variables[var].dimensions)
    for att in inf.variables[var].ncattrs():
        outvar.setncattr(att,inf.variables[var].getncattr(att))
    if var in ['latitude','longitude']:
        outvar[:]=inf.variables[var][:]
    elif var in ['time']:
        outvar[:]=inf.variables[var][0]
    else:
        outvar[:]=inf.variables[var][0,:]


outvar=outf.createVariable('area','float32',('time','y','x'),fill_value=-1e20)
outvar.long_name='Grid Cell Area'
outvar.coordinates='latitude longitude'
outvar.units='m2'
outvar[:]=area


inf.close()

