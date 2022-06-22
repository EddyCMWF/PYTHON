#!/bin/env python2.7

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

TRENDY_DIR = '/prj/CLIFFTOP/TRENDY/'
latlonfile=TRENDY_DIR+'S3/CRU-NCEP_1p1.update.vn4.6.S3.Annual.cVeg.nc'
area_file=TRENDY_DIR+'Area.nc'

inf=nc.Dataset(latlonfile,'r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]

lons_2d,lats_2d=np.meshgrid(lons,lats)

# calculate the area of the gid
rearth=6371000.
latdel = 1.25
londel = 1.875
area=rearth*rearth*np.radians(londel)*(np.sin(np.radians(lats_2d+latdel))-np.sin(np.radians(lats_2d)))

plt.imshow(area,origin='bottom')
plt.colorbar()
plt.show()

outf=nc.Dataset(area_file,'w')
outf.createDimension('latitude',len(inf.dimensions['latitude']))
outf.createDimension('longitude',len(inf.dimensions['longitude']))
outf.createDimension('time',1)

for var in ['latitude','longitude']:
    outvar=outf.createVariable(var,'float32',inf.variables[var].dimensions)
    for att in inf.variables[var].ncattrs():
        outvar.setncattr(att,inf.variables[var].getncattr(att))
    if var in ['latitude','longitude']:
        outvar[:]=inf.variables[var][:]
    elif var in ['time']:
        outvar[:]=inf.variables[var][0]
    else:
        outvar[:]=inf.variables[var][0,:]


outvar=outf.createVariable('area','float32',('time','latitude','longitude'),fill_value=-1e20)
outvar.long_name='Grid Cell Area'
outvar.coordinates='latitude longitude'
outvar.units='m2'
outvar[:]=area


inf.close()

