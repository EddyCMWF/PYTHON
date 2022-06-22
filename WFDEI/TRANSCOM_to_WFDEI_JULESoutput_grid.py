#!/bin/python

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

transcom_infile='/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5.nc'
transcom_outfile='/users/eow/edwcom/WFD_EI/TRANSCOM_Regions_WFDEI_landpoints.nc'

WFDEI_lmask_file='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
grindex = nc.Dataset(WFDEI_lmask_file).variables['land_index'][:]-1

grid_file = '/users/eow/edwcom/WFD_EI/EI-Halfdeg-land-elevation.nc'
# read in land_index
grinf=nc.Dataset(grid_file,'r')
grlat=grinf.variables['latitude'][:]
grlon=grinf.variables['longitude'][:]
grinf.close()


inf=nc.Dataset(transcom_infile,'r')
transcom_2d = inf.variables['transcom_regions'][:]
tr_lat=inf.variables['latitude'][:]
tr_lon=inf.variables['longitude'][:]
inf.close()


#TRAindex = range(360,720)
#for i in range(360):
#    TRAindex.append(i)
#transcom_2d = transcom_2d[:,TRAindex]

tr_lon[tr_lon>180.]-=360.
#tr_lon=tr_lon[TRAindex]

# loop round each
transcom_1d=np.zeros_like(grlat)
for i,lat,lon in zip(range(len(grlat)),grlat,grlon):
    x_index=np.where(tr_lon==lon)[0]
    y_index=np.where(tr_lat==lat)[0]
    transcom_1d[i]=transcom_2d[y_index,x_index]


# Check it converts back okay
plt.imshow(transcom_1d[grindex],origin='bottom')
plt.show()


outf=nc.Dataset(transcom_outfile,'w')

outf.createDimension('land',len(index_1d))

outvar=outf.createVariable('transcom_regions','float32',('land'))
outvar.missing_value=-999.9
outvar[:]=transcom_1d
outf.close()


