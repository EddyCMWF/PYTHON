#!/usr/bin/python2.7
#
# Program: WFDEI_rescaletoCS.py
# Purpose: rescale JULES CH4 output
# 
# Edward Comyn-Platt, 2015
#
######################################################
import numpy as np
import netCDF4 as nc
import plot_tools as PT

cs_file = '/users/eow/edwcom/WFD_EI/soil_cs_nofill.nc'
cs_ncname='carbMassMod'

JULES_infile='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/JULES_WFDEI_nti_NG-HWSD.monthly_wetl.nc'

JULES_outfile='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/JULES_WFDEI_nti_NG-HWSD_noCSfill.monthly_wetl.nc'

grid_file = '/users/eow/edwcom/WFD_EI/EI-Halfdeg-land-elevation.nc'

# read in lat/lon
grinf=nc.Dataset(grid_file,'r')
grlat=grinf.variables['latitude'][:]
grlon=grinf.variables['longitude'][:]
grinf.close()


inf= nc.Dataset(cs_file,'r')
new_SC = inf.variables[cs_ncname][:]
sc_lon=inf.variables['lon'][:]
sc_lat=inf.variables['lat'][:]
inf.close()


# loop round each lat lon pair
new_SC_1D=np.zeros_like(grlat)
for i,lat,lon in zip(range(len(grlat)),grlat,grlon):
    x_index=np.where(sc_lon==lon)[0]
    y_index=np.where(sc_lat==lat)[0]
    new_SC_1D[i]=new_SC[y_index,x_index]

mask = np.ones_like(new_SC_1D)
mask[np.isnan(new_SC_1D)] = 0.0

inf=nc.Dataset(JULES_infile,'r')
outf=nc.Dataset(JULES_outfile,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    for att in inf.variables[var].ncattrs():
        outvar.setncattr(str(att), \
                         inf.variables[var].getncattr(att) )

    if not (str(var) in ['fch4_wetl']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        print 'here'
        fch4_wetl = inf.variables[str(var)][:]
        fch4_wetl = fch4_wetl*mask
        outvar[:]=fch4_wetl

for att in inf.ncattrs():
    outf.setncattr(att,inf.getncattr(att))

outf.close()
inf.close()

