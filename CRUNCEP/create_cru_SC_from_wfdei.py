#!/usr/bin/python

import netCDF4 as nc
import numpy as np

CRUNCEP_DIR   = '/users/eow/edwcom/CRUNCEP/'
CRUNCEP_gfile = 'cru_ncep_land.nc'
CRUNCEP_HWSD_SC_outfile= 'CRUNCEP_HWSD_soilcarbon_landpoints_0p5.nc'

WFDEI_DIR     = '/users/eow/edwcom/WFD_EI/'
WFDEI_gfile   = 'EI-Halfdeg-land-elevation.nc'
WFDEI_HWSD_SCfile = 'soil_cs_nofill.nc'

WFDEI_index_file = '/users/eow/edwcom/CRUNCEP/WFDEI_to_CRUNCEP_index.dat'

WFDEI_SC_tempname='carbMassMod'

inf=open(WFDEI_index_file,'r')
index=np.array(inf.readlines(),dtype='int32')
inf.close()

Cinf  = nc.Dataset(CRUNCEP_DIR+CRUNCEP_gfile,'r')
Clats = Cinf.variables['Latitude'][:]
Clons = Cinf.variables['Longitude'][:]
Cinf.close()
Clatlon = zip(Clats,Clons) ## [ (clat,clon) for clat,clon in zip(Clats,Clons) ]

# Read in the SC data:
inf=nc.Dataset(WFDEI_DIR+WFDEI_HWSD_SCfile,'r')
new_SC = inf.variables[WFDEI_SC_tempname][:].squeeze()
sc_lon = inf.variables['lon'][:]
sc_lat = inf.variables['lat'][:]
inf.close()

# loop round each lat lon pair
new_SC_1D=np.zeros_like(Clats)+np.nan
for i,lat,lon in zip(range(len(Clats)),Clats,Clons):
    x_index=np.where(sc_lon==lon)[0]
    y_index=np.where(sc_lat==lat)[0]
    new_SC_1D[i]=new_SC[y_index,x_index]

badex   = np.where(np.isnan(new_SC_1D))[0]
#new_SC_1D[badex]=0.0
goodex = np.where(np.isfinite(new_SC_1D))[0]
print badex.shape
for pt in badex:
    temppoint = np.argmin(  ((Clats[pt]-Clats[goodex])**2) \
                          + ((Clons[pt]-Clons[goodex])**2) )
    new_point=goodex[temppoint]
    new_SC_1D[pt]=new_SC_1D[new_point]

print badex.shape

outf=nc.Dataset(CRUNCEP_DIR+CRUNCEP_HWSD_SC_outfile,'w')

outf.createDimension('land',len(new_SC_1D))

outvar=outf.createVariable('CS_HWSD','float32',('land'))
outvar.units='kgC m^-2'
outvar.longname='Soil Carbon from HWSD'
outvar[:]=new_SC_1D

outvar=outf.createVariable('lats','float32',('land'))
outvar.units=''
outvar.longname='Latitude'
outvar[:]=Clats

outvar=outf.createVariable('lons','float32',('land'))
outvar.units=''
outvar.longname='Longitude'
outvar[:]=Clons

outf.title='Soil Carbon from HWSD and HWSD-NCSCD for the NCEP-CRU land points'
outf.note='Data extracted from WFDEI data, gaps filled with closest neighbour value'
outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.history='Created 09/09/2015'

outf.close()




