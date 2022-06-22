#!/usr/bin/python

import netCDF4 as nc
import numpy as np

CRUNCEP_DIR   = '/users/eow/edwcom/CRUNCEP/'

CRUNCEP_indump = '/prj/ALANIS/jules_data/start_dumps/jasmin_global/NCEP_CRU/m045_GLOBAL_20110101_000000_dump_new.nc'

CRUNCEP_outdump = CRUNCEP_DIR+'CRUNCEP_startdump_from_WFDEIdump.nc'

WFDEI_DIR     = '/users/eow/edwcom/WFD_EI/'
WFDEI_gfile   = WFDEI_DIR+'EI-Halfdeg-land-elevation.nc'
WFDEI_indump  = WFDEI_DIR+'start_dumps/WFDEI_startdump_withSC_ECP20150731.nc'

WFDEI_index_file = '/users/eow/edwcom/CRUNCEP/WFDEI_to_CRUNCEP_index.dat'

WFDEI_SC_tempname='field1397'

inf=open(WFDEI_index_file,'r')
index=np.array(inf.readlines(),dtype='int32')
inf.close()

Cinf  = nc.Dataset(CRUNCEP_indump,'r')
Clats = Cinf.variables['lat'][:]
Clons = Cinf.variables['lon'][:]
Cinf.close()
Clatlon = zip(Clats,Clons) ## [ (clat,clon) for clat,clon in zip(Clats,Clons) ]

Winf=nc.Dataset(WFDEI_gfile,'r')
Wlats = Winf.variables['latitude'][:]
Wlons = Winf.variables['longitude'][:]
Winf.close()
Wlatlon = zip(Wlats,Wlons) ## [ (wlat,wlon) for wlat,wlon in zip(Wlats,Wlons) ]

# Replace missing index values with closest neighbour
for ind in np.where(index==-9999.0)[0]:
    distarray = ((Wlats-Clats[ind])**2) + ((Wlons-Clons[ind])**2)
    index[ind] = np.argmin( distarray )


# Open the WFDEI dump:
Winf=nc.Dataset(WFDEI_indump,'r')

# OPen the CRUNCEP out dump file.
Coutf=nc.Dataset(CRUNCEP_outdump,'w')

# Copy dimensions:
for dim in Winf.dimensions:
    if (dim!='land'):
        Coutf.createDimension( str(dim),len(Winf.dimensions[dim]) )
    else:
        Coutf.createDimension( str(dim),len(index) )

for var in Winf.variables:
    outvar= Coutf.createVariable( str(var), \
                                  Winf.variables[var].dtype, \
                                  Winf.variables[var].dimensions )

    if (len(Winf.variables[str(var)][:].shape)==1):
        outdata=Winf.variables[str(var)][index]
    elif (len(Winf.variables[str(var)][:].shape)==2):
        outdata=Winf.variables[str(var)][:,index]
    elif (len(Winf.variables[str(var)][:].shape)==3):
        outdata=Winf.variables[str(var)][:,:,index]
        
    outvar[:]=outdata

Coutf.title='CRUNCEP start dump reconstructed from WFDEI start dump'
Coutf.note='Gaps filled with closest neighbour value'
Coutf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
Coutf.history='Created 30/09/2015'

Coutf.close()


Winf.close()

