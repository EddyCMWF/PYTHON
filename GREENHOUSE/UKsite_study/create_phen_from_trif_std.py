#!/bin/env python
import netCDF4 as nc
import os

npfts=5
lai=[5,4,2,4,2]
canht=[19.01,16.38,0.79,1.26,1.00]
nscpools=1
cs=[8.0]
clay=0.23

Phenolfile='Brattleby_phen.dump.20130601.0.nc'
TRIFFIDfile='Brattleby_std.dump.20130601.0.nc'


# Reopen phenol as readonly
Pf=nc.Dataset(Phenolfile,'w')
# Open TRIFFID file in write mode
Tf=nc.Dataset(TRIFFIDfile,'r')
# Copy dimensions over, extending scpool dimension to specified amount
for dim in Tf.dimensions:
    if str(dim)!='scpool':
        Pf.createDimension(str(dim),len(Tf.dimensions[dim]))
    else:
        Pf.createDimension(str(dim),nscpools)

# Copy variables over, inserting new cs values
for var in Tf.variables:
    outvar=Pf.createVariable(str(var),'float32',Tf.variables[var].dimensions)
    if str(var)!='cs':
        outvar[:]=Tf.variables[var][:]
    else:
        outvar[:]=cs

Tf.close()
Pf.close()
