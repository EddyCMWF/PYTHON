#!/usr/bin/env python
##################################################################################
#
# Program: fix_n96_land_frac.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: fix the n96 land frac file to double precision for JULES run
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import pylab as plt

infile='/prj/jules_water/data/cru_ncep/n96/ancil/trendy_2013/HadGEM2ES_Ancil.nc'
soilfile='/users/eow/edwcom/CRUNCEP/n96/ancil/HadGEM2ES_Soil_Ancil.nc'
outfile='/users/eow/edwcom/CRUNCEP/n96/ancil/HadGEM2ES_Ancil_FixedFrac.nc'

param_name='field1391'
lsm_name='lsm'

inf=nc.Dataset(infile,'r')
outf=nc.Dataset(outfile,'w')

frac_in=inf.variables[param_name][:]
lsm_in=inf.variables[lsm_name][:]

#create LSmask
LSmask=np.ones_like(lsm_in.data)
LSmask[np.where((lsm_in.data>1.0) | (lsm_in.data<=0.))]=0.
LSmask_inv=np.ones_like(LSmask)-LSmask

frac_tot=np.sum(frac_in.data,axis=0,dtype='float64')

frac_out=np.zeros_like(frac_in.data,dtype='float64')+frac_in.data
for tile in range(frac_out.shape[0]):
    frac_out[tile,:,:]=frac_out[tile,:,:]/frac_tot
    frac_out[tile,:,:]=frac_out[tile,:,:]*LSmask

# apply sea mask
new_frac_tot=np.sum(frac_out,axis=0,dtype='float64')

icemask=np.ones_like(new_frac_tot)

icemask[np.where( (frac_out[8,:,:]>0.) & \
                  (frac_out[8,:,:]<=1.) )]=0

ice=np.ones_like(icemask)
ice=ice-icemask

for tile in range(8):
    frac_out[tile,:,:]=frac_out[tile,:,:]*icemask

frac_out[8,:,:]=ice

#frac_out[7,:,:]=frac_out[7,:,:]+LSmask_inv

#frac_out[np.where(frac_in.mask==True)]=.0
#frac_out=np.ma.masked_equal(frac_out,-9999.0)
new_frac_tot=np.sum(frac_out,axis=0,dtype='float64')

outf.createDimension('pseudo',frac_out.shape[0])
outf.createDimension('latitude',frac_out.shape[1])
outf.createDimension('longitude',frac_out.shape[2])
outvar=outf.createVariable(param_name,'float64',('pseudo','latitude','longitude'))
for att in inf.variables[param_name].ncattrs():
    if (str(att)!='_FillValue'):
        outvar.setncattr(str(att),inf.variables[param_name].getncattr(att))

outvar[:]=frac_out

outvar=outf.createVariable('land_frac','float64',('latitude','longitude'))
outvar[:]=frac_out[7,:,:]

outf.note='Land Fractions fixed, now they sum to 1'
outf.author='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.close()

soilf=nc.Dataset(soilfile,'a')
sm_sat=soilf.variables['field332']
temp=sm_sat[:]*icemask
sm_sat[:]=temp
soilf.close()

