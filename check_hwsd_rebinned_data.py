#!/usr/bin/python
#
# Python module to plot the hwsd dat on EMEP grid
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plot
import plot_map_ECP as PM
#
# Input arguments:
if (len(sys.argv)>1):
    INFILE         = sys.argv[1]
else:
    INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_rebin_emep_50km_grid_missing.nc'
INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_rebinned.nc'
INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd.nc'

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_rebin_emep_50km_grid.nc'

INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_rebin_emep_50km_grid_missing.nc'
INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsdREBIN_emep4ukEUROPE_grid_missing.nc'

INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsdREBIN_emep4ukEUROPE_grid.nc'


INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsdREBIN_emep4uk_grid_missing.nc'

INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_emep4uk_grid_missing.nc'

INFILE         = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/hwsd_emep4uk_grid.nc'

INFILE='/users/eow/edwcom/EMEP/hwsd2emep/EMEP4UK_soilparams_hwsd_bc.nc'

inf=nc.Dataset(INFILE,'r')

lons = inf.variables['lon'][:]
lats = inf.variables['lat'][:]
mu   = inf.variables['mu'][:]

mask = np.ones_like(mu.data)
mask[np.where(mu.mask)]=0.0

inf.close()

#mu.shape=(len(mu[0,:,0]),len(mu[0,0,:]))
# Normal orientaiton:
#temp = (mu.reshape( len(mu[0,:,0]),len(mu[0,0,:]) ) ).transpose(1,0)[:,::-1]

plot.imshow( mu.reshape( len(mu[0,:,0]),len(mu[0,0,:]) ) , origin='bottom')
plot.colorbar()
plot.show()

temp=mu.data.reshape( len(mu[0,:,0]),len(mu[0,0,:]) )
plot.imshow(temp,origin='bottom')
plot.colorbar()
plot.show()

plot.imshow( mask.reshape( len(mu[0,:,0]),len(mu[0,0,:]) ) , origin='bottom' )
plot.show()

PM.plot_map(data,lons,lats,DATA_RANGE=[0,25000],NLEVELS=200,MPL_CBAR='spectral')

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/EMEP4UK_EUROPE_gridfile.nc'
#INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/hwsd_30arcsec_rebinned.nc'
inf_EMEP=nc.Dataset(INFILE,'r')

dims = inf_EMEP.variables['grid_dims'][:]
mask = inf_EMEP.variables['grid_imask'][:]
lats = inf_EMEP.variables['grid_center_lat'][:]
lons = inf_EMEP.variables['grid_center_lon'][:]
inf.close()

plot.imshow(mask.reshape(dims[::-1]),origin='bottom')
plot.show()

plot.imshow(mask.reshape(dims[::-1])[700:1100,4000:4500])
#plot.imshow(lats.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.show()



INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/EMEP4UK_EUROPE_soilparams_hwsd_bc.nc'

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/EMEP4UK_EUROPE_soilparams_hwsd_vg.nc'


inf=nc.Dataset(INFILE,'r')

lons = inf.variables['lon'][:]
lats = inf.variables['lat'][:]

hcap = inf.variables['hcap'][:]

plot



sathh=inf.variables['sathh'][:]
b=inf.variables['b'][:]
hcap=inf.variables['hcap'][:]
hcon=inf.variables['hcon'][:]
satcon=inf.variables['satcon'][:]
vcrit=inf.variables['vcrit'][:]
vsat=inf.variables['vsat'][:]
vwilt=inf.variables['vwilt'][:]
cs=inf.variables['cs'][:]

plot.figure(1,figsize=(15,15))

plot.subplot(3,3,1)
plot.subplot(3,3,2)
plot.subplot(3,3,3)


plot.imshow(b[0,:,:].reshape(b.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('bexp')
plot.colorbar()
plot.show()

plot.figure(1,figsize=(10,10))
plot.imshow(cs[0,:,:].reshape(cs.shape[1:]),origin='bottom',cmap='YlGn',vmax=35,vmin=0)
plot.title('cs')
plot.colorbar()
plot.show()

plot.figure(1,figsize=(10,10))
plot.imshow(hcap[0,:,:].reshape(hcap.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('hcap')
plot.colorbar()
plot.show()


plot.figure(1,figsize=(10,10))
plot.imshow(hcon[0,:,:].reshape(hcon.shape[1:]),origin='bottom',cmap='coolwarm',vmin=0,vmax=0.35)
plot.title('hcon')
plot.colorbar()
plot.show()

plot.subplot(3,3,4)
plot.imshow(satcon[0,:,:].reshape(satcon.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('satcon')
plot.colorbar()

plot.subplot(3,3,5)
plot.imshow(sathh[0,:,:].reshape(sathh.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('sathh')
plot.colorbar()

plot.subplot(3,3,6)
plot.imshow(vcrit[0,:,:].reshape(vcrit.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('vcrit')
plot.colorbar()

plot.subplot(3,3,7)
plot.imshow(vsat[0,:,:].reshape(vsat.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('vsat')
plot.colorbar()

plot.subplot(3,3,8)
plot.imshow(vwilt[0,:,:].reshape(vwilt.shape[1:]),origin='bottom',cmap='coolwarm')
plot.title('vwilt')
plot.colorbar()


plot.show()
