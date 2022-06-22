#!/usr/bin/python
#
# Python module to plot the hwsd dat on EMEP grid
#
# Edard Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plot


INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/EMEP4UK_EUROPE_gridfile.nc'

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/EMEP4UK_gridfile.nc'

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/hwsd_30arcsec_rebinned.nc'

INFILE = '/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/hwsd_30arcsec_EMEP4UK.nc'

INFILE='/users/eow/edwcom/EMEP/topo2emep/x1y1/input_data/grid_files/EMEP4UK_EUROPE_gridfile.nc'

INFILE='/users/eow/edwcom/EMEP/topo2emep/x1y1/input_data/topidx.nc'

INFILE='/users/eow/edwcom/EMEP/chess2emep/input_data/grid_files/chess_1km.nc'

INFILE='/users/eow/edwcom/EMEP/chess2emep/input_data/grid_files/EMEP4UK_gridfile.nc'


inf=nc.Dataset(INFILE,'r')

dims = inf.variables['grid_dims'][:]
mask = inf.variables['grid_imask'][:]
lats = inf.variables['grid_center_lat'][:]
lons = inf.variables['grid_center_lon'][:]
latcors = inf.variables['grid_corner_lat'][:]
loncors = inf.variables['grid_corner_lon'][:]
inf.close()

plot.figure(1,figsize=(20,10))
plot.subplot(1,3,1)
plot.imshow(lats.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.subplot(1,3,2)
plot.imshow(lons.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.subplot(1,3,3)
plot.imshow(mask.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.show()

vmin,vmax=-0.01,0.01
vmin,vmax=-0.5,0.5

plot.figure(1,figsize=(25,10))
plot.subplot(3,4,1)
plot.imshow(latcors[:,0].reshape(dims[::-1])-lats.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,2)
plot.imshow(latcors[:,1].reshape(dims[::-1])-lats.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,3)
plot.imshow(latcors[:,2].reshape(dims[::-1])-lats.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,4)
plot.imshow(latcors[:,3].reshape(dims[::-1])-lats.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,5)
plot.imshow(loncors[:,0].reshape(dims[::-1])-lons.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,6)
plot.imshow(loncors[:,1].reshape(dims[::-1])-lons.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,7)
plot.imshow(loncors[:,2].reshape(dims[::-1])-lons.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,8)
plot.imshow(loncors[:,3].reshape(dims[::-1])-lons.reshape(dims[::-1]),origin='bottom',vmin=vmin,vmax=vmax)
plot.colorbar(fraction=0.04)
plot.subplot(3,4,9)
plot.imshow(mask.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.subplot(3,4,10)
plot.imshow(lats.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.subplot(3,4,11)
plot.imshow(lons.reshape(dims[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)


plot.show()



log_file='/users/eow/edwcom/code/shell_scripts/EMEP/scrip2.log'


log_file='/users/eow/edwcom/EMEP/chess2emep/scrip.log'
l_file=open(log_file)
l_data=l_file.readlines()
l_file.close()
l_data=l_data[6:-2]

index1=[]
index2=[]
extremity=[]
mapnum=[]
for line in l_data:
    index1.append(line.split()[5])
    index2.append(line.split()[6])
    extremity.append(line.split()[-1])
    mapnum.append(line.split()[1])

mapnum=np.array(mapnum,dtype='int')
extremity=np.array(extremity,dtype='float32')
index1=np.array(index1,dtype='int64')
index2=np.array(index2,dtype='int64')

index1=index1[mapnum==2]
index2=index2[mapnum==1]

index2=np.sort( np.array( list(set(index2)), dtype='int' ) )
index1=np.sort( np.array( list(set(index1)), dtype='int' ) )


#index2=[]
#for line in l_data:
#    index2.append(line.split()[6])
#
#index2=np.sort( np.array( list(set(index2)), dtype='int' ) )


for pt in index1:
    print lats[pt], latcors[pt,:]
    print lons[pt], loncors[pt,:]

mask1=mask_BU.copy()

mask1[index1]=-1
plot.imshow(mask1.reshape(dims1[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.show()

mask2=mask.copy()
mask2[index2]=-1
plot.imshow(mask2.reshape(dims2[::-1]),origin='bottom')
plot.colorbar(fraction=0.04)
plot.show()


# the histogram of the data with histtype='step'
n, bins, patches = plot.hist(extremity, 50, normed=1, histtype='stepfilled')
plot.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

## add a line showing the expected distribution
#y = P.normpdf( bins, mu, sigma)
#l = P.plot(bins, y, 'k--', linewidth=1.5)

