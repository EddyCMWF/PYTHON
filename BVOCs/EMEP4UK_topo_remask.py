#!/usr/bin/python
#
########################################################################
#
# Program: EMEP50km_restitch.py
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: To restitch remapped EMEP data
#
########################################################################
#
import numpy as np
import netCDF4 as nc
import pylab as plt
import plot_map_ECP as PM
#
#
#
EMEP_dir='/users/eow/edwcom/EMEP/'
#
param_subdir='topo2emep/'
infile=EMEP_dir+param_subdir+'topidx_emep4uk_grid_raw.nc'
param_dataname='mu'
#
Gridfile=EMEP_dir+'EMEP4UK/EMEP4UK_gridfile.nc'
#
outfile=EMEP_dir+param_subdir+'topidx_emep4uk_grid.nc'
outmv  = -9999.0
land_fillvalue = 3.0
land_fillvalue_std = 2.0
#
# Read Gridfile
Grid_inf    = nc.Dataset(Gridfile,'r')
dims        = Grid_inf.variables['grid_dims'][:]
i_mast      = Grid_inf.variables['i'][:]
j_mast      = Grid_inf.variables['j'][:]
i_EMEP_mast = Grid_inf.variables['grid_EMEP_i_coord'][:].reshape(dims[::-1])
j_EMEP_mast = Grid_inf.variables['grid_EMEP_j_coord'][:].reshape(dims[::-1])
lat_mast    = Grid_inf.variables['grid_center_lat'][:].reshape(dims[::-1])
lon_mast    = Grid_inf.variables['grid_center_lon'][:].reshape(dims[::-1])
mask_mast   = Grid_inf.variables['grid_imask'][:].reshape(dims[::-1])
Grid_inf.close()
#
# Read remapped data file
inf=nc.Dataset(infile,'r')
topidx     = inf.variables[param_dataname][:]
topidx_std = inf.variables[param_dataname+'_std'][:]
inf.close()

topidx[mask_mast==0]     = outmv
topidx_std[mask_mast==0] = outmv

########################################################################################
# Fill bad land points with closest good land point value
# index bad land points
badland_x, badland_y=np.where((mask_mast==1) & (topidx.data<0.))
nbadpoints=len(badland_x)
# create boolean array of where there are good land points
uber_mask = (mask_mast == 1) &  \
            (topidx.data> 0.)
# loop through bad land points
k=1
while (nbadpoints>0):
    print k, nbadpoints
    for cnt in range(nbadpoints):
        x=badland_x[cnt]
        y=badland_y[cnt]
        # find closest good point(s) by searching ever increasing square of data
        radius=k               # "radius" of square in pixels
        #while not any( (uber_mask[max(x-radius,0):min(x+radius+1,dims[1]),          \
            #                          max(y-radius,0):min(y+radius+1,dims[0])].flat )   ) :
        #    radius=radius+1
    
        # extract mask for surrounding area upto radius
        mask_temp=uber_mask[max(x-radius,0):min(x+radius+1,dims[1]),  \
                            max(y-radius,0):min(y+radius+1,dims[0])]
        if any(mask_temp.flat):
            # extract lat, lon and topidx and topidx_std for surrounding area upto radius
            # and immedietely apply mask
            lat_temp=lat_mast[max(x-radius,0):min(x+radius+1,dims[1]),  \
                              max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            lon_temp=lon_mast[max(x-radius,0):min(x+radius+1,dims[1]),  \
                              max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            topidx_temp=topidx[max(x-radius,0):min(x+radius+1,dims[1]),  \
                               max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            topidx_std_temp=topidx_std[max(x-radius,0):min(x+radius+1,dims[1]),  \
                                       max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            
            # locate closest good land point using lats and lons and fill gap with it
            topidx[x,y] = topidx_temp[ np.argmin( ((lat_temp-lat_mast[x,y])**2) + \
                                                  ((lon_temp-lon_mast[x,y])**2)   )  ]
            topidx_std[x,y] = topidx_std_temp[ np.argmin( ((lat_temp-lat_mast[x,y])**2) + \
                                                          ((lon_temp-lon_mast[x,y])**2)   )  ]
        
    badland_x, badland_y=np.where((mask_mast==1) & (topidx.data<0.))
    nbadpoints=len(badland_x)
    k=k+1
        
        

topidx=np.ma.masked_equal(topidx,outmv)
topidx_std=np.ma.masked_equal(topidx_std,outmv)


# Write to file
print  "opening file: "+outfile
outf=nc.Dataset(outfile,'w',format='NETCDF4_CLASSIC')

# Dimensions
outf.createDimension('x',dims[0])
outf.createDimension('y',dims[1])

# Variables
# EMEP i coordinate:
out_iemep=outf.createVariable('i_EMEP','i',('x',))
out_iemep.longname='EMEP i coordinate'
out_iemep.units='unitless'
out_iemep[:]=i_EMEP_mast[0,:]
# EMEP j coordinate:
out_jemep=outf.createVariable('j_EMEP','i',('y',))
out_jemep.longname='EMEP j coordinate'
out_jemep.units='unitless'
out_jemep[:]=j_EMEP_mast[:,0]
# Grid Centre Longitude:
out_lon=outf.createVariable('cen_lon','float64',('y','x'))
out_lon.longname='Centre Longitude of pixel'
out_lon.units='Degrees East'
out_lon[:]=lon_mast
# Grid Centre Latitude:
out_lat=outf.createVariable('cen_lat','float64',('y','x'))
out_lat.longname='Centre Latitude of pixel'
out_lat.units='Degrees North'
out_lat[:]=lat_mast
# topidx:
out_topidx=outf.createVariable('timean','float64',('y','x'),fill_value=outmv)
out_topidx.longname='Topogrphic Index'
out_topidx.units='Unitless'
out_topidx[:]=topidx
# topidx_STD:
out_topidxSTD=outf.createVariable('tistd','float64',('y','x'),fill_value=outmv)
out_topidxSTD.longname='Topogrphic Index Standard Deviation'
out_topidxSTD.units='Unitless'
out_topidxSTD[:]=topidx_std

outf.title="Topographic index data remapped to the EMEP4UK 5km grid"
outf.close()


# Plot data
PM.plot_map(topidx,lon_mast,lat_mast, \
            DATA_RANGE=[0.0,8.0], \
            MAP_TYPE='Mesh', MPL_CBAR='Oranges_r', NLEVELS=17, \
            WIDTH=12, HEIGHT=8, PLOT_LABEL='Topographic Index', CBAR_ORIENTATION='vertical', \
            PLOT_TITLE='Topographic Index on EMEP4UK 5km grid', FONTSIZES=[12,12,12,18], \
            iDISPLAY='N', FILE_PLOT=EMEP_dir+param_subdir+'topidx_emep4uk.png', \
            LATDEL=2, LONDEL=2, RESOLUTION='h', PROJECTION='stere')

# Plot data
PM.plot_map(topidx_std,lon_mast,lat_mast, \
            DATA_RANGE=[0.0,4.0], \
            MAP_TYPE='Mesh', MPL_CBAR='Oranges_r', NLEVELS=17, \
            WIDTH=12, HEIGHT=8, PLOT_LABEL='Topographic Index STD', CBAR_ORIENTATION='vertical', \
            PLOT_TITLE='Standard Deviation of Topographic Index on EMEP4UK 5km grid', FONTSIZES=[12,12,12,18], \
            iDISPLAY='N', FILE_PLOT=EMEP_dir+param_subdir+'topidx_std_emep4uk.png', \
            LATDEL=2, LONDEL=2, RESOLUTION='h', PROJECTION='stere')


#plt.subplot(1,3,1)
#plt.imshow(EMEP_topidx_mast,vmin=0)
#plt.colorbar(fraction=0.04)
#plt.subplot(1,3,2)
#plt.imshow(EMEP_topidx_std_mast,vmin=0)
#plt.colorbar(fraction=0.04)
#plt.subplot(1,3,3)
#plt.imshow(mask_mast)
#plt.colorbar(fraction=0.04)
#
#
#plt.show()
##



