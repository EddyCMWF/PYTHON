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
import plot_map_ECP as PM
#
#
#
EMEP_dir='/users/eow/edwcom/EMEP/'
#
param_subdir='topo2emep/'
param_filename='topidx_emep4ukEUROPE_grid_section.nc'
param_dataname='mu'
#
Master_Gridfile=EMEP_dir+'EMEP4UK/EMEP4UK_EUROPE_gridfile.nc'
#
section_list_file=EMEP_dir+'EMEP4UK/sections/EMEP4UK_EUROPE_gridfile_list.txt'
#
outfile=EMEP_dir+param_subdir+'topidx_emep4ukEUROPE_grid.nc'
outmv  = -9999.0
land_fillvalue = 3.0
land_fillvalue_std = 2.0
#
#
Mast_inf=nc.Dataset(Master_Gridfile,'r')
dims        = Mast_inf.variables['grid_dims'][:]
i_mast      = Mast_inf.variables['i'][:]
j_mast      = Mast_inf.variables['j'][:]
i_EMEP_mast = Mast_inf.variables['grid_EMEP_i_coord'][:].reshape(dims[::-1])
j_EMEP_mast = Mast_inf.variables['grid_EMEP_j_coord'][:].reshape(dims[::-1])
lat_mast    = Mast_inf.variables['grid_center_lat'][:].reshape(dims[::-1])
lon_mast    = Mast_inf.variables['grid_center_lon'][:].reshape(dims[::-1])
mask_mast   = Mast_inf.variables['grid_imask'][:].reshape(dims[::-1])
Mast_inf.close()
#
#
EMEP_topidx_mast=np.zeros_like(mask_mast)-9999.0
EMEP_topidx_std_mast=np.zeros_like(mask_mast)-9999.0
#
list=open(section_list_file).readlines()
#
for section in list:
    section_dir=EMEP_dir+param_subdir+section[:-1]+'/'
    Sect_inf=nc.Dataset(section_dir+param_filename,'r')
    topidx_sect     = Sect_inf.variables[param_dataname][:]
    topidx_std_sect = Sect_inf.variables[param_dataname+'_std'][:]
    Sect_inf.close()
    
    Section_Gridfile=EMEP_dir+param_subdir+section[:-1]+'/input_data/grid_files/EMEP4UK_EUROPE_gridfile.nc'
    SectGF_inf  = nc.Dataset(Section_Gridfile,'r')
    dims_sect   = SectGF_inf.variables['grid_dims'][:]
    i_EMEP_sect = SectGF_inf.variables['grid_EMEP_i_coord'][:].reshape(dims_sect[::-1])
    j_EMEP_sect = SectGF_inf.variables['grid_EMEP_j_coord'][:].reshape(dims_sect[::-1])
    #lat_sect    = SectGF_inf.variables['grid_center_lat'][:].reshape(dims_sect[::-1])
    #lon_sect    = SectGF_inf.variables['grid_center_lon'][:].reshape(dims_sect[::-1])
    SectGF_inf.close()
    
    x_SP, y_SP = np.where( (i_EMEP_mast==i_EMEP_sect[0,0]) & (j_EMEP_mast==j_EMEP_sect[0,0]) )
    x_EP, y_EP = np.where( (i_EMEP_mast==i_EMEP_sect[-1,-1]) & (j_EMEP_mast==j_EMEP_sect[-1,-1]) )
    
    EMEP_topidx_mast[x_SP:x_EP+1,y_SP:y_EP+1] = topidx_sect
    EMEP_topidx_std_mast[x_SP:x_EP+1,y_SP:y_EP+1] = topidx_std_sect


EMEP_topidx_mast[mask_mast==0]     = outmv
EMEP_topidx_std_mast[mask_mast==0] = outmv

########################################################################################
# Fill bad land points with closest good land point value
# index bad land points
badland_x, badland_y=np.where((mask_mast==1) & (EMEP_topidx_mast<0.))
nbadpoints=len(badland_x)

# create boolean array of where there are good land points
uber_mask = (mask_mast == 1) &  \
            (EMEP_topidx_mast > 0.)
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
            topidx_temp=EMEP_topidx_mast[max(x-radius,0):min(x+radius+1,dims[1]),  \
                                         max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            topidx_std_temp=EMEP_topidx_std_mast[max(x-radius,0):min(x+radius+1,dims[1]),  \
                                                 max(y-radius,0):min(y+radius+1,dims[0])][mask_temp]
            
            # locate closest good land point using lats and lons and fill gap with it
            EMEP_topidx_std_mast[x,y] = topidx_temp[ np.argmin( ((lat_temp-lat_mast[x,y])**2) + \
                                                                ((lon_temp-lon_mast[x,y])**2)   )  ]
            EMEP_topidx_mast[x,y]     = topidx_std_temp[ np.argmin( ((lat_temp-lat_mast[x,y])**2) + \
                                                                    ((lon_temp-lon_mast[x,y])**2)   )  ]
    
    badland_x, badland_y=np.where((mask_mast==1) & (EMEP_topidx_mast<0.))
    nbadpoints=len(badland_x)
    k=k+1


EMEP_topidx_mast=np.ma.masked_equal(EMEP_topidx_mast,outmv)
EMEP_topidx_std_mast=np.ma.masked_equal(EMEP_topidx_std_mast,outmv)


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
out_topidx[:]=EMEP_topidx_mast
# topidx_STD:
out_topidxSTD=outf.createVariable('tistd','float64',('y','x'),fill_value=outmv)
out_topidxSTD.longname='Topogrphic Index Standard Deviation'
out_topidxSTD.units='Unitless'
out_topidxSTD[:]=EMEP_topidx_std_mast

outf.title="Topographic index data remapped to the EMEP4UK-EUROPE grid"
outf.close()


# Plot data
PM.plot_map(EMEP_topidx_mast,lon_mast,lat_mast, \
            DATA_RANGE=[0.0,8.0], \
            MAP_TYPE='Mesh', MPL_CBAR='Oranges_r', NLEVELS=17, \
            WIDTH=12, HEIGHT=8, PLOT_LABEL='Topographic Index', CBAR_ORIENTATION='vertical', \
            PLOT_TITLE='Topographic Index on EMEP4UK-Europe grid', FONTSIZES=[12,12,12,18], \
            iDISPLAY='N', FILE_PLOT=EMEP_dir+param_subdir+'topidx_emep4ukEUROPE.png', \
            LATDEL=10, LONDEL=10, RESOLUTION='h', PROJECTION='stere')

# Plot data
PM.plot_map(EMEP_topidx_std_mast,lon_mast,lat_mast, \
            DATA_RANGE=[0.0,4.0], \
            MAP_TYPE='Mesh', MPL_CBAR='Oranges_r', NLEVELS=17, \
            WIDTH=12, HEIGHT=8, PLOT_LABEL='Topographic Index STD', CBAR_ORIENTATION='vertical', \
            PLOT_TITLE='Standard Deviation of Topographic Index on EMEP4UK-Europe grid', FONTSIZES=[12,12,12,18], \
            iDISPLAY='N', FILE_PLOT=EMEP_dir+param_subdir+'topidx_std_emep4ukEUROPE.png', \
            LATDEL=10, LONDEL=10, RESOLUTION='h', PROJECTION='stere')


