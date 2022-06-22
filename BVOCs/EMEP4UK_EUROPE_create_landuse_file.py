#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To remap data based on SCRIP output
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
#import pylab as plt
#import time

EMEP4UK_EUROPE_gridfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_gridfile.nc'

EMEP4UK_LCfile='/users/eow/edwcom/EMEP/EMEP4UK/Landuse_PS_5km_LC.nc'

outfile = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_Landuse.nc'


#Extract EMEP4UK-EUROPE grid data from file
infGRID=nc.Dataset(EMEP4UK_EUROPE_gridfile,'r')
grid_dims_EU=infGRID.variables['grid_dims'][::-1]
grid_i_EMEP_EU=np.array(np.round(infGRID.variables['grid_EMEP_i_coord'][:].reshape(grid_dims_EU)[0,:]*10.),dtype='int')/10.
grid_j_EMEP_EU=np.array(np.round(infGRID.variables['grid_EMEP_j_coord'][:].reshape(grid_dims_EU)[:,0]*10.),dtype='int')/10.
grid_lat_EU=infGRID.variables['grid_center_lat'][:].reshape(grid_dims_EU)
grid_lon_EU=infGRID.variables['grid_center_lon'][:].reshape(grid_dims_EU)

grid_imask=infGRID.variables['grid_imask'][:].reshape(grid_dims_EU)
j_temp=infGRID.variables['grid_EMEP_j_coord'][:].reshape(grid_dims_EU)
i_temp=infGRID.variables['grid_EMEP_i_coord'][:].reshape(grid_dims_EU)

infGRID.close()


# open Landuse file
infLC=nc.Dataset(EMEP4UK_LCfile,'r')
LC_lat=infLC.variables['lat'][:]
LC_lon=infLC.variables['lon'][:]
LC_i=((infLC.variables['i'][:]+2.5)/50.)+7.9
LC_j=((infLC.variables['j'][:]+2.5)/50.)+110.

LC_i_grid, LC_j_grid = np.meshgrid(LC_i,LC_j)

#LC_i=((infLC.variables['i'][:]+2.5)/50.)+7.9
#LC_j=((infLC.variables['j'][:]+2.5)/50.)+8.

# Create index arrays for i and j (add 0.0001 to round 0.5 upwards)
index_i=np.round(LC_i_grid-np.amin(grid_i_EMEP_EU)+0.00001).flatten()
index_j=np.round(LC_j_grid-np.amin(grid_j_EMEP_EU)+0.00001).flatten()

index_ji     = ((index_j*grid_dims_EU[1])+index_i)
index_ji     = index_ji[index_ji<len(grid_lat_EU.flatten())]
uni_index_ji = set( index_ji )
uber_index=[ np.where(index_ji==idx)[0]  for idx in uni_index_ji ]


# Open out file to write data
outf=nc.Dataset(outfile,'w',format='NETCDF4_CLASSIC')
outf.createDimension('i',len(grid_i_EMEP_EU))
outf.createDimension('j',len(grid_j_EMEP_EU))

out_i_EMEP=outf.createVariable('i_EMEP','float32',('i',))
out_i_EMEP.long_name='EMEP grid i coordinate'
out_i_EMEP.units='unitless'
out_i_EMEP[:]=grid_i_EMEP_EU

out_j_EMEP=outf.createVariable('j_EMEP','float32',('j',))
out_j_EMEP.long_name='EMEP grid j coordinate'
out_j_EMEP.units='unitless'
out_j_EMEP[:]=grid_j_EMEP_EU

out_lat=outf.createVariable('lat','float64',('j','i'))
out_lat.long_name='latitude'
out_lat.units='degrees_north'
out_lat[:]=grid_lat_EU

out_lon=outf.createVariable('lon','float64',('j','i'))
out_lon.long_name='longitude'
out_lon.units='degrees_east'
out_lon[:]=grid_lon_EU
            
infnames=['LC:CF:EMEP','LC:DF:EMEP','LC:NF:EMEP','LC:BF:EMEP','LC:TC:EMEP', \
          'LC:MC:EMEP','LC:RC:EMEP','LC:SNL:EMEP','LC:GR:EMEP','LC:MS:EMEP', \
          'LC:WE:EMEP','LC:TU:EMEP','LC:DE:EMEP','LC:W:EMEP','LC:ICE:EMEP', \
          'LC:U:EMEP','LC:IAM_VEG:EMEP','OLD:IAM_MF:EMEP','OLD:IAM_DF:EMEP' ]

outfnames=['CF','DF','NF','BF','TC','MC','RC','SNL','GR','MS','WE','TU','DE','W','ICE','U',\
           'IAM_VEG','OLD_IAM_MF','OLD_IAM_DF']

outf_longnames=[ 'T/B Conif',     \
                 'T/B Decid',     \
                 'Med. Needle',   \
                 'Med Broadleaf', \
                 'T/B crop',      \
                 'Med. crop',     \
                 'Root crop',     \
                 'Moorland',      \
                 'Grass',         \
                 'Med. scrub',    \
                 'Wetlands',      \
                 'Tundra',        \
                 'Desert',        \
                 'Water',         \
                 'Ice',           \
                 'Urban',         \
                 'Generic crop',  \
                 'Generic DF',    \
                 'Generic MF'     ]


#indata=[]
outdata=[]
#Calculate Totals whilst reading data
TOTAL=np.zeros_like(grid_lat_EU)
for var in range(len(infnames)):
    print infnames[var], outfnames[var], outf_longnames[var]
    # Read in data as an array
    indata_temp=(infLC.variables[infnames[var]][:])
    # For each 50km grid cell (in uber_index) mean the appropriate points and append to list
    outdata_temp=[]
    for l in uber_index:
        outdata_temp.append(np.mean(indata_temp.flat[l]))
    del indata_temp
    # Create np array for outdata
    outdata_temp2=np.zeros_like(grid_lat_EU).flatten()+(infLC.variables[infnames[var]]._FillValue)
    # insert mean points in approriate locations
    outdata_temp2[list(uni_index_ji)]=np.array(outdata_temp)
    del outdata_temp
    # reshape to correct dimensions and append outdata to list
    outdata.append(outdata_temp2.reshape(grid_dims_EU))
    #
    # Do not include the IAM components in Totalling (final 3)
    if not (outf_longnames[var][:7]=='Generic'):
        TOTAL+=outdata_temp2.reshape(grid_dims_EU)
    del outdata_temp2

# set bad points to fill value
TOTAL[TOTAL<0]=infLC.variables[infnames[var]]._FillValue
TOTAL=np.ma.masked_equal(TOTAL,infLC.variables[infnames[var]]._FillValue)


out_var=outf.createVariable('Total','float64',('j','i'),fill_value=infLC.variables['Total']._FillValue)
out_var.long_name='Total of non IAM components (should equal 1)'
out_var.units='Fractional Cover'
out_var[:]=TOTAL

for var in range(len(infnames)):
    out_var=outf.createVariable(outfnames[var],'float64',('j','i'),fill_value=infLC.variables[infnames[var]]._FillValue)
    out_var.long_name=outf_longnames[var]
    out_var.units='Fractional Cover'
    out_var[:]=np.ma.masked_equal(outdata[var],infLC.variables[infnames[var]]._FillValue)


outf.Conventions='CF-1.0'
outf.projection='Stereographic'

outf.close()

infLC.close()


#n, bins, patches = plt.hist(TOTAL[TOTAL>0.9].flat, 50, normed=1, histtype='stepfilled')
#plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
#plt.show()
##
#
#
#
#
#
#plt.subplot(2,4,1)
#plt.imshow(i_temp,origin='bottom',vmin=0,vmax=120)
#plt.colorbar()
#plt.subplot(2,4,2)
#plt.imshow(j_temp,origin='bottom',vmin=0,vmax=120)
#plt.colorbar()
#plt.subplot(2,4,3)
#plt.imshow(grid_lat_EU,origin='bottom',vmin=20,vmax=90)
#plt.colorbar()
#plt.subplot(2,4,4)
#plt.imshow(grid_lon_EU,origin='bottom')
#plt.colorbar()
##
#
#plt.subplot(2,4,5)
#plt.imshow(LC_i_grid,origin='bottom',vmin=0,vmax=120)
#plt.colorbar()
#plt.subplot(2,4,6)
#plt.imshow(LC_j_grid,origin='bottom',vmin=0,vmax=120)
#plt.colorbar()
#plt.subplot(2,4,7)
#plt.imshow(LC_lat,origin='bottom',vmin=20,vmax=90)
#plt.colorbar()
#plt.subplot(2,4,8)
#plt.imshow(LC_lon,origin='bottom')
#plt.colorbar()
#
#
#plt.show()
#
#
#plt.subplot(1,2,1)
#plt.imshow(grid_imask,origin='bottom',vmin=0,vmax=1,interpolation='none')
#plt.colorbar()
#plt.subplot(1,2,2)
#plt.imshow(indata_orig,origin='bottom',vmin=0,vmax=1,interpolation='none')
#plt.colorbar()
#plt.show()
#
