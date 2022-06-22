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
import argparse
#import pylab as plt

EMEP4UK_gridfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_gridfile.nc'

EMEP4UK_LCfile='/users/eow/edwcom/EMEP/EMEP4UK/Landuse_PS_5km_LC.nc'

outfile = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Landuse.nc'

#Extract EMEP4UK grid data from file
infGRID=nc.Dataset(EMEP4UK_gridfile,'r')
grid_dims=infGRID.variables['grid_dims'][::-1]
grid_i_EMEP=np.array(np.round(infGRID.variables['grid_EMEP_i_coord'][:].reshape(grid_dims)[0,:]*10.),dtype='int')/10.
grid_j_EMEP=np.array(np.round(infGRID.variables['grid_EMEP_j_coord'][:].reshape(grid_dims)[:,0]*10.),dtype='int')/10.
grid_lat=infGRID.variables['grid_center_lat'][:].reshape(grid_dims)
grid_lon=infGRID.variables['grid_center_lon'][:].reshape(grid_dims)
grid_imask=infGRID.variables['grid_imask'][:].reshape(grid_dims)
infGRID.close()

# open Landuse file
infLC=nc.Dataset(EMEP4UK_LCfile,'r')
LC_lat=infLC.variables['lat'][:]
LC_lon=infLC.variables['lon'][:]
LC_i=infLC.variables['i'][:]
LC_j=infLC.variables['j'][:]


# Because someone has decided to use yet another slightly different EMEP type grid 
# the best option is to find the minimum of the lat/lon corner coordinates
temp=((LC_lon-grid_lon[0,0])**2)+((LC_lat-grid_lat[0,0])**2)
LL_cor=[np.argmin(temp)/len(infLC.dimensions['i'])+1,\
        np.argmin(temp)%len(infLC.dimensions['i'])-1 ]

temp=((LC_lon-grid_lon[-1,-1])**2)+((LC_lat-grid_lat[-1,-1])**2)
UR_cor=[np.argmin(temp)/len(infLC.dimensions['i'])+1,\
        np.argmin(temp)%len(infLC.dimensions['i'])-1 ]
del temp

LC_lat_sect=LC_lat[LL_cor[0]:UR_cor[0]+1,LL_cor[1]:UR_cor[1]+1]
LC_lon_sect=LC_lon[LL_cor[0]:UR_cor[0]+1,LL_cor[1]:UR_cor[1]+1]

LC_i_sect=((LC_i[LL_cor[1]:UR_cor[1]+1]+2.5)/50.)+7.9
LC_j_sect=((LC_j[LL_cor[0]:UR_cor[0]+1]+2.5)/50.)+110.



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
     

#Calculate Totals array as data need to be normalised:
Totals = np.zeros_like(LC_lon_sect)  
for var in range(len(infnames)):
    if not (outf_longnames[var][:7]=='Generic'):
        Totals+=infLC.variables[infnames[var]][LL_cor[0]:UR_cor[0]+1,LL_cor[1]:UR_cor[1]+1]

# Open out file to write data
outf=nc.Dataset(outfile,'w',format='NETCDF4_CLASSIC')
outf.createDimension('i',len(grid_i_EMEP))
outf.createDimension('j',len(grid_j_EMEP))

out_i_EMEP=outf.createVariable('i_EMEP','float32',('i',))
out_i_EMEP.long_name='EMEP grid i coordinate'
out_i_EMEP.units='unitless'
out_i_EMEP[:]=grid_i_EMEP

out_j_EMEP=outf.createVariable('j_EMEP','float32',('j',))
out_j_EMEP.long_name='EMEP grid j coordinate'
out_j_EMEP.units='unitless'
out_j_EMEP[:]=grid_j_EMEP

out_lat=outf.createVariable('lat','float64',('j','i'))
out_lat.long_name='latitude'
out_lat.units='degrees_north'
out_lat[:]=LC_lat_sect

out_lon=outf.createVariable('lon','float64',('j','i'))
out_lon.long_name='longitude'
out_lon.units='degrees_east'
out_lon[:]=LC_lon_sect
 
# Save each land cover type
for var in range(len(infnames)):
    print infnames[var], outfnames[var], outf_longnames[var]
    out_var=outf.createVariable(outfnames[var],'float64',('j','i'),fill_value=infLC.variables[infnames[var]]._FillValue)
    out_var.long_name=outf_longnames[var]
    out_var.units='Fractional Cover'
    # extract section of data from infile and divide by total to normalise
    outdata=infLC.variables[infnames[var]][LL_cor[0]:UR_cor[0]+1,LL_cor[1]:UR_cor[1]+1]/Totals
    out_var[:]=outdata
    del outdata


out_var=outf.createVariable('Total','float64',('j','i'),fill_value=-9999.0)
out_var.long_name="Total of non IAM components (should equal 1)"
out_var.units='Fractional Cover'
out_var[:]=Totals/Totals

outf.Conventions='CF-1.0'
outf.projection='Stereographic'

outf.close()



infLC.close()
