#!/usr/bin/env python
##################################################################################
#
# Program: EMEP4UK_create_LandFrac_file.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To output an Landfrac for EMEP4UK for JULES run
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
#import plot_tools as PT

PFT_Snames=['CF','DF','NF','BF','TC','MC','RC','SNL','GR','MS','WE','TU','DE','W','ICE','U']  #,'IAM_VEG','OLD_IAM_MF','OLD_IAM_DF']
PFT_indices=range(len(PFT_Snames))
PFT_note    = ''
for pft_sh in PFT_Snames:
    PFT_note+=pft_sh+' '


EMEP4UK_DIR='/users/eow/edwcom/EMEP/EMEP4UK/'


EMEP4UK_landusefile=EMEP4UK_DIR+'EMEP4UK_Landuse.nc'
EMEP4UK_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_LandFrac.nc.noIAM'


EMEP4UK_EURO_landusefile=EMEP4UK_DIR+'EMEP4UK_EUROPE_Landuse.nc'
EMEP4UK_EURO_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_LandFrac.nc.noIAM'



#fetch landcover fractions from EMEP4UK landuse files
inf_EMUK_LU = nc.Dataset(EMEP4UK_landusefile,'r')



#fetch data from EMEP4UK landuse file
inf_EMUK = nc.Dataset(EMEP4UK_landusefile,'r')
EMEP4UK_lats=inf_EMUK.variables['lat'][:]
EMEP4UK_lons=inf_EMUK.variables['lon'][:]
EMEP4UK_i=inf_EMUK.variables['i_EMEP'][:]
EMEP4UK_j=inf_EMUK.variables['j_EMEP'][:]
EMEP4UK_LF_data=[]
for PFT in PFT_Snames:
    EMEP4UK_LF_data.append(inf_EMUK.variables[PFT][:])

inf_EMUK.close()
EMEP4UK_LF_data=np.array(EMEP4UK_LF_data)


#fetch data from EMEP4UK landuse file
inf_EMUK_EU = nc.Dataset(EMEP4UK_EURO_landusefile,'r')
EMEP4UK_EURO_lats=inf_EMUK_EU.variables['lat'][:]
EMEP4UK_EURO_lons=inf_EMUK_EU.variables['lon'][:]
EMEP4UK_EURO_i=inf_EMUK_EU.variables['i_EMEP'][:]
EMEP4UK_EURO_j=inf_EMUK_EU.variables['j_EMEP'][:]
EMEP4UK_EURO_LF_data=[]
for PFT in PFT_Snames:
    EMEP4UK_EURO_LF_data.append(inf_EMUK_EU.variables[PFT][:])

inf_EMUK_EU.close()
EMEP4UK_EURO_LF_data=np.array(EMEP4UK_EURO_LF_data)


#PT.plot_map(EMEP4UK_EURO_LF_data[1,:,:], \
#            EMEP4UK_EURO_lons, \
#            EMEP4UK_EURO_lats, \
#            DATA_RANGE=[0.0,1.0], \
#            MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=11,CBAR_ORIENTATION='vertical', \
#            WIDTH=12, HEIGHT=8, CBAR_LABEL='Water Fraction',\
#            PLOT_TITLE='Water Fraction', FONTSIZES=[12,12,12,18], \
#            iDISPLAY='Y', \
#            LATDEL=10., LONDEL=10., RESOLUTION='h', 
#            PROJECTION='stere')
#
#PT.plot_map(EMEP4UK_LF_data[13,:,:], \
#            EMEP4UK_lons, \
#            EMEP4UK_lats, \
#            DATA_RANGE=[0.0,1.0], \
#            MAP_TYPE='Mesh', MPL_CBAR='RdYlBu_r', NLEVELS=11,CBAR_ORIENTATION='vertical', \
#            WIDTH=12, HEIGHT=8, CBAR_LABEL='Water Fraction',\
#            PLOT_TITLE='Water Fraction', FONTSIZES=[12,12,12,18], \
#            iDISPLAY='Y', \
#            LATDEL=10., LONDEL=10., RESOLUTION='h', 
#            PROJECTION='stere')



#Write out EMEP4UK data to file
outf=nc.Dataset(EMEP4UK_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('x',len(EMEP4UK_i))
outf.createDimension('y',len(EMEP4UK_j))
outf.createDimension('pseudo',len(PFT_Snames))

#create and store dimension variables
out_i_EMEP4UK=outf.createVariable('x','float32',('x',))
out_i_EMEP4UK.long_name='EMEP grid i coordinate'
out_i_EMEP4UK.units='unitless'
out_i_EMEP4UK[:]=EMEP4UK_i

out_j_EMEP4UK=outf.createVariable('y','float32',('y',))
out_j_EMEP4UK.long_name='EMEP grid j coordinate'
out_j_EMEP4UK.units='unitless'
out_j_EMEP4UK[:]=EMEP4UK_j

out_pseudo_EMEP4UK=outf.createVariable('pseudo','i',('pseudo',))
out_pseudo_EMEP4UK.long_name='Land Surface Type'
out_pseudo_EMEP4UK.units='unitless'
out_pseudo_EMEP4UK.note=PFT_note
out_pseudo_EMEP4UK[:]=PFT_indices

# Create and store PFT_fractions and LS mask variables
out_lai_EMEP4UK=outf.createVariable('Land_Frac','float32',('pseudo','y','x'))
out_lai_EMEP4UK.long_name='Landuse Fraction for each PFT'
out_lai_EMEP4UK.units='Fraction of Gridcell'
out_lai_EMEP4UK[:]=EMEP4UK_LF_data

out_LS_EMEP4UK=outf.createVariable('LSmask','float32',('y','x'))
out_LS_EMEP4UK.long_name='Land-Sea mask'
out_LS_EMEP4UK.units='Fraction of Gridcell'
out_LS_EMEP4UK[:]=(np.ones_like(EMEP4UK_LF_data[PFT_Snames=='W',:,:])-EMEP4UK_LF_data[PFT_Snames=='W',:,:])

#Global Attributes
outf.title='Land Fraction for each PFT on the EMEP4UK 5km grid'
outf.note='Reference: Simpson et al., 2012, The EMEP MSC-W chemical transport model: technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()



#Write out EMEP4UK data to file
outf=nc.Dataset(EMEP4UK_EURO_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('x',len(EMEP4UK_EURO_i))
outf.createDimension('y',len(EMEP4UK_EURO_j))
outf.createDimension('pseudo',len(PFT_Snames))

#create and store dimension variables
out_i_EMEP4UK_EU=outf.createVariable('x','float32',('x',))
out_i_EMEP4UK_EU.long_name='EMEP grid i coordinate'
out_i_EMEP4UK_EU.units='unitless'
out_i_EMEP4UK_EU[:]=EMEP4UK_EURO_i

out_j_EMEP4UK_EU=outf.createVariable('y','float32',('y',))
out_j_EMEP4UK_EU.long_name='EMEP grid j coordinate'
out_j_EMEP4UK_EU.units='unitless'
out_j_EMEP4UK_EU[:]=EMEP4UK_EURO_j

out_pseudo_EMEP4UK_EU=outf.createVariable('pseudo','i',('pseudo',))
out_pseudo_EMEP4UK_EU.long_name='Land Surface Type'
out_pseudo_EMEP4UK_EU.units='unitless'
out_pseudo_EMEP4UK_EU.note=PFT_note
out_pseudo_EMEP4UK_EU[:]=PFT_indices

# Create and store lai and canht variables
out_lai_EMEP4UK_EU=outf.createVariable('Land_Frac','float32',('pseudo','y','x'))
out_lai_EMEP4UK_EU.long_name='Landuse Fraction for each PFT'
out_lai_EMEP4UK_EU.units='Fraction of Gridcell'
out_lai_EMEP4UK_EU[:]=EMEP4UK_EURO_LF_data

out_LS_EMEP4UK=outf.createVariable('LSmask','float32',('y','x'))
out_LS_EMEP4UK.long_name='Land-Sea mask'
out_LS_EMEP4UK.units='Fraction of Gridcell'
out_LS_EMEP4UK[:]=(np.ones_like(EMEP4UK_EURO_LF_data[PFT_Snames=='W',:,:])-EMEP4UK_EURO_LF_data[PFT_Snames=='W',:,:])

#Global Attributes
outf.title='Land Fraction for each PFT on the EMEP4UK-EUROPE 50km grid'
outf.note='Reference: Simpson et al., 2012, The EMEP MSC-W chemical transport model: technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()

