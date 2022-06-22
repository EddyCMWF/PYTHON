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
#import pylab as plt

PFT_Snames=['CF','DF','NF','BF','TC','MC','RC','SNL','GR','MS','WE','TU','DE','W','ICE','U']  #,'IAM_VEG','OLD_IAM_MF','OLD_IAM_DF']
PFT_indices=range(len(PFT_Snames))
PFT_note    = ''
for pft_sh in PFT_Snames:
    PFT_note+=pft_sh+' '

longnames=[ 'T/B Conif',     \
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
            'Urban'          ]
    #'Generic crop',  \
    #'Generic DF',    \
    #'Generic MF'     ]

PFT_Snames_OUT=['CF','DF','NF','BF','TC','MC','RC','SNL','GR','MS','SOIL','W','ICE','U']
PFT_indices_OUT=range(len(PFT_Snames_OUT))
PFT_note_OUT  = ''
for pft_sh in PFT_Snames_OUT:
    PFT_note_OUT+=pft_sh+' '

W_INDEX=11

EMEP4UK_DIR='/users/eow/edwcom/EMEP/EMEP4UK/'


EMEP4UK_landusefile=EMEP4UK_DIR+'EMEP4UK_Landuse.nc'
WRFfile_UK='/users/eow/edwcom/EMEP/EMEP4UK/WRF_DATA/wrfout_d03_2001-12-30.nc'
EMEP4UK_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_LandFrac.nc.noIAM.oneSOIL'


EMEP4UK_EURO_landusefile=EMEP4UK_DIR+'EMEP4UK_EUROPE_Landuse.nc'
WRFfile_EURO='/users/eow/edwcom/EMEP/EMEP4UK/WRF_DATA/wrfout_d01_2001-12-30.nc'
EMEP4UK_EURO_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_LandFrac.nc.noIAM.oneSOIL'



########################################################################################################

#fetch LS mask from a sample drive file
WRFin = nc.Dataset(WRFfile_UK,'r')
WRF_LSmask  = WRFin.variables['XLAND'][0,:,:]
WRFin.close()

#fetch landcover fractions from EMEP4UK landuse files
inf_EMUK_LU = nc.Dataset(EMEP4UK_landusefile,'r')
    
       
#fetch data from EMEP4UK landuse file
inf_EMUK = nc.Dataset(EMEP4UK_landusefile,'r')
EMEP4UK_lats=inf_EMUK.variables['lat'][:]
EMEP4UK_lons=inf_EMUK.variables['lon'][:]
EMEP4UK_i=inf_EMUK.variables['i_EMEP'][:]
EMEP4UK_j=inf_EMUK.variables['j_EMEP'][:]
UK_raw_totals=inf_EMUK.variables['Total'][:]

# extract each fraction for each surface type and append to list
EMEP4UK_LF_data=[]
for PFT_num in range(len(PFT_Snames)):
    PFT=PFT_Snames[PFT_num]
    EMEP4UK_LF_data.append(inf_EMUK.variables[PFT][:])
    #plotdata=np.ma.masked_less_equal(inf_EMUK.variables[PFT][:],0)    
    #PT.plot_map(plotdata,EMEP4UK_lons,EMEP4UK_lats, \
    #            DATA_RANGE=[0,1], NLEVELS=11, \
    #            MPL_CBAR='YlOrRd',CBAR_LABEL='Fractional Cover',PLOT_TITLE=longnames[PFT_num],
    #            FILE_PLOT=EMEP4UK_DIR+'plots/EMEP4UK_'+PFT+'.png', \
    #            LATDEL=2., LONDEL=2., RESOLUTION='h',PROJECTION='stere', \
    #            LON_RANGE=[-14,11],LAT_RANGE=[51,57] )

inf_EMUK.close()

# Compile the soil classes into one class
# first elements 0-9 are the PFTs and can be copied straight over:
EMEP4UK_LF_array_oneSOIL=EMEP4UK_LF_data[:10]
# new 10th element is the combined soil, i.e. WE + TU + DE = SOIL
EMEP4UK_LF_array_oneSOIL.append( np.sum(np.array(EMEP4UK_LF_data[10:13]),axis=0) )
#remaining elements are the same so copied into remaining slots
for i in range(13,16):
    EMEP4UK_LF_array_oneSOIL.append(EMEP4UK_LF_data[i])

#Create Land Frac, array of ones minus Water
EMEP4UK_LSmask=(np.ones_like(EMEP4UK_LF_array_oneSOIL[W_INDEX])-EMEP4UK_LF_array_oneSOIL[W_INDEX])
#Apply WRF LSmask
EMEP4UK_LSmask[WRF_LSmask==2]=0.

#Remask
EMEP4UK_LSmask[EMEP4UK_LF_array_oneSOIL[W_INDEX]<-1.]=-9999.

EMEP4UK_LSmask=np.ma.masked_equal(EMEP4UK_LSmask,-9999)

# Convert data from list to array
EMEP4UK_LF_array_oneSOIL=np.array(EMEP4UK_LF_array_oneSOIL)
# change fill_filue to -9999.0
EMEP4UK_LF_array_oneSOIL[EMEP4UK_LF_array_oneSOIL<-1.]=-9999.0
EMEP4UK_LF_array_oneSOIL=np.ma.masked_equal(EMEP4UK_LF_array_oneSOIL,-9999.)



########################################################################################################

#fetch LS mask from a sample drive file
WRFin = nc.Dataset(WRFfile_EURO,'r')
WRF_LSmask  = WRFin.variables['XLAND'][0,:,:]
WRFin.close()

#fetch data from EMEP4UK landuse file
inf_EMUK_EU = nc.Dataset(EMEP4UK_EURO_landusefile,'r')
EMEP4UK_EURO_lats=inf_EMUK_EU.variables['lat'][:]
EMEP4UK_EURO_lons=inf_EMUK_EU.variables['lon'][:]
EMEP4UK_EURO_i=inf_EMUK_EU.variables['i_EMEP'][:]
EMEP4UK_EURO_j=inf_EMUK_EU.variables['j_EMEP'][:]
EURO_raw_totals=inf_EMUK_EU.variables['Total'][:]
EMEP4UK_EURO_LF_data=[]
for PFT_num in range(len(PFT_Snames)):
    PFT=PFT_Snames[PFT_num]
    EMEP4UK_EURO_LF_data.append(inf_EMUK_EU.variables[PFT][:])
    #plotdata=np.ma.masked_less_equal(inf_EMUK_EU.variables[PFT][:],0)
    #PT.plot_map(plotdata,EMEP4UK_EURO_lons,EMEP4UK_EURO_lats, \
    #            DATA_RANGE=[0,1], NLEVELS=11, \
    #            MPL_CBAR='YlOrRd',CBAR_LABEL='Fractional Cover',PLOT_TITLE=longnames[PFT_num],
    #            FILE_PLOT=EMEP4UK_DIR+'plots/EMEP4UK_EURO_'+PFT+'.png', \
    #            LATDEL=10., LONDEL=10., RESOLUTION='h',PROJECTION='stere' )

inf_EMUK_EU.close()

# Compile the soil classes into one class

# first elements 0-9 are the PFTs and can be copied straight over:
EMEP4UK_EURO_LF_array_oneSOIL=EMEP4UK_EURO_LF_data[:10]
# new 10th element is the combined soil, i.e. WE + TU + DE = SOIL
EMEP4UK_EURO_LF_array_oneSOIL.append( np.sum(np.array(EMEP4UK_EURO_LF_data[10:13]),axis=0) )
#remaining elements are the same so copied into remaining slots
for i in range(13,16):
    EMEP4UK_EURO_LF_array_oneSOIL.append(EMEP4UK_EURO_LF_data[i])

#Create Land Frac, array of ones minus Water
EMEP4UK_EURO_LSmask=(np.ones_like(EMEP4UK_EURO_LF_array_oneSOIL[W_INDEX].data)-EMEP4UK_EURO_LF_array_oneSOIL[W_INDEX].data)
#Apply WRF LSmask
EMEP4UK_EURO_LSmask[WRF_LSmask==2]=0.

# remask with fill_value of -9999.0
EMEP4UK_EURO_LSmask[EMEP4UK_EURO_LF_array_oneSOIL[W_INDEX].mask]=-9999.
EMEP4UK_EURO_LSmask=np.ma.masked_equal(EMEP4UK_EURO_LSmask,-9999)

# Convert data from list to array
EMEP4UK_EURO_LF_array_oneSOIL=np.array(EMEP4UK_EURO_LF_array_oneSOIL)

# change fill_filue to -9999.0
EMEP4UK_EURO_LF_array_oneSOIL[EMEP4UK_EURO_LF_array_oneSOIL<-1.]=-9999.0
EMEP4UK_EURO_LF_array_oneSOIL=np.ma.masked_equal(EMEP4UK_EURO_LF_array_oneSOIL,-9999.)



########################################################################################################

#Write out EMEP4UK data to file
outf=nc.Dataset(EMEP4UK_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('west_east',len(EMEP4UK_i))
outf.createDimension('south_north',len(EMEP4UK_j))
outf.createDimension('pseudo',len(PFT_Snames_OUT))

#create and store dimension variables
out_i_EMEP4UK=outf.createVariable('i_EMEP','float32',('west_east'))
out_i_EMEP4UK.long_name='EMEP grid i coordinate'
out_i_EMEP4UK.units='unitless'
out_i_EMEP4UK[:]=EMEP4UK_i

out_j_EMEP4UK=outf.createVariable('j_EMEP','float32',('south_north',))
out_j_EMEP4UK.long_name='EMEP grid j coordinate'
out_j_EMEP4UK.units='unitless'
out_j_EMEP4UK[:]=EMEP4UK_j

out_pseudo_EMEP4UK=outf.createVariable('pseudo','i',('pseudo',))
out_pseudo_EMEP4UK.long_name='Land Surface Type'
out_pseudo_EMEP4UK.units='unitless'
out_pseudo_EMEP4UK.note=PFT_note_OUT
out_pseudo_EMEP4UK.note2='Wetlands, Tundra and Desert combined into soil classification'
out_pseudo_EMEP4UK[:]=PFT_indices_OUT

# Create and store PFT_fractions and LS mask variables
out_lf_EMEP4UK=outf.createVariable('Land_Frac','float32',('pseudo','south_north','west_east'),fill_value=-9999.)
out_lf_EMEP4UK.long_name='Landuse Fraction for each PFT'
out_lf_EMEP4UK.units='Fraction of Gridcell'
out_lf_EMEP4UK[:]=EMEP4UK_LF_array_oneSOIL

out_LS_EMEP4UK=outf.createVariable('LSmask','float32',('south_north','west_east'),fill_value=-9999.)
out_LS_EMEP4UK.long_name='Land-Sea mask'
out_LS_EMEP4UK.units='Fraction of Gridcell'
out_LS_EMEP4UK[:]=EMEP4UK_LSmask

#Global Attributes
outf.title='Land Fraction for each PFT on the EMEP4UK 5km grid'
outf.note='Reference: Simpson et al., 2012, The EMEP MSC-W chemical transport model: technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()




########################################################################################################

#Write out EMEP4UK data to file
outf=nc.Dataset(EMEP4UK_EURO_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('west_east',len(EMEP4UK_EURO_i))
outf.createDimension('south_north',len(EMEP4UK_EURO_j))
outf.createDimension('pseudo',len(PFT_Snames_OUT))

#create and store dimension variables
out_i_EMEP4UK_EU=outf.createVariable('i_EMEP','float32',('west_east',))
out_i_EMEP4UK_EU.long_name='EMEP grid i coordinate'
out_i_EMEP4UK_EU.units='unitless'
out_i_EMEP4UK_EU[:]=EMEP4UK_EURO_i

out_j_EMEP4UK_EU=outf.createVariable('j_EMEP','float32',('south_north',))
out_j_EMEP4UK_EU.long_name='EMEP grid j coordinate'
out_j_EMEP4UK_EU.units='unitless'
out_j_EMEP4UK_EU[:]=EMEP4UK_EURO_j

out_pseudo_EMEP4UK=outf.createVariable('pseudo','i',('pseudo',))
out_pseudo_EMEP4UK.long_name='Land Surface Type'
out_pseudo_EMEP4UK.units='unitless'
out_pseudo_EMEP4UK.note=PFT_note_OUT
out_pseudo_EMEP4UK.note2='Wetlands, Tundra and Desert combined into soil classification'
out_pseudo_EMEP4UK[:]=PFT_indices_OUT

# Create and store lai and canht variables
out_lai_EMEP4UK_EU=outf.createVariable('Land_Frac','float32',('pseudo','south_north','west_east'),fill_value=-9999.)
out_lai_EMEP4UK_EU.long_name='Landuse Fraction for each PFT'
out_lai_EMEP4UK_EU.units='Fraction of Gridcell'
out_lai_EMEP4UK_EU[:]=EMEP4UK_EURO_LF_array_oneSOIL

out_LS_EMEP4UK=outf.createVariable('LSmask','float32',('south_north','west_east'),fill_value=-9999.)
out_LS_EMEP4UK.long_name='Land-Sea mask'
out_LS_EMEP4UK.units='Fraction of Gridcell'
out_LS_EMEP4UK[:]=EMEP4UK_EURO_LSmask

#Global Attributes
outf.title='Land Fraction for each PFT on the EMEP4UK-EUROPE 50km grid'
outf.note='Reference: Simpson et al., 2012, The EMEP MSC-W chemical transport model: technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()


########################################################################################################
