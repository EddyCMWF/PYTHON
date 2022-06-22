#!/usr/bin/env python
##################################################################################
#
# Program: EMEP4UK_create_LAI_file.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To output an LAI seasonality for EMEP4UK based on the EMEP method
#           see Simpson et al 2012, The EMEP MSC-W chemical transport model
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
import EMEP_PFT_tools as EPT

#Create list of PFT Short_names
PFT_sh_names=['CF','DF','NF','BF','TC','MC','RC','SNL','GR','MS']   #,'IAM_CR','IAM_DF','IAM_MF']
PFT_indices =[ 1  , 2  , 3  , 4  , 5  , 6  , 7  , 8   , 9  , 10 ]   #,  17    , 18     , 19     ]

PFT_note    = ''
for pft_sh in PFT_sh_names:
    PFT_note=PFT_note+pft_sh+' '

LAIF=EPT.LAI_functions()
CHF =EPT.CH_functions()

EMEP4UK_DIR='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_base_output/'

EMEP4UK_basefile=EMEP4UK_DIR+'EMEP4UK_Base_emep_4.3_2001_day.nc'
EMEP4UK_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_LAI.nc.noIAM'

EMEP_basefile=EMEP4UK_DIR+'EMEP4UK_Base_EUROPE_emep_4.3_2001_day.nc'
EMEP_outfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_LAI.nc.noIAM'

#fetch lat and lon data from EMEP4UK base files
inf_EMUK = nc.Dataset(EMEP4UK_basefile,'r')
EMEP4UK_lats=inf_EMUK.variables['lat'][:]
EMEP4UK_lons=inf_EMUK.variables['lon'][:]
EMEP4UK_i=inf_EMUK.variables['i'][:]
EMEP4UK_j=inf_EMUK.variables['j'][:]
EMEP4UK_time=inf_EMUK.variables['time'][:]
inf_EMUK.close()


#fetch lat and lon data from EMEP4UK-EUROPE base files
inf_EMUK_EU = nc.Dataset(EMEP_basefile,'r')
EMEP4UK_EURO_lats=inf_EMUK_EU.variables['lat'][:]
EMEP4UK_EURO_lons=inf_EMUK_EU.variables['lon'][:]
EMEP4UK_EURO_i=inf_EMUK_EU.variables['i'][:]
EMEP4UK_EURO_j=inf_EMUK_EU.variables['j'][:]
inf_EMUK_EU.close()


#Create Day of Year array based on EMEP4UK met files
dayofyear=EMEP4UK_time-EMEP4UK_time[0]+0.5
# append 365th day
dayofyear=np.append(dayofyear,dayofyear[-1]+1)

EMEP4UK_LAI = np.array( [ LAIF.CF(dayofyear,EMEP4UK_lats),      \
                          LAIF.DF(dayofyear,EMEP4UK_lats),      \
                          LAIF.NF(dayofyear,EMEP4UK_lats),      \
                          LAIF.BF(dayofyear,EMEP4UK_lats),      \
                          LAIF.TC(dayofyear,EMEP4UK_lats),      \
                          LAIF.MC(dayofyear,EMEP4UK_lats),      \
                          LAIF.RC(dayofyear,EMEP4UK_lats),      \
                          LAIF.SNL(dayofyear,EMEP4UK_lats),     \
                          LAIF.GR(dayofyear,EMEP4UK_lats),      \
                          LAIF.MS(dayofyear,EMEP4UK_lats)      ] \
                        ).transpose(1,0,2,3)        # swap PFT and time dimensions so that dims:
                                                    # [ time, pft, i, j ]


                          # Re insert for IAM types
                          #LAIF.IAM_CR(dayofyear,EMEP4UK_lats),  \
                          #LAIF.IAM_DF(dayofyear,EMEP4UK_lats),  \
                          #LAIF.IAM_MF(dayofyear,EMEP4UK_lats) \

EMEP4UK_CH = np.array( [ [ CHF.CF(EMEP4UK_lats),             \
                           CHF.DF(EMEP4UK_lats),             \
                           np.zeros_like(EMEP4UK_lats)+8.,  \
                           np.zeros_like(EMEP4UK_lats)+15., \
                           np.zeros_like(EMEP4UK_lats)+1.,  \
                           np.zeros_like(EMEP4UK_lats)+2.,  \
                           np.zeros_like(EMEP4UK_lats)+1.,  \
                           np.zeros_like(EMEP4UK_lats)+0.5, \
                           np.zeros_like(EMEP4UK_lats)+0.3, \
                           np.zeros_like(EMEP4UK_lats)+2.  ] \
                         for i in range(len(dayofyear))]  )

                           # Re insert for IAM types
                           #np.zeros_like(EMEP4UK_lats)+1.,  \
                           #np.zeros_like(EMEP4UK_lats)+20., \
                           #np.zeros_like(EMEP4UK_lats)+8. \

EMEP4UK_EURO_LAI = np.array( [ LAIF.CF(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.DF(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.NF(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.BF(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.TC(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.MC(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.RC(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.SNL(dayofyear,EMEP4UK_EURO_lats),     \
                               LAIF.GR(dayofyear,EMEP4UK_EURO_lats),      \
                               LAIF.MS(dayofyear,EMEP4UK_EURO_lats)      ] \
                             ).transpose(1,0,2,3)        # swap PFT and time dimensions so that dims:
                                                         # [ time, pft, i, j ]


                            # Re insert for IAM types
                            #   LAIF.IAM_CR(dayofyear,EMEP4UK_EURO_lats),  \
                            #   LAIF.IAM_DF(dayofyear,EMEP4UK_EURO_lats),  \
                            #   LAIF.IAM_MF(dayofyear,EMEP4UK_EURO_lats)   \


EMEP4UK_EURO_CH = np.array( [ [ CHF.CF(EMEP4UK_EURO_lats),             \
                                CHF.DF(EMEP4UK_EURO_lats),             \
                                np.zeros_like(EMEP4UK_EURO_lats)+8.,  \
                                np.zeros_like(EMEP4UK_EURO_lats)+15., \
                                np.zeros_like(EMEP4UK_EURO_lats)+1.,  \
                                np.zeros_like(EMEP4UK_EURO_lats)+2.,  \
                                np.zeros_like(EMEP4UK_EURO_lats)+1.,  \
                                np.zeros_like(EMEP4UK_EURO_lats)+0.5, \
                                np.zeros_like(EMEP4UK_EURO_lats)+0.3, \
                                np.zeros_like(EMEP4UK_EURO_lats)+2.] \
                              for i in range(len(dayofyear))]  )

                            # Re insert for IAM types
                            #,  \
                            #    np.zeros_like(EMEP4UK_EURO_lats)+1.,  \
                            #    np.zeros_like(EMEP4UK_EURO_lats)+20., \
                            #    np.zeros_like(EMEP4UK_EURO_lats)+8. 

#Write out EMEP4UK data to file
outf=nc.Dataset(EMEP4UK_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('west_east',len(EMEP4UK_i))
outf.createDimension('south_north',len(EMEP4UK_j))
outf.createDimension('pseudo',len(PFT_sh_names))
outf.createDimension('Time',len(dayofyear))

#create and store dimension variables
out_i_EMEP4UK=outf.createVariable('i_EMEP','float32',('west_east',))
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
out_pseudo_EMEP4UK.note=PFT_note
out_pseudo_EMEP4UK[:]=PFT_indices

out_time_EMEP4UK=outf.createVariable('time','i',('Time',))
out_time_EMEP4UK.long_name='Day of Year'
out_time_EMEP4UK.units='days'
out_time_EMEP4UK[:]=dayofyear

# Create and store lai and canht variables
out_lai_EMEP4UK=outf.createVariable('lai','float32',('Time','pseudo','south_north','west_east'))
out_lai_EMEP4UK.long_name='Leaf Area Index'
out_lai_EMEP4UK.units='m^2/m^2'
out_lai_EMEP4UK[:]=EMEP4UK_LAI

out_ch_EMEP4UK=outf.createVariable('canht','float32',('Time','pseudo','south_north','west_east'))
out_ch_EMEP4UK.long_name='Canopy Height'
out_ch_EMEP4UK.units='metres'
out_ch_EMEP4UK[:]=EMEP4UK_CH

#Global Attributes
outf.title='Vegetation climatology for EMEP4UK 5km grid'
outf.note='Based on EMEP model LAI and Canopy Height calculations'
outf.note1='Reference: Simpson et al., 2012, The EMEP MSCW chemical transport model technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()




#Write out EMEP4UK_EURO data to file
outf=nc.Dataset(EMEP_outfile,'w',format='NETCDF4_CLASSIC')

#create dimensions
outf.createDimension('west_east',len(EMEP4UK_EURO_i))
outf.createDimension('south_north',len(EMEP4UK_EURO_j))
outf.createDimension('pseudo',len(PFT_sh_names))
outf.createDimension('Time',len(dayofyear))

#create and store dimension variables
out_i_EMEP4UK_EURO=outf.createVariable('i_EMEP','float32',('west_east',))
out_i_EMEP4UK_EURO.long_name='EMEP grid i coordinate'
out_i_EMEP4UK_EURO.units='unitless'
out_i_EMEP4UK_EURO[:]=EMEP4UK_EURO_i

out_j_EMEP4UK_EURO=outf.createVariable('j_EMEP','float32',('south_north',))
out_j_EMEP4UK_EURO.long_name='EMEP grid j coordinate'
out_j_EMEP4UK_EURO.units='unitless'
out_j_EMEP4UK_EURO[:]=EMEP4UK_EURO_j

out_pseudo_EMEP4UK_EURO=outf.createVariable('pseudo','i',('pseudo',))
out_pseudo_EMEP4UK_EURO.long_name='Land Surface Type'
out_pseudo_EMEP4UK_EURO.units='unitless'
out_pseudo_EMEP4UK_EURO.note=PFT_note
out_pseudo_EMEP4UK_EURO[:]=PFT_indices

out_time_EMEP4UK_EURO=outf.createVariable('time','i',('Time',))
out_time_EMEP4UK_EURO.long_name='Day of Year'
out_time_EMEP4UK_EURO.units='days'
out_time_EMEP4UK_EURO[:]=dayofyear

# Create and store lai and canht variables
out_lai_EMEP4UK_EURO=outf.createVariable('lai','float32',('Time','pseudo','south_north','west_east'))
out_lai_EMEP4UK_EURO.long_name='Leaf Area Index'
out_lai_EMEP4UK_EURO.units='m^2/m^2'
out_lai_EMEP4UK_EURO[:]=EMEP4UK_EURO_LAI

out_ch_EMEP4UK_EURO=outf.createVariable('canht','float32',('Time','pseudo','south_north','west_east'))
out_ch_EMEP4UK_EURO.long_name='Canopy Height'
out_ch_EMEP4UK_EURO.units='metres'
out_ch_EMEP4UK_EURO[:]=EMEP4UK_EURO_CH

#Global Attributes
outf.title='Vegetation climatology for EMEP4UK-Europe 50km grid'
outf.note='Based on EMEP model LAI and Canopy Height calculations'
outf.note1='Reference: Simpson et al., 2012, The EMEP MSCW chemical transport model technical description'
outf.history='Created by Edward Comyn-Platt (March 2015)'

#Close file
outf.close()


