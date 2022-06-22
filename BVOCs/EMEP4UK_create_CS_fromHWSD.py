#!/usr/bin/env python
##################################################################################
#
# Program: EMEP4UK_create_CS_fromHWSD.py     
# Author: Edward Comyn-Platt, 02/2013
#
#Purpose: To output an Soil Carbon Map from HWSD on the EMEP4UK grid
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys
#import plot_tools as PT
import pylab as plt


hwsd2EMEP4UK_DIR='/users/eow/edwcom/EMEP/hwsd2emep/'
EMEP4UK_hwsd_MU_file=hwsd2EMEP4UK_DIR+'input_data/hwsd_emep4uk_grid.nc'
OUTPUT_file=hwsd2EMEP4UK_DIR+'EMEP4UK_HWSD_SC.nc'
HWSD_LUT='/users/eow/edwcom/data/HWSD/HWSD_DATA.txt'

T_weight=0.3
S_weight=1.-T_weight

########################################################################################################

# read in the relevent parameters from HWSD_LUT
HWSD_inf = open(HWSD_LUT,'r').readlines()
HWSD_hdrs = HWSD_inf.pop(0).replace('"','').split(',')

LUT = { 'mu':[],   \
        'share':[], \
        'soilT':[], \
        'T_OC':[], \
        'T_BD':[], \
        'T_RD':[], \
        'T_gr':[], \
        'S_OC':[], \
        'S_BD':[], \
        'S_RD':[], \
        'S_gr':[], \
        }

mu_index=HWSD_hdrs.index('MU_GLOBAL')
share_index=HWSD_hdrs.index('SHARE')
soilT_index=HWSD_hdrs.index('SU_SYM90')
TOC_index=HWSD_hdrs.index('T_OC')
TBD_index=HWSD_hdrs.index('T_BULK_DENSITY')
TRD_index=HWSD_hdrs.index('T_REF_BULK_DENSITY')
Tgr_index=HWSD_hdrs.index('T_GRAVEL')
SOC_index=HWSD_hdrs.index('S_OC')
SBD_index=HWSD_hdrs.index('S_BULK_DENSITY')
SRD_index=HWSD_hdrs.index('S_REF_BULK_DENSITY')
Sgr_index=HWSD_hdrs.index('S_GRAVEL')

for line in HWSD_inf:
    line_split=line.replace('"','').split(',')
    LUT['mu'].append(float(line_split[mu_index]))
    LUT['share'].append(float(line_split[share_index]))
    LUT['soilT'].append(line_split[soilT_index][:2])
    LUT['T_OC'].append(float(line_split[TOC_index]))
    LUT['T_BD'].append(float(line_split[TBD_index]))
    LUT['T_RD'].append(float(line_split[TRD_index]))
    LUT['T_gr'].append(float(line_split[Tgr_index]))
    LUT['S_OC'].append(float(line_split[SOC_index]))
    LUT['S_BD'].append(float(line_split[SBD_index]))
    LUT['S_RD'].append(float(line_split[SRD_index]))
    LUT['S_gr'].append(float(line_split[Sgr_index]))


# Convert lists to arrays
LUT['mu']=np.array(LUT['mu'])
LUT['share']=np.ma.masked_equal(LUT['share'],-9999)
LUT['soilT']=np.array(LUT['soilT'])
LUT['T_OC']=np.ma.masked_equal(LUT['T_OC'],-9999)
LUT['T_BD']=np.ma.masked_equal(LUT['T_BD'],-9999)
LUT['T_RD']=np.ma.masked_equal(LUT['T_RD'],-9999)
LUT['T_gr']=np.ma.masked_equal(LUT['T_gr'],-9999)
LUT['S_OC']=np.ma.masked_equal(LUT['S_OC'],-9999)
LUT['S_BD']=np.ma.masked_equal(LUT['S_BD'],-9999)
LUT['S_RD']=np.ma.masked_equal(LUT['S_RD'],-9999)
LUT['S_gr']=np.ma.masked_equal(LUT['S_gr'],-9999)


# Calcualte Top and sub-surface soil carbon amounts using the Refernce Density
TCS_RD=(LUT['T_OC']/100.)*LUT['T_RD']*1000.*(1-(0.01*LUT['T_gr']))
SCS_RD=(LUT['S_OC']/100.)*LUT['S_RD']*1000.*(1-(0.01*LUT['S_gr']))

# Calcualte Top and sub-surface soil carbon amounts using the Bulk Density
TCS_BD=(LUT['T_OC']/100.)*LUT['T_BD']*1000.*(1-(0.01*LUT['T_gr']))
SCS_BD=(LUT['S_OC']/100.)*LUT['S_BD']*1000.*(1-(0.01*LUT['S_gr']))

# create LUT elements for CS based on reference density calc
LUT['T_CS']=TCS_RD
LUT['S_CS']=SCS_RD

# replace ay histosols or andosols with the bulk density calculation
index = (LUT['soilT']=='HS')|(LUT['soilT']=='AN')
LUT['T_CS'][index]=TCS_BD[index]
LUT['S_CS'][index]=SCS_BD[index]


############################################################################################
#  Read in the emep4uk hwsd mapping unit data
inf=nc.Dataset(EMEP4UK_hwsd_MU_file,'r')
EMEP4UK_hwsd_MU = inf.variables['mu'][:].squeeze()
# Create list of unique Mapping Unit values
unique_MUs = list( set( list( EMEP4UK_hwsd_MU.data.flatten() ) ) )

# Create array to store CS data in 
EMEP4UK_hwsd_CS = np.zeros_like(EMEP4UK_hwsd_MU.data)+EMEP4UK_hwsd_MU.fill_value

for MU in unique_MUs:
    index=np.where( LUT['mu']==MU )[0]
    if len(index)==1:
        # only one soil type so compute CS from single top and subsurface pair
        EMEP4UK_hwsd_CS[EMEP4UK_hwsd_MU==MU] = \
                (LUT['T_CS'][index]*T_weight) +\
                (LUT['S_CS'][index]*S_weight) 

    elif len(index)>1:
        # more than one type so sum the weighted compoenents 
        #              (divide by 100 as in percent)
        EMEP4UK_hwsd_CS[EMEP4UK_hwsd_MU==MU] = \
                (np.sum(LUT['T_CS'][index] * LUT['share'][index])*T_weight/100) + \
                (np.sum(LUT['S_CS'][index] * LUT['share'][index])*S_weight/100.) 

    #else:
        # No else statement required as all tiles filled with fill_value

EMEP4UK_hwsd_CS=np.ma.masked_equal(EMEP4UK_hwsd_CS,EMEP4UK_hwsd_MU.fill_value)


##############################################################################################
# Write out to 2d file
outf=nc.Dataset(OUTPUT_file,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

copy_vars=['lat','lon']
for var in copy_vars:
    outvar=outf.createVariable(var,'float32',inf.variables[var].dimensions)
    outvar[:]=inf.variables[var][:]

outvar=outf.createVariable('cs','float32',('x','y'),fill_value=EMEP4UK_hwsd_CS.fill_value)
outvar.units='kgC m^-2'
outvar.long_name='Soil Carbon'
outvar[:]=EMEP4UK_hwsd_CS


outf.title='HWSD Soil Carbon for the EMEP4UK grid'
outf.method1='CS = Organic Carbon Fraction * Ref_Bulk_Density'
outf.method2='For histosols and andesols: CS = Organic Carbon Fraction * Bulk_Density'
outf.note1='Top soil fraction = '+str(T_weight)
outf.owner='Edward Comyn-Platt, edwcom@ceh.ac.uk'

outf.close()
inf.close()














