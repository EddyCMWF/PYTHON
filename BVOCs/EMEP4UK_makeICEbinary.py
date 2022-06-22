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

ICE_type=12
WATER_type=11
NO_ICE_types=[0,1,2,3,4,5,6,7,8,9,10,11,13]

EMEP4UK_DIR='/users/eow/edwcom/EMEP/EMEP4UK/'

# NO ICE IN UK SO ONLY DONE FOR EUROPE
#EMEP4UK_infile = EMEP4UK_DIR+'EMEP4UK_LandFrac.nc.noIAM.oneSOIL'
#EMEP4UK_outfile= EMEP4UK_infile+'.binaryICE'

EMEP4UK_EURO_infile = EMEP4UK_DIR+'EMEP4UK_EUROPE_LandFrac.nc.noIAM.oneSOIL'
EMEP4UK_EURO_outfile= EMEP4UK_EURO_infile+'.binaryICE'


inf_EURO=nc.Dataset(EMEP4UK_EURO_infile)
Land_Frac=inf_EURO.variables['Land_Frac'][:]

ICE=Land_Frac[ICE_type,:,:]
WATER=Land_Frac[WATER_type,:,:]
NO_ICE=np.delete(Land_Frac,ICE_type,0)

ICE_index=np.where(ICE>0.)

NEW_Land_Frac=Land_Frac.copy()


##################################################################################
# set cells with ICE>40% to 100% ICE
###################################
index=np.where(ICE>0.4)
# first set all indexed pixels to zero
NEW_Land_Frac[:,index[0],index[1]]=0.
# Then set indexed ice to 1.
NEW_Land_Frac[ICE_type,index[0],index[1]]=1.
##################################################################################



##################################################################################
# set cells with ICE>20% and ICE+Water>50% to 100% ICE
###################################
index=np.where( (ICE>0.2)&(WATER+ICE>0.5) )
# first set all indexed pixels to zero
NEW_Land_Frac[:,index[0],index[1]]=0.
# Then set indexed ice to 1.
NEW_Land_Frac[ICE_type,index[0],index[1]]=1.
##################################################################################


##################################################################################
# set remaining cells with ICE>0% to ratio of other surface types
###################################
index=np.where( (NEW_Land_Frac[ICE_type,:,:]>0.) & (NEW_Land_Frac[ICE_type,:,:]!=1.) )

#factor=np.sum(NO_ICE,axis=0)
factor=1.-NEW_Land_Frac[ICE_type,:,:]

for ltype in NO_ICE_types:
    NEW_Land_Frac[ltype,index[0],index[1]] = NEW_Land_Frac[ltype,index[0],index[1]] /  \
                                             factor[index[0],index[1]]

NEW_Land_Frac[ICE_type,index[0],index[1]] = 0.


##################################################################################

outf=nc.Dataset(EMEP4UK_EURO_outfile,'w')

for dim in inf_EURO.dimensions:
    outf.createDimension(str(dim),len(inf_EURO.dimensions[dim]))


#dimvars=['i_EMEP','j_EMEP','pseudo','LSmask']
for var in inf_EURO.variables:
    outvar=outf.createVariable(str(var),         \
              inf_EURO.variables[var].dtype,     \
              inf_EURO.variables[var].dimensions )
    for att in inf_EURO.variables[var].ncattrs():
        outvar.setncattr(str(att), inf_EURO.variables[var].getncattr(str(att)) )
    if (str(var)=='Land_Frac'):
        outvar[:]=NEW_Land_Frac
    else:
        outvar[:]=inf_EURO.variables[var][:]

for att in inf_EURO.ncattrs():
    if (str(att)=='history'):
        outf.setncattr(str(att),inf_EURO.getncattr(str(att)) + \
              '; Modified May 2015, Ice made binary for JULES requirements')
    else:
        outf.setncattr(str(att),inf_EURO.getncattr(str(att)))


outf.close()

