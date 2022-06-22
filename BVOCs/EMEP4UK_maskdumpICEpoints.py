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
from datetime import datetime as dt
#import plot_tools as PT
#import pylab as plt

ICE_type=12

LandFrac_file= '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_LandFrac.nc.noIAM.oneSOIL.binaryICE'

cruncep_dump_infile  = '/users/eow/edwcom/EMEP/cruncep2emep/cruncep_emep4ukEUROPE_startdump_oneSOIL.nc'
cruncep_dump_outfile = cruncep_dump_infile+'.ICEmasked'

wfdei_dump_infile  = '/users/eow/edwcom/EMEP/wfd2emep/wfdei_emep4ukEUROPE_startdump_oneSOIL.nc'
wfdei_dump_outfile = wfdei_dump_infile+'.ICEmasked'

infiles=[wfdei_dump_infile,cruncep_dump_infile]
outfiles=[wfdei_dump_outfile,cruncep_dump_outfile]


LF_inf=nc.Dataset(LandFrac_file)
ICE_mask=LF_inf.variables['Land_Frac'][ICE_type,:,:]
LF_inf.close()

for infile,outfile in zip(infiles,outfiles):
    inf=nc.Dataset(infile,'r')
    outf=nc.Dataset(outfile,'w')
    
    for dim in inf.dimensions:
        outf.createDimension(str(dim),len(inf.dimensions[dim]))
    
    for var in inf.variables:
        outvar=outf.createVariable( str(var), \
                                    inf.variables[var].dtype, \
                                    inf.variables[var].dimensions )
       
        for att in inf.variables[var].ncattrs():
            outvar.setncattr(str(att), \
                             inf.variables[var].getncattr(att) )
        
        indata=inf.variables[var][:]
        outdata=indata.copy()
        
        if (var=='vsat'):
            for i in range(indata.shape[0]):
                tempdata=indata[i,:,:]
                tempdata[ICE_mask==1]=0.
                outdata[i,:,:]=tempdata
        
        outvar[:]=outdata
    
    for att in inf.ncattrs():
        if (att=='history'):
            outf.setncattr(att,inf.getncattr(att) + \
                           '\n'+dt.now().strftime("%c") + \
                           ': ECP applied Ices mask for JULES applications' )
        else:
            outf.setncattr(att,inf.getncattr(att))
    
    outf.close()




