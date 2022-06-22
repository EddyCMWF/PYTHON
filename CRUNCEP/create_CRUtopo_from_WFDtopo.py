#/usr/bin/python

import netCDF4 as nc
import numpy as np

WFD_topo_file  = '/users/eow/edwcom/WFD_EI/topoidx_WFDEI_0p5_lp_global.nc'

CRU_topo_file  = '/users/eow/edwcom/CRUNCEP/topoidx_CRUNCEP_0p5_lp_global.nc'

CRU_index_file = '/users/eow/edwcom/CRUNCEP/WFDEI_to_CRUNCEP_index.dat'

Cind_inf = open(CRU_index_file,'r')
CRU_index = Cind_inf.readlines()


Winf = nc.Dataset(WFD_topo_file,'r')
Coutf = nc.Dataset(CRU_topo_file,'w')

#only one dimension ('land') so coded to only put in dimension of lenght of CRU index
for dim in Winf.dimensions:
    Coutf.createDimension( str(dim),len(CRU_index) )

for var in Winf.variables:
    outvar= Coutf.createVariable( str(var), \
                                  Winf.variables[var].dtype, \
                                  Winf.variables[var].dimensions )
    
    for att in Winf.variables[var].ncattrs():
        outvar.setncattr(str(att), \
                         Winf.variables[var].getncattr(att) )  
    
    outdata=Winf.variables[var][:]
    outvar[:]=outdata[CRU_index]


Coutf.setncattr('title','Aggregate topographic index data for CRU-NCEP')
Coutf.setncattr('institution','CEH - Wallingford')
Coutf.setncattr('source','Simon Dadson and Toby Marthews, U. Oxford')
Coutf.setncattr('contact', 'E. Comyn-Platt (edwcom@ceh.ac.uk)')
Coutf.setncattr('note','Extracted from WFDEI data provided by T. Marthews')

Coutf.close()


