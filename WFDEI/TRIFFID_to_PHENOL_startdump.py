#!/bin/env python

import netCDF4 as nc
import numpy as np

DUMP_DIR='/users/eow/edwcom/WFD_EI/'

TRIF_file=DUMP_DIR+'J4.6_WFDEI_TRIFFID_VG.dump.spin1.19800101.0.nc'
PHEN_file=DUMP_DIR+'J4.6_WFDEI_PHENOL_VG.dump.spin1.19800101.0.nc'
print(TRIF_file)
inf=nc.Dataset(TRIF_file,'r')
outf=nc.Dataset(PHEN_file,'w')

for dim in inf.dimensions:
    if str(dim)!='scpool':
        outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))
    else:
        outf.createDimension('scpool',1)

for var in inf.variables:
    invar=inf.variables[var]
    outvar=outf.createVariable(str(var),'float32',invar.dimensions)
    if str(var)!='cs':
        outvar[:]=invar[:]
    else:
        outvar[:]=np.sum(invar[:],axis=0)

outf.close()
inf.close()



