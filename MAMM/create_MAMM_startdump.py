#!/usr/bin/python
#
#
#

import netCDF4 as nc
import numpy as np
from datetime import datetime as dt



global_dump_file='/users/eow/edwcom/MAMM/start_dumps/m045_GLOBAL_20110101_000000_dump_new.nc'

MAMM_domain_file='/users/eow/edwcom/MAMM/domains/JULES_NCEP_CRU_scan_subdomains.dat'
OUT_domain_file='/users/eow/edwcom/MAMM/domains/JULES_NCEP_CRU_scan_subdomain.dat'

output_MAMM_dumpfile='/users/eow/edwcom/MAMM/start_dumps/m045_MAMM_20110101_000000_dump.nc'


#open MAMM domain file and read index, first column of data
dom_index = np.loadtxt(MAMM_domain_file)[:,0]
dom_index = dom_index.astype('int32')

dom_lon  = list(np.loadtxt(MAMM_domain_file)[:,2])
dom_lat  = list(np.loadtxt(MAMM_domain_file)[:,3])

out_dom_f=open(OUT_domain_file,'w')

for lat,lon in zip(dom_lat,dom_lon):
    out_dom_f.write('   '+str(lat)+'    '+str(lon)+'\n')

out_dom_f.close()


inf=nc.Dataset(global_dump_file,'r')

outf=nc.Dataset(output_MAMM_dumpfile,'w')

for dim in inf.dimensions:
    if (str(dim)=='land'):
        outf.createDimension(str(dim),len(dom_index))
    else:
        outf.createDimension(str(dim),len(inf.dimensions[dim]))




for var in inf.variables:
    print str(var)
    in_data  = inf.variables[var][:]
    in_dtype = inf.variables[var].dtype
    in_dims  = inf.variables[var].dimensions 
    if 'land' in in_dims:
        land_loc = in_dims.index('land')
        out_data = in_data.take(dom_index,axis=land_loc)
    else:
        out_data = in_data
    
    outvar=outf.createVariable( str(var), \
                                inf.variables[var].dtype, \
                                inf.variables[var].dimensions )
    for att in inf.variables[var].ncattrs():
        outvar.setncattr(str(att), \
                         inf.variables[var].getncattr(att) )    
    outvar[:]=out_data


for att in inf.ncattrs():
    if (att=='history'):
        outf.setncattr(att,inf.getncattr(att) + \
                       '\n'+dt.now().strftime("%c") + \
                       ': ECP applied Ices mask for JULES applications' )
    else:
        outf.setncattr(att,inf.getncattr(att))

outf.close()

inf.close()





