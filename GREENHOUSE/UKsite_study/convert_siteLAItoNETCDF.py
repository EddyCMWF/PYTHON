#!/bin/env python3.5

from datetime import datetime as dt
import netCDF4 as nc
import numpy as np

#raw_file='/users/global/albmar/CarboBiocrop/Brattleby/data/Arable_Field09__LAI+height.dat'
DIR='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
raw_file=DIR+'Brattleby_LAI+height.dat'
outfile=DIR+'Brattleby_LAI_Canht_data.nc.new'
template_file=DIR+'Brattleby_LAI_Canht_data.nc'


inlines=open(raw_file,'r').readlines()
headers=inlines.pop(0).replace(' ','').split(',')

headers=['Date','lai','lai_sd','canht']

data_dict={ hdr:[] for hdr in headers }
for line in inlines:
    split=line.replace(' ','').replace('\n','').split(',')
    for hdr,val in zip(headers,split):
        if hdr=='Date':
            data_dict[hdr].append(val)
        else:
            data_dict[hdr].append(float(val))

date_str=data_dict['Date']
date_obj=[dt.strptime(date,'%d-%b-%Y') for date in date_str]
date_num=nc.date2num(date_obj,units='days since 2012-01-01')
del data_dict['Date']

outf=nc.Dataset(outfile,'w')
tempf=nc.Dataset(template_file,'r')

outf.createDimension('date',len(date_num))

outvar=outf.createVariable('date','float32',('date'),fill_value=-9999)
outvar.units='days since 2012-01-01'
outvar[:]=date_num

for var in data_dict:
    tempvar=tempf.variables[var]
    outvar=outf.createVariable(var,'float32',('date'),fill_value=-9999)
    for att in tempvar.ncattrs():
        outvar.setncattr(att,tempvar.getncattr(att))
    outvar[:]=data_dict[var]

outf.close()
tempf.close()


