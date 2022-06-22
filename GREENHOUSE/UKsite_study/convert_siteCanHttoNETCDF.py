#!/bin/env python3.5

from datetime import datetime as dt
import netCDF4 as nc
import numpy as np
import pandas as pd


DIR='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
raw_file=DIR+'ConCrop_Dynamic_Metadata_ALL.txt'
in_file=DIR+'Brattleby_LAI_Canht_data.nc'
out_file=DIR+'Brattleby_LAI_Canht_data.nc.new'

inlines=open(raw_file,'r').readlines()
headers=inlines.pop(0).replace(' ','').split(',')

headers=['Date','time','canht']

data_dict={ hdr:[] for hdr in headers }
for line in inlines:
    split=line.replace(' ','').replace('\n','').split(',')
    for hdr,val in zip(headers,split):
        if hdr=='canht':
            data_dict[hdr].append(float(val))
        else:
            data_dict[hdr].append(val)

date_str=data_dict['Date']
date_obj=[dt.strptime(date,'%Y-%m-%d') for date in date_str]
date_num=nc.date2num(date_obj,units='days since 2012-01-01')
del data_dict['Date']
del data_dict['time']
panda1=pd.DataFrame(data_dict,index=date_obj)


inf=nc.Dataset(in_file,'r')
indate=nc.num2date(inf.variables['date'][:],
                    units=inf.variables['date'].units)

in_dict={}
for var in ['lai','lai_sd']:
    in_dict[var]=inf.variables[var][:]

panda2=pd.DataFrame(in_dict,index=indate)

oldCHpanda=pd.DataFrame({'canht':inf.variables['canht'][:]},index=indate)
oldCHindex= [tindex not in panda1.index.date for tindex in oldCHpanda.index.date]

oldCHpanda=oldCHpanda[oldCHindex]

panda1=pd.concat([panda1,oldCHpanda]).sort_index()

panda=pd.concat([panda1,panda2],axis=1)

date_obj=[dt.combine(date,dt.min.time()) for date in panda.index.date]
date_num=nc.date2num(date_obj,units='days since 2012-01-01')


outf=nc.Dataset(out_file,'w')

outf.createDimension('date',len(panda.index))

outvar=outf.createVariable('date','float32',('date'),fill_value=-9999)
outvar.units='days since 2012-01-01'
outvar[:]=date_num

for var in panda:
    tempvar=inf.variables[var]
    outvar=outf.createVariable(var,'float32',('date'),fill_value=-9999)
    for att in tempvar.ncattrs():
        outvar.setncattr(att,tempvar.getncattr(att))
    outvar[:]=panda[var].values

outf.close()
inf.close()


