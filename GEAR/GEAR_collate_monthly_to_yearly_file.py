#!/bin/env python
#
#
import netCDF4 as nc
import numpy as np

in_dir='/prj/uk_rain/grids_generator/OUTPUT/temporal/GB/daily/netCDF/'
out_dir='/users/eow/edwcom/GEAR/'
file_tag='CEH_GEAR_daily_GB_'
data_names=['rainfall_amount','min_dist']
time_name='time'

year=2015
year_str=str(year)

outfile=out_dir+file_tag+year_str+'.nc'

outdata={}

for month in range(1,13):
    month_str='%02d'%month
    infile=in_dir+file_tag+year_str+month_str+'.nc'
    inf=nc.Dataset(infile,'r')
    
    if month==1:
        for data_name in data_names:
            outdata[data_name]=inf.variables[data_name][:]
        outtime=inf.variables[time_name][:]
    else:
        for data_name in data_names:
            outdata[data_name]=np.append(outdata[data_name],inf.variables[data_name][:],axis=0)
        outtime=np.append(outtime,inf.variables[time_name][:],axis=0)
    
    inf.close()

out_t_length=outdata[data_names[0]].shape[0]

# re open January infile for metadata
infile=in_dir+file_tag+year_str+'01.nc'
inf=nc.Dataset(infile,'r')

outf=nc.Dataset(outfile,'w')

# copy over dimensions
for dim in inf.dimensions:
    if str(dim)=='time':
        outf.createDimension(str(dim),out_t_length)
    else:
        outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    for ncattr in inf.variables[str(var)].ncattrs():
        outvar.setncattr(str(ncattr),inf.variables[str(var)].getncattr(str(ncattr)))
    #
    if not str(var) in data_names+[time_name]:
        print str(var)
        outvar[:]=inf.variables[str(var)][:]
    elif str(var) in data_names:
        print 'writing '+str(var)
        outvar[:]=outdata[str(var)]
    elif str(var)==time_name:
        print 'writing '+time_name
        outvar[:]==time_name

for ncattr in inf.ncattrs():
    outf.setncattr(str(ncattr),inf.getncattr(str(ncattr)))


outf.close()
inf.close()


