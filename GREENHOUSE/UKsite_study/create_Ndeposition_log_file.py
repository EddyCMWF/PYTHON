#/bin/env python
# Create a daily depostion file od N based on the farm log file provided

import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import pandas as pd
import datetime as dt

template_file='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/output/Crichton/Jvn4.3.1-E-F/Jvn4.3.1-E-F_crich.day.nc'
outfile='/prj/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Mark_Data/Crichton_Managed_plus_AtmosBG_Ndeposition'

# background atmospheric deposition from accimp on WFDEI grid
# units = kg m^-2 s^-1
Background_N_dep = 2.84e-11

N_con_fact = 1e-4 * (1./86400.)

# Managed N and dates
Managed_N_dt = np.array( [ \
                           dt.date(2014,3,4),\
                           dt.date(2014,3,12),\
                           dt.date(2014,5,21),\
                           dt.date(2014,6,3),\
                           dt.date(2014,7,4),\
                           dt.date(2014,7,7),\
                           dt.date(2015,3,6),\
                           dt.date(2015,3,12),\
                           dt.date(2015,5,14),\
                           dt.date(2015,6,10),\
                           dt.date(2015,7,10),\
                           dt.date(2015,7,13),\
                           dt.date(2016,3,6),\
                           dt.date(2016,3,12),\
                           ] , dtype='datetime64[D]'  )

Managed_N_qt = np.array( [ 36.3*N_con_fact , \
                           50.6*N_con_fact , \
                           25.5*N_con_fact , \
                           51.7*N_con_fact , \
                           51.7*N_con_fact , \
                           25.5*N_con_fact , \
                           36.3*N_con_fact , \
                           50.6*N_con_fact , \
                           36.3*N_con_fact , \
                           51.7*N_con_fact , \
                           36.3*N_con_fact , \
                           25.5*N_con_fact , \
                           36.3*N_con_fact , \
                           50.6*N_con_fact , \
                           ]  )

# Get time vector from template file
inf=nc.Dataset(template_file,'r')
start_time=nctime.num2date( inf.variables['time'][0], \
                          units=inf.variables['time'].units )
end_time=nctime.num2date( inf.variables['time'][-1]+86400, \
                          units=inf.variables['time'].units )
time_array=inf.variables['time'][:] 
time_units=inf.variables['time'].units
inf.close()

date_array=np.arange(start_time.date(),end_time.date(),dtype='datetime64')

N_deposition=np.zeros_like(date_array,dtype='float32')+Background_N_dep

for N_day,N_Q in zip(Managed_N_dt,Managed_N_qt):
    N_deposition[date_array==N_day]+=N_Q

#write to netcdf
outf=nc.Dataset(outfile+'.nc','w')
outf.createDimension('time',len(date_array))
outf.createDimension('land',1)
outvar=outf.createVariable('time','float32',('time'))
outvar.units=time_units
outvar[:]=time_array
outvar=outf.createVariable('N_dep','float32',('time','land'))
outvar.units="kg N m-2 s-1"
outvar[:]=N_deposition
outf.title='Total N deposition (Atmospheric+Farm Managed) for the Crichton Field Site'
outf.note='Constructed from farmers log file provided by Mark Lee'
outf.owener='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.close()


#write to ascii
outf=open(outfile+'.dat','w')
outf.write('# Total N deposition (Atmospheric+Farm Managed) for the Crichton Field Site\n')
outf.write('# Daily temporal resolution, starting 2014-12-04, ending 2016-03-30\n')
for N_dep in N_deposition:
    outf.write('%15.4e\n'%(N_dep))
outf.close()

