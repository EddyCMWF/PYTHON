#/bin/env python
# Create a daily depostion file od N based on the farm log file provided

import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import pandas as pd
import datetime as dt
from maths_tools import DateTimeTools as DTT

outfile='/prj/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Mark_Data/Crichton_Management_Log_Data'

# background atmospheric deposition from accimp on WFDEI grid
# units = kg m^-2 s^-1
Background_N_dep = 2.84e-11

# Conversion factor for N units (kgN per ha per day to kgN per m^2 per sec)
N_con_fact = 1e-4 * (1./86400.)

# Managed N and dates
Managed_N    = [ ( dt.datetime(2014,3,4), 36.3*N_con_fact ) ,\
                 ( dt.datetime(2014,3,12),50.6*N_con_fact ) ,\
                 ( dt.datetime(2014,5,21),25.5*N_con_fact ) ,\
                 ( dt.datetime(2014,6,3), 51.7*N_con_fact ) ,\
                 ( dt.datetime(2014,7,4), 51.7*N_con_fact ) ,\
                 ( dt.datetime(2014,7,7), 25.5*N_con_fact ) ,\
                 ( dt.datetime(2015,3,6), 36.3*N_con_fact ) ,\
                 ( dt.datetime(2015,3,12),50.6*N_con_fact ) ,\
                 ( dt.datetime(2015,5,14),36.3*N_con_fact ) ,\
                 ( dt.datetime(2015,6,10),51.7*N_con_fact ) ,\
                 ( dt.datetime(2015,7,10),36.3*N_con_fact ) ,\
                 ] 

# Conver sion factor for Grass Harvest (tons to kg)
GH_con_fact = 1e3
# Managed Grass Harvest and dates
Managed_GH   = [ ( dt.datetime(2014,5,14), dt.datetime(2014,5,16), 149.580*GH_con_fact ) ,\
                 ( dt.datetime(2014,6,28), dt.datetime(2014,6,29),  94.0*GH_con_fact ) ,\
                 ( dt.datetime(2014,8,11), dt.datetime(2014,8,12), 110.000*GH_con_fact ) ,\
                 ( dt.datetime(2015,5,13), dt.datetime(2015,5,14),  62.260*GH_con_fact ) ,\
                 ( dt.datetime(2015,7,1),  dt.datetime(2015,7,2),   90.940*GH_con_fact ) ,\
                ]


# Conver sion factor for Grass Harvest (L per ha to L per m^2)
Herb_con_fact = 1e-4
# Managed Herbicide and dates
Managed_Herb = [ ( dt.datetime(2014,6,10), 2.0*Herb_con_fact ) ,\
                ]

# Managed Cows in field and dates
# ( in_date, out_date, N cows, N bulls )
Managed_cows = [ ( dt.datetime(2014,9,9)  , dt.datetime(2014,9,12),   2, 0 ) , \
                 ( dt.datetime(2014,9,12) , dt.datetime(2014,9,24),  22, 1 ) , \
                 ( dt.datetime(2014,10,10), dt.datetime(2014,11,20), 22, 0 ) , \
                 ]


start_time = dt.datetime(2014,1,1)
end_time   = dt.datetime(2015,8,1)

date_array = np.array(DTT.DTarange( start_time, end_time ) )

time_units='seconds since 2014-01-01 00:00:00'
time_array=nctime.date2num(date_array,units=time_units)

N_deposition = np.zeros_like(date_array,dtype='float32')
GH_flag      = np.zeros_like(date_array,dtype='float32')
GH_rate      = np.zeros_like(date_array,dtype='float32')
Herbicide    = np.zeros_like(date_array,dtype='float32')
cows         = np.zeros_like(date_array,dtype='float32')
bulls        = np.zeros_like(date_array,dtype='float32')

for N_day in Managed_N:
    print N_day[0], np.where(date_array==N_day[0])
    N_deposition[date_array==N_day[0]]+=N_day[1]

for GH_slot in Managed_GH:
    sp=np.where(date_array==GH_slot[0])[0]
    ep=np.where(date_array==GH_slot[1])[0]
    GH_flag[sp:ep+1]+=1
    harvest_rate=GH_slot[2]/ float(ep-sp+1)
    GH_rate[sp:ep]+=harvest_rate

for H_day in Managed_Herb:
    Herbicide[date_array==H_day[0]]+=H_day[1]

for cow_slot in Managed_cows:
    sp=np.where(date_array==cow_slot[0])[0]
    ep=np.where(date_array==cow_slot[1])[0]
    cows[sp:ep+1]+=cow_slot[2]
    bulls[sp:ep+1]+=cow_slot[3]


#write to netcdf
outf=nc.Dataset(outfile+'.nc','w')
outf.createDimension('time',len(date_array))
outf.createDimension('land',1)
outvar=outf.createVariable('time','float32',('time'))
outvar.units=time_units
outvar[:]=time_array

outvar=outf.createVariable('N_addition','float32',('time','land'))
outvar.units="kg N m-2 s-1"
outvar[:]=N_deposition

outvar=outf.createVariable('GH_flag','float32',('time','land'))
outvar.units="-"
outvar[:]=GH_flag

outvar=outf.createVariable('GH_rate','float32',('time','land'))
outvar.units="kg day^-1"
outvar.note='Fresh weight'
outvar[:]=GH_rate

outvar=outf.createVariable('Herbicide','float32',('time','land'))
outvar.units="L m^-2"
outvar[:]=Herbicide

outvar=outf.createVariable('cows','float32',('time','land'))
outvar.units="cow"
outvar[:]=cows

outvar=outf.createVariable('bulls','float32',('time','land'))
outvar.units="bull"
outvar[:]=bulls

outf.title='Digitised Farm Management Log 2014-03-04 to 2015-07-10'
outf.note='Constructed from farmers log file provided by Mark Lee'
outf.owener='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.close()


