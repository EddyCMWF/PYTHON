#!/bin/env python

import netCDF4 as nc
import sys

def change_time_units(file,new_units,time_name='time'):

    # open file in ammend mode
    amf=nc.Dataset(file,'a')
    
    amvar=amf.variables[time_name]
    
    #read in time as object
    intime_obj=nc.num2date(amvar[:],units=amvar.units)
    # Calculat out time as number from new units
    outtime_num = nc.date2num(intime_obj,units=new_units)

    # ammend the variable
    amvar.setncattr('units',new_units)
    amvar.setncattr('long_name','Time in '+new_units)
    amvar[:] = outtime_num



if __name__ =='__main__':

    if '-time_name' in sys.argv:
        index=sys.argv.index('-time_name')
        temp=sys.argv.pop(index)
        time_name=sys.argv.pop(index)
    else:
        time_name='time'

    file=sys.argv[1]
    new_units=sys.argv[2]
    

    change_time_units(file,new_units,time_name=time_name)

