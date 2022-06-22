#!/usr/bin/python2.7
#
#  
# Edward Comyn-Platt, 2015
#
######################################################
import numpy as np
import netCDF4 as nc

in_SD='/home/users/ecomynplatt/JULES/SC_simulator_points/JULES_vn4.3.1/CRUNCEP_0p5/dumps/pre_run_output.dump.spin1.19900101.0.nc'

LAI_file='/group_workspaces/jasmin/jules/data/cru_ncep/0.5deg/ancil/v3.1/data/veg_func_igbp_cruncep_0p5deg_capUM6.6_time.nc'

param_names=['lai','canht']

out_SD='/home/users/ecomynplatt/JULES/SC_simulator_points/JULES_vn4.3.1/CRUNCEP_0p5/dumps/JULES_v431_CRUNCEP_GLOBAL.dump.nc'


# read in frac data:
inf_LAI=nc.Dataset(LAI_file,'r')
param_dict={}
for param in param_names:
    param_dict[param]=inf_LAI.variables[param][:].squeeze()
    print param_dict[param].shape
inf_LAI.close()

inf=nc.Dataset(in_SD,'r')

outf=nc.Dataset(out_SD,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

outf.createDimension('pft',5)
        
for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    outvar[:]=inf.variables[str(var)][:]

for param in param_names:
    outvar=outf.createVariable(param,'float32',('pft','land'))
    outvar[:]=param_dict[param]

outf.close()
inf.close()






