#!/usr/bin/python2.7
#
#  
# Edward Comyn-Platt, 2015
#
######################################################
import numpy as np
import netCDF4 as nc

in_SD='/home/users/ecomynplatt/JULES/SC_simulator_points/JULES_vn4.3.1/WFDEI_0p5/dumps/JULES_v42_WFD-EI-GPCC_Zinke_global_DD_newtopo.dump.spin8.19800101.0.nc'

frac_file='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.veg.fracNew.nc'
frac_name='field1391'

LAI_file='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.veg.funcNew3.3.nc'

LAI_name='field1392'
canht_name='field1393'

out_SD='/home/users/ecomynplatt/JULES/SC_simulator_points/JULES_vn4.3.1/WFDEI_0p5/dumps/JULES_WFD-EI-GPCC_Zinke_global_DD_newtopo.19800101.startdump.nc'


# read in frac data:
inf_frac=nc.Dataset(frac_file,'r')
frac_data=inf_frac.variables[frac_name][:]
inf_frac.close()

#read in LAI and canht data
inf_LAI=nc.Dataset(LAI_file,'r')
LAI_data=inf_LAI.variables[LAI_name][:]
canht_data=inf_LAI.variables[canht_name][:]
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

outvar=outf.createVariable('frac','float32',('tile','land'))
outvar[:]=frac_data

outvar=outf.createVariable('lai','float32',('pft','land'))
outvar[:]=LAI_data

outvar=outf.createVariable('canht','float32',('pft','land'))
outvar[:]=canht_data


outf.close()
inf.close()






