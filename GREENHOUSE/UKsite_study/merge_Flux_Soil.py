#!/bin/env python3.5

import netCDF4 as nc

Flux_file='./ConCrop_Flux_data.nc'
SOIL_file='./ConCrop_Soil_data.nc'
fill_value=-9999.
OUTFILE='ConCrop_FullFlux_data.nc'

# Open Outfile
outf=nc.Dataset(OUTFILE,'w')


# First Copy Flux data straight over
Finf=nc.Dataset(Flux_file,'r')
#Dimensions
for dim in Finf.dimensions:
    outf.createDimension(dim,len(Finf.dimensions[dim]))

# Now copy vars and their attributes
copied_vars=[]
for var in Finf.variables:
    copied_vars.append(str(var))
    if '_FillValue' in Finf.variables[var].ncattrs():
        outvar=outf.createVariable(var,'float',Finf.variables[var].dimensions,\
                          fill_value=fill_value )
    else:
        outvar=outf.createVariable(var,'float',Finf.variables[var].dimensions)

    for att in Finf.variables[var].ncattrs():
        if att!='_FillValue':
            outvar.setncattr(att,Finf.variables[var].getncattr(att))

    outvar[:]=Finf.variables[var][:]

for att in Finf.ncattrs():
    outf.setncattr(att,Finf.getncattr(att))

Finf.close()

# now the soil file
Sinf=nc.Dataset(SOIL_file,'r')
# create new dimensions for depth and position
outf.createDimension('z',2)
outf.createDimension('y',2)

outvar=outf.createVariable('y','float32','y')
outvar.long_name='horizontal postion'
outvar.units='-'
outvar[:]=Sinf.variables['x']
outvar=outf.createVariable('z','float32','z')
outvar.long_name='depth'
outvar.units='m'
outvar[:]=Sinf.variables['z']

for var in ['T_soil','SMC']:
    outvar=outf.createVariable(var,'float32',('time','y','z'), \
                                fill_value=fill_value    )
        
    for att in Sinf.variables[var].ncattrs():
        if att!='_FillValue':
            outvar.setncattr(att,Sinf.variables[var].getncattr(att))
    
    outdata=Sinf.variables[var][:]
    outdata[outdata<-9000]=fill_value

    outvar[:]=outdata

var='downward_heat_flux_at_ground_level'
outvar=outf.createVariable(var,'float32',('time','y'), \
                           fill_value=fill_value    )
   
for att in Sinf.variables[var].ncattrs():
    if att!='_FillValue':
        outvar.setncattr(att,Sinf.variables[var].getncattr(att))

outvar[:]=Sinf.variables[var][:]
    
Sinf.close()

outf.close()
        



