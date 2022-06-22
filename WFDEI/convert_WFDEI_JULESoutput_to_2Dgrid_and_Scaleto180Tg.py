#!/usr/bin/python
#
# Python script to convert CRUNCEP JULES fch4_wetl output onto 2D grid with annual ncdf files
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# April 2015
#
# Contains
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
#
gridfile='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
#
data_dir='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/'
#
infile='JULES_v42_WFD-EI-GPCC_Zinke_global_DD_newtopo.monthly_wetl.nc'
#
outfile='JULES_v42_WFD-EI-GPCC_Zinke-SC_Marthews-TI_global_gridded.monthly.nc'
#
area_file='/users/eow/edwcom/WFD_EI/wfdei-land-area.nc'

#params=['fch4_wetl','fwetl']
fch4_conv_factor=(16.01/12.01)*1e-9
fch4_scale_factor= 180./300.    # Based on scaling the 2000 total global emissions to 180 Tg y-1
time_name='time'
fill_value=-999.9
#

#grinf=nc.Dataset(gridfile,'r')
#index=grinf.variables['land_index'][:]
#grinf.close()

#area_inf=nc.Dataset(area_file,'r')
#AREA=area_inf.variables['area'][:].squeeze()
#area_inf.close()

#AREA_1D=np.zeros(67209)
#AREA_1D[index[np.where(index.mask==False)]-1] = AREA[np.where(index.mask==False)]

# Year to scale data to
#YEAR=2000
#
#YEAR_obj = [ nctime.num2date(0,    \
#                units='seconds since '+str(YEAR)+'-01-01 00:00:00', \
#                calendar='standard'), \
#             nctime.num2date(0,    \
#                units='seconds since '+str(YEAR)+'-12-01 00:00:00', \
#                calendar='standard') ]


# Repeat for NTI data
inf=nc.Dataset(data_dir+infile,'r')
outf=nc.Dataset(data_dir+outfile,'w')

#inf_time_obj=nctime.num2date(inf.variables['time'][:], \
#                             units=inf.variables['time'].units, \
#                             calendar='standard')

#SP=np.where(inf_time_obj==YEAR_obj[0])[0]
#EP=np.where(inf_time_obj==YEAR_obj[1])[0]


# Create Dimensions
outf.createDimension('x',len(inf.dimensions['x']))
outf.createDimension('y',len(inf.dimensions['y']))
outf.createDimension('time',len(inf.dimensions['time']))

#copy dimension variables straight over
for param in ['latitude','longitude','time']:
    outvar=outf.createVariable(param, inf.variables[param].dtype, inf.variables[param].dimensions)
    for attr in inf.variables[param].ncattrs():
        outvar.setncattr( str(attr),inf.variables[param].getncattr(str(attr)) )
    outvar[:]=inf.variables[param][:]

#Copy 'fwetl' straight over
param = 'fwetl'
outvar=outf.createVariable(param, inf.variables[param].dtype, inf.variables[param].dimensions)
for attr in inf.variables[param].ncattrs():
    outvar.setncattr( str(attr),inf.variables[param].getncattr(str(attr)) )
outvar[:]=inf.variables[param][:]


#Scale 'fch4_wetl' before copying
param = 'fch4_wetl'
outvar= outf.createVariable(param, inf.variables[param].dtype, inf.variables[param].dimensions)
for attr in  inf.variables[param].ncattrs():
    if str(attr)=='units':
        outvar.setncattr( str(attr),'kg CH4 m^-2 s^-1')
    else:
        outvar.setncattr( str(attr),inf.variables[param].getncattr(str(attr)) )

# Add extra note stating that the JULES output has been rescaled
outvar.setncattr('Correction_Scale_Factor',str(fch4_scale_factor) )
outvar[:]=inf.variables[param][:]*fch4_conv_factor*fch4_scale_factor

outf.title='JULES output: Methane emissions from wetlands'
outf.owner='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.metdata='WFDEI at 0.5 degrees'
outf.SCdata='Zinke/IGBP soil carbon'
outf.topography='Marthews (2014)'
outf.soilproperties='Van Genuchten'
outf.note='fch4_wetl has been scaled to give total annual emissions for 2000  ~180 Tg'

outf.close()

inf.close()





    






    
    
        


        
    
