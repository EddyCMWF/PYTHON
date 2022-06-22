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
import glob
#
param='fch4_wetl'
param_conv_factor=(16.01/12.01)*1e-9
time_name='time'
fill_value=-999.9


gridfile='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
#
data_dir='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/'
#
files = glob.glob(data_dir+'*.nc')

for i in range(len(files)):
    print i,' - ', files[i][len(data_dir):]

iFILE = input('\nSelect a file number: ')
infile = files[iFILE]
out_tag = files[iFILE][len(data_dir):-15]+param+'.monthly_gridded.'

#
# Open gridfile and read in index
grinf=nc.Dataset(gridfile,'r')
index=grinf.variables['land_index'][:]
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
mask = np.ones_like(index)
#

# Convert lons to 0-360
# create lon conversion index
lon_con_ind= np.append(np.arange(360,720),np.arange(0,360)).flatten()
#
# apply to lons
lons=lons[lon_con_ind]
# add 360 onto all lons less than 0
lons[np.where(lons<0)]=lons[np.where(lons<0)]+360.
#

# Repeat for NTI data
inf       = nc.Dataset(infile,'r')
#var_dict  = inf.variables
#var_names = []
#for var in var_dict:
#    var_names.append(str(var))
#iVAR = 99
#while (iVAR>=0):
#    for i in range(len(var_names)):
#        print i,' - ', var_names[i]
#    iVAR = input('\nSelect a variable number (negative number if no more required): ')
#    if iVAR<0:
#        break
#    scale_factor = input('Enter Scale Factor Conversion: ')


data_in = inf.variables[param][:].squeeze()

# extract time and create time object
time_in   = inf.variables[time_name][:]
time_units= inf.variables[time_name].units
time_obj  = nctime.num2date(time_in, \
                            units=time_units, \
                            calendar='standard')
inf.close()
#

# Apply index to final dimension of data
#   !this should be 2 at the moment but if statement 
#   !useful for other applications
if (len(data_in.shape)==1):
    print '1D input data'
    data_out=data_in[index-1]*mask
elif (len(data_in.shape)==2):
    print '2D input data'
    data_out=data_in[:,index-1]*mask
elif (len(data_in.shape)==3):
    print '3D input data'
    data_out=data_in[:,:,index-1]*mask
#
# Apply conversion factor
data_out.data[data_out.mask==True]=fill_value
data_out.fill_value=fill_value
data_out=data_out*param_conv_factor
#

# Apply lon_con_ind to final dimension of data
#   !this should be 3 at the moment but if statement 
#   !useful for other applications
if (len(data_out.shape)==2):
    print '2D output data'
    data_out=data_out[:,lon_con_ind]
elif (len(data_out.shape)==3):
    print '3D output data'
    data_out=data_out[:,:,lon_con_ind]
elif (len(data_out.shape)==4):
    print '4D output data'
    data_out=data_out[:,:,:,lon_con_ind]

# now NTI 
for year in range(time_obj[0].year,time_obj[-1].year):
    year_str=str(year)
    #
    # get start and end indices
    start_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year)+'-02-01 00:00:00',
          calendar='standard')
    sp=np.where(time_obj==start_time_obj)[0]
    end_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year+1)+'-02-01 00:00:00',
          calendar='standard')
    ep=np.where(time_obj==end_time_obj)[0]

    #
    # Get year slice
    if year==2014:
       out_data=data_out[sp:,:,:]
       out_time=time_in[sp:]
    else:
       out_data=data_out[sp:ep,:,:]
       out_time=time_in[sp:ep]
    #
    # Open output file
    outfile = data_dir+'gridded/'+out_tag+str(year)+'.nc'
    print outfile
    outf=nc.Dataset(outfile,'w')
    #
    # Create Dimensions
    outf.createDimension('time',out_data.shape[0])
    outf.createDimension('latitude',out_data.shape[1])
    outf.createDimension('longitude',out_data.shape[2])
    #
    #
    outvar=outf.createVariable('time','float32',('time'))
    outvar.units=time_units
    outvar[:]=out_time
    #
    outvar=outf.createVariable('latitude','float32',('latitude'))
    outvar.units='degrees_north'
    outvar[:]=lats
    #
    outvar=outf.createVariable('longitude','float32',('longitude'))
    outvar.units='degrees_east'
    outvar[:]=lons
    #
    outvar=outf.createVariable('fch4wetl','float32',('time','latitude','longitude'), \
                                 fill_value=fill_value )
    outvar.units='kgCH4 m-2 s-1'
    outvar[:]=out_data
    #
    outf.title='fch4_wetl from wfdei-gpcc jules-v4.2 with Marthews TI'
    outf.author='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
    outf.close()
#

    

