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
gridfile='/users/eow/edwcom/WFDEI/wfdei-land-mask.nc'
#
data_dir='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/'
#
oti_file_in='JULES_v42_WFD-EI-GPCC_MPI_global_DougDiag.monthly_wetl.nc'
nti_file_in='JULES_v42_WFD-EI-GPCC_MPI_global_DougDiag_newtopo.monthly_wetl.nc'
#
oti_out_tag='JULES_v42_WFD-EI-GPCC_MPI_global_DougDiag_gridded.fch4_wetl.monthly.'
nti_out_tag='JULES_v42_WFD-EI-GPCC_MPI_global_DougDiag_nti_gridded.fch4_wetl.monthly.'
#
#
param='fch4_wetl'
param_conv_factor=(16.01/12.01)*1e-9
time_name='time'
fill_value=-999.9
#
# Open gridfile and read in index
grinf=nc.Dataset(gridfile,'r')
index=grinf.variables['land_index'][:]
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
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
# Open OTI file and read in variable data
OTinf=nc.Dataset(data_dir+oti_file_in,'r')
oti_data_in=OTinf.variables[param][:].squeeze()
oti_time_in=OTinf.variables[time_name][:]
oti_time_units=OTinf.variables[time_name].units
OTinf.close()
#
# Create time object
oti_time_obj=nctime.num2date(oti_time_in, \
                            units=oti_time_units, \
                            calendar='standard')
#
# Apply index to final dimension of data
#   !this should be 2 at the moment but if statement 
#   !useful for other applications
if (len(oti_data_in.shape)==1):
    print '1D input oti data'
    oti_data_out=oti_data_in[index-1]
elif (len(oti_data_in.shape)==2):
    print '2D input oti data'
    oti_data_out=oti_data_in[:,index-1]
elif (len(oti_data_in.shape)==3):
    print '3D input oti data'
    oti_data_out=oti_data_in[:,:,index-1]
#
#
# Apply conversion factor
oti_data_out=oti_data_out*param_conv_factor
#
# Apply index mask to data
oti_data_out[:,index.mask==True]=fill_value
oti_data_out=np.ma.masked_equal(oti_data_out,fill_value)
#
#
# Apply lon_con_ind to final dimension of data
#   !this should be 3 at the moment but if statement 
#   !useful for other applications
if (len(oti_data_out.shape)==2):
    print '2D output oti data'
    oti_data_out=oti_data_out[:,lon_con_ind]
elif (len(oti_data_out.shape)==3):
    print '3D output oti data'
    oti_data_out=oti_data_out[:,:,lon_con_ind]
elif (len(oti_data_out.shape)==4):
    print '4D output oti data'
    oti_data_out=oti_data_out[:,:,:,lon_con_ind]
#
#
# Repeat for NTI data
NTinf=nc.Dataset(data_dir+nti_file_in,'r')
nti_data_in=NTinf.variables[param][:].squeeze()
nti_time_in=NTinf.variables[time_name][:]
nti_time_units=NTinf.variables[time_name].units
NTinf.close()
#
# Create time object
nti_time_obj=nctime.num2date(nti_time_in, \
                            units=nti_time_units, \
                            calendar='standard')
#
# Apply index to final dimension of data
#   !this should be 2 at the moment but if statement 
#   !useful for other applications
if (len(nti_data_in.shape)==1):
    print '1D input nti data'
    nti_data_out=nti_data_in[index-1]
elif (len(nti_data_in.shape)==2):
    print '2D input nti data'
    nti_data_out=nti_data_in[:,index-1]
elif (len(nti_data_in.shape)==3):
    print '2D input nti data'
    nti_data_out=nti_data_in[:,:,index-1]
#
# Apply conversion factor
nti_data_out=nti_data_out*param_conv_factor
#
# Apply index mask to data
nti_data_out[:,index.mask==True]=fill_value
nti_data_out=np.ma.masked_equal(nti_data_out,fill_value)
#
# Apply lon_con_ind to final dimension of data
#   !this should be 3 at the moment but if statement 
#   !useful for other applications
if (len(nti_data_out.shape)==2):
    print '2D output oti data'
    nti_data_out=nti_data_out[:,lon_con_ind]
elif (len(nti_data_out.shape)==3):
    print '3D output oti data'
    nti_data_out=nti_data_out[:,:,lon_con_ind]
elif (len(nti_data_out.shape)==4):
    print '4D output oti data'
    nti_data_out=nti_data_out[:,:,:,lon_con_ind]
#
# now loop round each year of data and produce output
# OTI first
for year in range(oti_time_obj[0].year,oti_time_obj[-1].year):
    year_str=str(year)
    #
    # get start and end indices
    start_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year)+'-02-01 00:00:00',
          calendar='standard')
    sp=np.where(oti_time_obj==start_time_obj)[0]
    end_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year+1)+'-02-01 00:00:00',
          calendar='standard')
    ep=np.where(oti_time_obj==end_time_obj)[0]
    #
    # Get year slice
    if year==2012:
        out_data=oti_data_out[sp:,:,:]
        out_time=oti_time_in[sp:]
    else:
        out_data=oti_data_out[sp:ep,:,:]
        out_time=oti_time_in[sp:ep]
    #
    # Open output file
    outf=nc.Dataset(data_dir+oti_out_tag+str(year)+'.nc','w')
    #
    # Create Dimensions
    outf.createDimension('time',out_data.shape[0])
    outf.createDimension('latitude',out_data.shape[1])
    outf.createDimension('longitude',out_data.shape[2])
    #
    #
    outvar=outf.createVariable('time','float32',('time'))
    outvar.units=oti_time_units
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
    outf.title='fch4_wetl from wfdei-gpcc jules-v4.2 with hydro1k TI'
    outf.author='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
    outf.close()
#
# now NTI 
for year in range(nti_time_obj[0].year,nti_time_obj[-1].year):
    year_str=str(year)
    #
    # get start and end indices
    start_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year)+'-02-01 00:00:00',
          calendar='standard')
    sp=np.where(oti_time_obj==start_time_obj)[0]
    end_time_obj=nctime.num2date(0, 
          units='seconds since '+str(year+1)+'-02-01 00:00:00',
          calendar='standard')
    ep=np.where(nti_time_obj==end_time_obj)[0]
    #
    # Get year slice
    if year==2012:
       out_data=nti_data_out[sp:,:,:]
       out_time=nti_time_in[sp:]
    else:
       out_data=nti_data_out[sp:ep,:,:]
       out_time=nti_time_in[sp:ep]
    #
    # Open output file
    outf=nc.Dataset(data_dir+nti_out_tag+str(year)+'.nc','w')
    #
    # Create Dimensions
    outf.createDimension('time',out_data.shape[0])
    outf.createDimension('latitude',out_data.shape[1])
    outf.createDimension('longitude',out_data.shape[2])
    #
    #
    outvar=outf.createVariable('time','float32',('time'))
    outvar.units=nti_time_units
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

    






    
    
        


        
    
