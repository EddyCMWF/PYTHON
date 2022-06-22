#!/bin/env python

import numpy as np
import netCDF4 as nc
import iris
from maths_tools import DateTimeTools as DTTs

indir='/users/eow/garr/Work/Data/CLIFFTOP/OZONE/INPUT/'
outdir='/users/eow/edwcom/CLIFFTOP/OZONE/'

lower_ch4=1285
upper_ch4=2062
ch4_range=upper_ch4-lower_ch4

month_strings=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

lower_files= [ indir+'ch4_'+str(lower_ch4)+'.pm'+month+'.pp' 
                for month in month_strings ]

upper_files= [ indir+'ch4_'+str(upper_ch4)+'.pm'+month+'.pp' 
                for month in month_strings ]


lower_cube=iris.load_cubes(lower_files)[0]
lower_cube.units='ppbv'
upper_cube=iris.load_cubes(upper_files)[0]
upper_cube.units='ppbv'


m_cube=(upper_cube-lower_cube)/ch4_range

c_cube=lower_cube-(m_cube*lower_ch4)

m_array=np.append(m_cube.data[:,:,96:],m_cube.data[:,:,:96],axis=2)
c_array=np.append(c_cube.data[:,:,96:],c_cube.data[:,:,:96],axis=2)

# Create a time object for output
time=DTTs.DTarange_months(DTTs.datetime.datetime(2000,1,1),DTTs.datetime.datetime(2000,12,1))

lats=m_cube.dim_coords[1].points
lons_o=m_cube.dim_coords[2].points
lons=np.append(lons_o[96:]-360,lons_o[:96])

# Write out full data to netCDF
outfname=outdir+'Ozone_gradient_and_intercept_wrtCH4.nc'
outf=nc.Dataset(outfname,'w')

outf.createDimension('tstep',12)
outf.createDimension('x',192)
outf.createDimension('y',145)

time_units='seconds since 2000-01-01 00:00:00'
outvar=outf.createVariable('time','float32',('tstep'))
outvar.units=time_units
outvar.long_name='time at start of integration period'
outvar[:]=nc.date2num(time,units=time_units)

outvar=outf.createVariable('lat','float32',('y'))
outvar.units='Degrees North'
outvar.long_name='Latitude'
outvar[:]=lats

outvar=outf.createVariable('lon','float32',('x'))
outvar.units='Degrees East'
outvar.long_name='Longitude'
outvar[:]=lons

outvar=outf.createVariable('Ozone_m','float32',('tstep','y','x'))
outvar.units='ppbO3 / ppbCH4'
outvar.long_name='Gradient of Ozone wrt Methane (m)'
outvar[:]=m_array

outvar=outf.createVariable('Ozone_c','float32',('tstep','y','x'))
outvar.units='ppbO3'
outvar.long_name='y-axis intercept of Ozone wrt Methane (c)'
outvar[:]=c_array

outf.Title='Parameters for calculating Ozone as a function of Methane in the form y=mx+c'
outf.Purpose='For interactive ozone damage in JULES-IMOGEN simulations with interactive methane'
outf.Theory='Based on Bill Collins suggestion'
outf.Owner='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
outf.DataSource='Bill Collins simulations at 1285ppb and 2062ppb of CH4'

outf.close()


###################################################################
#  THIS IS NOT FINISHED, RESOLUTION NEEDS TO BE CHANGE TO 2.5x3.75
# Write out imogen grid to netCDF
yindex=(lats>=-55)&(lats<=82.5)
lats_im=lats[yindex]
ny_im=len(lats_im)

outfname=outdir+'Ozone_gradient_and_intercept_wrtCH4_imogengrid.nc'
outf=nc.Dataset(outfname,'w')

outf.createDimension('tstep',12)
outf.createDimension('x',192)
outf.createDimension('y',ny_im)

time_units='seconds since 2000-01-01 00:00:00'
outvar=outf.createVariable('time','float32',('tstep'))
outvar.units=time_units
outvar.long_name='time at start of integration period'
outvar[:]=nc.date2num(time,units=time_units)

outvar=outf.createVariable('lat','float32',('y'))
outvar.units='Degrees North'
outvar.long_name='Latitude'
outvar[:]=lats_im

outvar=outf.createVariable('lon','float32',('x'))
outvar.units='Degrees East'
outvar.long_name='Longitude'
outvar[:]=lons

outvar=outf.createVariable('Ozone_m','float32',('tstep','y','x'))
outvar.units='ppbO3 / ppbCH4'
outvar.long_name='Gradient of Ozone wrt Methane (m)'
outvar[:]=m_array[:,yindex,:]

outvar=outf.createVariable('Ozone_c','float32',('tstep','y','x'))
outvar.units='ppbO3'
outvar.long_name='y-axis intercept of Ozone wrt Methane (c)'
outvar[:]=c_array[:,yindex,:]

outf.Title='Parameters for calculating Ozone as a function of Methane in the form y=mx+c on the Imogen Grid'
outf.Purpose='For interactive ozone damage in JULES-IMOGEN simulations with interactive methane'
outf.Theory='Based on Bill Collins suggestion'
outf.Owner='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
outf.DataSource='Bill Collins simulations at 1285ppb and 2062ppb of CH4'

outf.close()




