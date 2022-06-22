
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

garr_file=outdir+'HadGEM_O3coeffs_imogen.nc'
ginf=nc.Dataset(garr_file,'r')
ozone_m=ginf.variables['ozone_m'][:]
ozone_c=ginf.variables['ozone_c'][:]
ginf.close()





