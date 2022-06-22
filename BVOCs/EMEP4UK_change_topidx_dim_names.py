#!/usr/bin/python2.7

import netCDF4 as nc
import numpy as np

infile='/users/eow/edwcom/EMEP/topo2emep/topidx_emep4uk_grid.nc'
outfile='/users/eow/edwcom/EMEP/topo2emep/topidx_emep4uk_grid_newdims.nc'

inf=nc.Dataset(infile,'r')
outf=nc.Dataset(outfile,'w')

#Create dimensions
outf.createDimension('west_east',len(inf.dimensions['x']))
outf.createDimension('south_north',len(inf.dimensions['y']))

#i_EMEP
outvar=outf.createVariable('i_EMEP','float32',('west_east'))
outvar.longname="EMEP i coordinate"
outvar.units="unitless"
outvar[:]=np.arange(34.4,56.31,0.1)

#j_EMEP
outvar=outf.createVariable('j_EMEP','float32',('south_north'))
outvar.longname="EMEP j coordinate"
outvar.units="unitless"
outvar[:]=np.arange(31.3,58.21,0.1)

#cen_lon
outvar=outf.createVariable('cen_lon','float64',('south_north','west_east'))
outvar.longname="Centre Longitude of Grid-Cell"
outvar.units="Degrees East"
outvar[:]=inf.variables['cen_lon'][:]

#cen_lat
outvar=outf.createVariable('cen_lat','float64',('south_north','west_east'))
outvar.longname="Centre Latitude of Grid-Cell"
outvar.units="Degrees North"
outvar[:]=inf.variables['cen_lat'][:]

#timean
outvar=outf.createVariable('timean','float64',('south_north','west_east'),fill_value=-9999.)
outvar.longname="Topographic index mean"
outvar.units="Unitless"
outvar[:]=inf.variables['timean'][:]

#tistd
outvar=outf.createVariable('tistd','float64',('south_north','west_east'),fill_value=-9999.)
outvar.longname="Topographic index standard deviation"
outvar.units="Unitless"
outvar[:]=inf.variables['timean'][:]

outf.title= "Topographic index data remapped to the EMEP4UK 5km grid" 

outf.close()
inf.close()



infile='/users/eow/edwcom/EMEP/topo2emep/topidx_emep4ukEUROPE_grid.nc'
outfile='/users/eow/edwcom/EMEP/topo2emep/topidx_emep4ukEUROPE_grid_newdims.nc'

inf=nc.Dataset(infile,'r')
outf=nc.Dataset(outfile,'w')

#Create dimensions
outf.createDimension('west_east',len(inf.dimensions['x']))
outf.createDimension('south_north',len(inf.dimensions['y']))

#i_EMEP
outvar=outf.createVariable('i_EMEP','float32',('west_east'))
outvar.longname="EMEP i coordinate"
outvar.units="unitless"
outvar[:]=inf.variables['i_EMEP'][:]

#j_EMEP
outvar=outf.createVariable('j_EMEP','float32',('south_north'))
outvar.longname="EMEP j coordinate"
outvar.units="unitless"
outvar[:]=inf.variables['j_EMEP'][:]

#cen_lon
outvar=outf.createVariable('cen_lon','float64',('south_north','west_east'))
outvar.longname="Centre Longitude of Grid-Cell"
outvar.units="Degrees East"
outvar[:]=inf.variables['cen_lon'][:]

#cen_lat
outvar=outf.createVariable('cen_lat','float64',('south_north','west_east'))
outvar.longname="Centre Latitude of Grid-Cell"
outvar.units="Degrees North"
outvar[:]=inf.variables['cen_lat'][:]

#timean
outvar=outf.createVariable('timean','float64',('south_north','west_east'),fill_value=-9999.)
outvar.longname="Topographic index mean"
outvar.units="Unitless"
outvar[:]=inf.variables['timean'][:]

#tistd
outvar=outf.createVariable('tistd','float64',('south_north','west_east'),fill_value=-9999.)
outvar.longname="Topographic index standard deviation"
outvar.units="Unitless"
outvar[:]=inf.variables['timean'][:]

outf.title= "Topographic index data remapped to the EMEP4UK 5km grid" 

outf.close()
inf.close()
