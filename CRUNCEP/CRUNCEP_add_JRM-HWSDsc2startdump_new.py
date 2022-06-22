#!/usr/bin/python2.7
#
# Program: PLOT_STARTDUMPS.py
# Purpose: Plot data from start dumps to check quality
# 
# Edward Comyn-Platt, 2015
#
######################################################
import numpy as np
import netCDF4 as nc

CRUNCEP_dir ='/users/eow/edwcom/CRUNCEP/'

in_SD=CRUNCEP_dir+'JULES_CRUNCEP_MPI_DougDiag.dump.spin10.19700101.0.nc'

SC_file='/users/eow/edwcom/From_JOEY/ROB_/prescribedb.nc'
SC_param_name='carb'       

out_SD=CRUNCEP_dir+'CRUNCEP_startdump_JRM-HWSDSC_20151120.nc'

grid_file = CRUNCEP_dir+'cru_ncep_land.nc'
# read in lat/lons
grinf=nc.Dataset(grid_file,'r')
grlat=grinf.variables['Latitude'][:]
grlon=grinf.variables['Longitude'][:]
grinf.close()


#get SC data to substitute in
inf=nc.Dataset(SC_file)
new_SC=inf.variables[SC_param_name][0,:,:]
inf.close()
sc_lon = np.arange(-179.75,180.,0.5)
sc_lat = np.arange(-89.75,90,0.5)

# loop round each lat lon pair
new_SC_1D=np.zeros_like(grlat)
for i,lat,lon in zip(range(len(grlat)),grlat,grlon):
    x_index=np.where(sc_lon==lon)[0]
    y_index=np.where(sc_lat==lat)[0]
    new_SC_1D[i]=new_SC[y_index,x_index]

badex   = np.where(np.isnan(new_SC_1D))[0]
goodex = np.where(np.isfinite(new_SC_1D))[0]
for pt in badex:
    temppoint = np.argmin(  ((grlat[pt]-grlat[goodex])**2) \
              + ((grlon[pt]-grlon[goodex])**2) )
   
    new_point=goodex[temppoint]
    new_SC_1D[pt]=new_SC_1D[new_point]



# first oti
inf=nc.Dataset(in_SD,'r')
outf=nc.Dataset(out_SD,'w')
for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))
    
for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    if not (str(var) in ['cs']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        outvar[:]=new_SC_1D
    #
#

outf.close()
inf.close()







