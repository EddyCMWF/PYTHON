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
import plot_tools as PT

in_SD_oti='/users/eow/edwcom/WFD_EI/start_dumps/WFDEI_startdump_withHWSDSC_ECP20150909.nc'
##in_SD_nti='/users/eow/edwcom/WFD_EI/start_dumps/WFDEI_startdump_withHWSDSC_ECP20150909.nc' 

SC_file='/users/eow/edwcom/From_JOEY/ROB_/prescribedb.nc'
SC_param_name='carb'       

outfile_oti='/users/eow/edwcom/WFD_EI/start_dumps/WFDEI_startdump_JRM-HWSDSC_20151119.nc'
##outfile_nti='/home/users/ecomynplatt/JULES/MPI_tstep_v4.2_global_WFD_EI/dumps/WFDEI_newtopo_startdump_JRM-HWSDSC_ECP20151119.nc'
l_plot_cs_map=True
plotname='/users/eow/edwcom/WFD_EI/SoilCarbonMap_JRM-HWSD.png'

grid_file = '/users/eow/edwcom/WFD_EI/EI-Halfdeg-land-elevation.nc'

# read in lat/lon
grinf=nc.Dataset(grid_file,'r')
grlat=grinf.variables['latitude'][:]
grlon=grinf.variables['longitude'][:]
grinf.close()

#get SC data to substitute in
inf=nc.Dataset(SC_file)
new_SC=inf.variables[SC_param_name][0,:,:]
# Ignore lats and lons in file cos Joey=Dufus
#sc_lon=inf.variables['lon'][:]
#sc_lat=inf.variables['lat'][:]i+0.5
inf.close()

sc_lon = np.arange(-179.75,180.,0.5)
sc_lat = np.arange(-89.75,90,0.5)

if l_plot_cs_map:
    sc_lon_2d,sc_lat_2d = np.meshgrid(sc_lon,sc_lat)
    clevels=[0,1,2.,4.,8.,12,16.,24,32.,48,64]
    colours=['white','palegoldenrod','darkgreen','black']
    PT.plot_map(new_SC,sc_lon_2d,sc_lat_2d,\
                CLEVELS=clevels,COLOURS=colours, \
                INTERPOLATE_COLOURS=True, \
                FILE_PLOT=plotname,
                CBAR_LABEL='$kg$C $m^{-2}$',\
                PLOT_TITLE='Soil Carbon - HWSD (JRM)',\
                CBAR_PAD=0.3)


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


#new_SC_1D[np.isnan(new_SC_1D)]=0

# first oti
inf=nc.Dataset(in_SD_oti,'r')
outf=nc.Dataset(outfile_oti,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    #for ncattr in inf.variables[str(var)].ncattrs():
    #    outvar.setncattr(str(ncattr),inf.variables[str(var)].getncattr(str(ncattr)))
    #
    if not (str(var) in ['cs']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        print 'here'
        outvar[:]=new_SC_1D

outf.close()
inf.close()





