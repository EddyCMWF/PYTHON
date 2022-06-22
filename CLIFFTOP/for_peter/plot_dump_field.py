#!/bin/env python3

import numpy as np
import netCDF4 as nc
import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value


infile=optional_argparse('-infile','./vn4.8_imogen.spinup_710.dump.1850_zeroDtemp_zeroWP.nc')

plotvar=optional_argparse('-plotvar','cv')

CMAP=optional_argparse('-cmap','RdBu_r')
PALETTE   = plt.cm.get_cmap(name=CMAP,lut=250)

MAP_TYPE=optional_argparse('map_type','Mesh')

grid_file=optional_argparse('-gridfile','./grid_info.nc')
# Read in grid index for converting JULES output to 2D
grinf=nc.Dataset(grid_file,'r')
grindex=grinf.variables['land_index'][:]
lats_2d=grinf.variables['latitude'][:]
lons_2d=grinf.variables['longitude'][:]
Area=grinf.variables['Area'][:]
grinf.close()

inf=nc.Dataset(infile,'r')

invar=inf.variables[plotvar]

if 'land' not in invar.dimensions:
    print('ERROR: No land dimension for variable: '+plotvar)
    quit()

if len(invar.dimensions)>1:
    print('ERROR: No functionality for dealing with "z" dimensions: '+plotvar)


plotdata = np.ma.masked_array(invar[:][grindex],mask=grindex.mask)
print(np.sum(plotdata*Area))

inf.close()


#PTs.plot_map(plotdata,lons_2d,lats_2d,
FIG,AXIS=plt.subplots(ncols=1,nrows=1,figsize=[15,8])
M=Basemap(projection='cyl',resolution='c',ax=AXIS)
xx,yy=M(lons_2d,lats_2d)
if (MAP_TYPE=='Mesh'):
    IMAGE    = M.pcolormesh(xx,yy,plotdata,cmap=PALETTE) #,norm=NORM)
elif (MAP_TYPE=='Contour'):
    IMAGE    = M.contourf(xx,yy,plotdata,CLEVELS,cmap=PALETTE) #,norm=NORM)

AXIS.set_title(plotvar)

M.drawcountries(linewidth=0.4)
M.drawcoastlines(linewidth=0.4)
M.drawmapboundary()

parallels = np.arange(-90.,90.,30.)
M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5,fontsize=12 )

meridians = np.arange(-180.,180.,45.)
M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5,fontsize=12 )

plt.colorbar(IMAGE,orientation='horizontal') 
plt.show()


