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
import matplotlib.pyplot as plt
import plot_tools as PT

l_plot_cs_map=True
fillval=-999.

WFDEI_DIR='/users/eow/edwcom/WFD_EI/'

SC_sources = ['Zinke','HWSD','HWSD_NCSCD']

SC_files = ['qrparm.soil_HWSD_cont_cosbyWFDEI.nc',             \
            'soil_HWSD_soilcarbon_landpoints_0p5.safe.nc',     \
            'qrparm.soil_merge_HWSD_NCSCD_cont_cosbyWFDEI.nc'  ]
SC_pnames  = [ 'field1397', 'field1397', 'field1397' ]

grid_file = 'wfdei-land-mask.nc'

# Read grid file
grinf=nc.Dataset(WFDEI_DIR+grid_file,'r')
grindex=grinf.variables['land_index'][:]-1
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
grimask=np.ones_like(grindex)
lons_2d,lats_2d = np.meshgrid(lons,lats)

for i in [0]:  #range(len(SC_sources)):
    # Read in SC data
    inf=nc.Dataset(WFDEI_DIR+SC_files[i],'r')
    soilcarbon    = inf.variables[SC_pnames[i]][:]
    inf.close()
    
    # Regrid and fill value
    soilcarbon_2d=np.ma.masked_array(soilcarbon[grindex],\
                                     mask=grindex.mask,  \
                                     fill_value=fillval  )
    
    if l_plot_cs_map:
        clevels=[0,1,2.,4.,8.,12,16.,24,32.,48.,64.]
        colours=['white','palegoldenrod','darkgreen','black']
        PT.plot_map(soilcarbon_2d,lons_2d,lats_2d,\
                    CLEVELS=clevels,COLOURS=colours, \
                    INTERPOLATE_COLOURS=True, \
                    FILE_PLOT=WFDEI_DIR+'SoilCarbonMap_'+SC_sources[i]+'.png',
                    CBAR_LABEL='$kg$C $m^{-2}$',\
                    PLOT_TITLE='Soil Carbon - '+SC_sources[i],\
                    CBAR_PAD=0.3)
            
    outf=nc.Dataset(WFDEI_DIR+'Soil_Carbon_'+SC_sources[i]+'_2D.nc','w')
    outf.createDimension('lat',grindex.shape[0])
    outf.createDimension('lon',grindex.shape[1])
    outvar=outf.createVariable('cs','float32',('lat','lon'),fill_value=fillval)
    outvar.units='kg C m^-2'
    outvar[:]=soilcarbon_2d.filled()
    
    outf.title='WFDEI soil carbon data on a 2D grid'
    outf.owner='Edward Comyn-Platt, edwcom@ceh.ac.uk'
    outf.source=SC_sources[i]

    outf.close()


