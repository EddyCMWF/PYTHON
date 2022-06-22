#!/bin/env python2.7

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.path 
import cartopy.crs as ccrs
from cartopy.io import shapereader
from shapely.geometry import Point
from shapely.prepared import prep
import matplotlib.ticker as mticker
import pandas as pd

from PlotTools import CMAPs
from SoilTools import HWSD

START_YEAR, END_YEAR = 1992, 2015
lat_range = [4,7.5]
latdel=0.5
lon_range = [115,119.5]
londel=0.5

in_LCfile_template = '/prj/CLIFFTOP/ESA_CCI_LC/ESACCI-LC-L4-LCCS-Map-300m-P1Y-YYYY-v2.0.7b.nc'
in_LCvar = 'lccs_class'

output_dir = '/users/eow/edwcom/LOMBOK/Regional_Analysis/'

# Read in the Sabah shape:
shpfilename = shapereader.natural_earth(resolution='10m',category='cultural',
                                        name='admin_1_states_provinces')    

shpReader = shapereader.Reader(shpfilename)
Regions = shpReader.records()
MalayIndo_Regions = []
for region in Regions:
    if region.attributes['admin'] in ['Malaysia','Indonesia','Brunei','Philippines']:
        MalayIndo_Regions.append(region.geometry)
    if region.attributes['woe_label']=='Sabah, MY, Malaysia':
        Sabah_Region = region

Sabah_Region_prep = [ prep(polygon) for polygon in Sabah_Region.geometry ]


##########################################################################################
# HWSD data
hwsd_vars_dict = {'T_OC':{'long_name':'Organic Carbon Content','units':'%','data_range':[0,30]}, 
                  'T_SAND':{'long_name':'Sand Content','units':'%','data_range':[0,50]}, 
                  'T_CLAY':{'long_name':'Clay Content','units':'%','data_range':[0,50]}, 
                  'T_SILT':{'long_name':'Silt Content','units':'%','data_range':[0,50]},
                  'TcarbMass':{'long_name':'Carbon Mass','units':'kg m$^{-2}$','data_range':[0,30]},
                  'TcarbMassRef':{'long_name':'Carbon Mass','units':'kg m$^{-2}$','data_range':[0,30]},
                  'TcarbMassMod':{'long_name':'Carbon Mass','units':'kg m$^{-2}$','data_range':[0,30]},
                  'T_BULK_DENSITY':{'long_name':'Bulk Density','units':'kg m$^{-2}$','data_range':[0,1.5]},
                  }

hwsd = HWSD.hwsd(latrangeOut=lat_range,lonrangeOut=lon_range,
                 latresOut='native',lonresOut='native',
                 output_vars=hwsd_vars_dict.keys())
hwsd.doStuff()

# Get Region Mask:
hwsd_points = []
for lat in hwsd.latsOut: hwsd_points = hwsd_points + [ Point((lon,lat)) for lon in hwsd.lonsOut ] 

hwsd_Region_mask = []
for region_polygon in Sabah_Region_prep:
    hwsd_Region_mask.append([region_polygon.contains(point) for point in hwsd_points])

hwsd_Region_mask = np.sum(hwsd_Region_mask,axis=0).astype(bool).reshape([hwsd.nyOut,hwsd.nxOut])

for hwsd_var in hwsd.output_vars:
    # Mask out the non-Sabah points:
    hwsd.outputXRData[hwsd_var].data[~hwsd_Region_mask] = np.nan
    # Give long name and units for plot:
    hwsd.outputXRData[hwsd_var].attrs['long_name'] = hwsd_vars_dict[hwsd_var]['long_name']
    hwsd.outputXRData[hwsd_var].attrs['units'] = hwsd_vars_dict[hwsd_var]['units']
    
    plot_vmin = hwsd_vars_dict[hwsd_var]['data_range'][0]
    plot_vmax = hwsd_vars_dict[hwsd_var]['data_range'][1]
    
    fig, ax = plt.subplots(figsize=(12,6),subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(lon_range+lat_range)
    
    img = hwsd.outputXRData[hwsd_var].plot.pcolormesh(vmin=plot_vmin,vmax=plot_vmax)
    
    ax.add_geometries(MalayIndo_Regions, ccrs.PlateCarree(), edgecolor='k',facecolor='none')
    ax.add_geometries(Sabah_Region.geometry, ccrs.PlateCarree(), edgecolor='k',facecolor='none')
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(),color='k',linestyle='--',draw_labels=True)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_range[0],lon_range[1]+1,londel))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_range[0],lat_range[1]+1,latdel))
    
    fig.suptitle('HWSD - '+hwsd_vars_dict[hwsd_var]['long_name'],fontsize=18)
    fig.savefig(output_dir + 'Sabah_HWSD_'+hwsd_var+'.png',bbox_inches='tight')
    plt.close()


##########################################################################################
# ESA-CCI land cover data
# Read in geo data from 2015 file:
inf_LC = nc.Dataset(in_LCfile_template.replace('YYYY','2015'),'r')

LC_nlats = len(inf_LC.dimensions['lat'])
LC_nlons = len(inf_LC.dimensions['lon'])

LC_latres = 180./LC_nlats
LC_lonres = 360./LC_nlons

LC_yrange = [ int((inf_LC.variables['lat'].valid_max-lat)/LC_latres)
                    for lat in lat_range[::-1] ]
LC_xrange = [ int((lon-inf_LC.variables['lon'].valid_min)/LC_lonres)
                    for lon in lon_range ]

LC_lons = inf_LC.variables['lon'][LC_xrange[0]:LC_xrange[1]]
LC_lats = inf_LC.variables['lat'][LC_yrange[0]:LC_yrange[1]]
inf_LC.close()

# Calculate area of each LC gridcell:
LC_lons_2d, LC_lats_2d = np.meshgrid(LC_lons,LC_lats)
R_earth = 6371.
LC_area_2d = (  (R_earth**2)*np.radians(LC_lonres)
              * (np.sin(np.radians(LC_lats_2d+LC_latres))-np.sin(np.radians(LC_lats_2d)))
              )

# Get Region Mask:
LC_points = [ Point(point) for point in zip(LC_lons_2d.ravel(), LC_lats_2d.ravel()) ]
Region_mask = []
for region_polygon in Sabah_Region_prep:
    Region_mask.append([region_polygon.contains(point) for point in LC_points])
Region_mask = np.sum(Region_mask,axis=0).astype(bool).reshape(LC_lons_2d.shape)

LC_AREAS_dictionary = {}
years = []
for i,year in enumerate(range(START_YEAR,END_YEAR+1)):
    years.append(year)
    inf_LC = nc.Dataset(in_LCfile_template.replace('YYYY',str(year)),'r')
    LC_subset = inf_LC.variables[in_LCvar][LC_yrange[0]:LC_yrange[1],LC_xrange[0]:LC_xrange[1]].data
    inf_LC.close()
    LC_subset = LC_subset.astype(int)
    
    LC_subset = np.ma.masked_array(LC_subset,mask=np.logical_not(Region_mask),fill_value=0)
    LC_subset.data[LC_subset.mask]=LC_subset.fill_value
    
    
    unique_LCs = np.unique(LC_subset.data[LC_subset.mask==False])
    unique_LCs = np.append(unique_LCs,unique_LCs[-1])
    n_uni_LC = len(unique_LCs)
    
    LC_total_areas = [ np.sum(LC_area_2d[(LC_subset.data==unique_LC)&(LC_subset.mask==False)])
                    for unique_LC in unique_LCs ]
    
    LC_index,LC_labels,LC_cmap,LC_norm=CMAPs.ESA_CCI_cmap_reduced(unique_LCs)
    
    for j in range(n_uni_LC-1):
        if LC_labels[j] in LC_AREAS_dictionary.keys():
            LC_AREAS_dictionary[LC_labels[j]].append(LC_total_areas[j])
        else:
            LC_AREAS_dictionary[LC_labels[j]] = [0]*i + [LC_total_areas[j]]
    
    LC_legend_labels = [LC_label+', %4i km$^2$'%(LC_area) 
                        for LC_label,LC_area in zip (LC_labels,LC_total_areas)]
    
    fig, ax = plt.subplots(figsize=(12,6),subplot_kw={'projection': ccrs.PlateCarree()})
    #ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lon_range+lat_range)
    
    img=ax.pcolormesh(LC_lons_2d,LC_lats_2d,LC_subset,cmap=LC_cmap,norm=LC_norm)
    
    ax.add_geometries(MalayIndo_Regions, ccrs.PlateCarree(), edgecolor='k',facecolor='none')
    ax.add_geometries(Sabah_Region.geometry, ccrs.PlateCarree(), edgecolor='k',facecolor='none')
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(),color='k',linestyle='--',draw_labels=True)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_range[0],lon_range[1]+1,londel))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_range[0],lat_range[1]+1,latdel))
    
    tic_locs = (np.array(LC_index)[:-1]+np.array(LC_index)[1:])/2.
    cbar = plt.colorbar(img,ticks=tic_locs)
    cbar.ax.yaxis.set_ticklabels(LC_legend_labels)
    #ax.outline_patch.set_linewidth(2)
    
    fig.suptitle('ESA-CCI Land Cover for Sabah Region '+str(year),fontsize=18)
    fig.savefig(output_dir + 'Sabah_Landcover_'+str(year)+'.png',bbox_inches='tight')
    plt.close()


LC_AREAS_pd = pd.DataFrame(LC_AREAS_dictionary,index=years)
save_LCs = LC_AREAS_pd.columns[LC_AREAS_pd.mean()>100]
LC_AREAS_pd.to_csv('Sabah_Landcover_'+str(START_YEAR)+'_'+str(END_YEAR)+'.csv',
                   float_format=' %10.2f',index_label='Year',
                   columns=save_LCs) # , header=[LC[:10] for LC in save_LCs])







