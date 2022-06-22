#!/bin/env python2.7

#
#
#import shapely 
#import cartopy as cpy
#import geopandas as gpd

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.path 
import cartopy.crs as ccrs
from cartopy.io import shapereader
from shapely.geometry import Point
from shapely.prepared import prep
import matplotlib.ticker as mticker

from PlotTools import CMAPs

lat_range = [4,7.5]
latdel=0.5
lon_range = [115,119.5]
londel=0.5


in_LCfile = '/prj/CLIFFTOP/ESA_CCI_LC/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.nc'
in_LCvar = 'lccs_class'

shpfilename = shapereader.natural_earth(resolution='10m',category='cultural',
                                        name='admin_1_states_provinces')    
       #                                 name='admin_0_countries')

shpReader = shapereader.Reader(shpfilename)
Regions = shpReader.records()

MalayIndo_Regions = []
for region in Regions:
    if region.attributes['admin'] in ['Malaysia','Indonesia','Brunei','Philippines']:
        MalayIndo_Regions.append(region.geometry)
    if region.attributes['woe_label']=='Sabah, MY, Malaysia':
        Sabah_Region = region


inf_LC = nc.Dataset(in_LCfile,'r')

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
LC_lons_2d, LC_lats_2d = np.meshgrid(LC_lons,LC_lats)
R_earth = 6371.
LC_area_2d = (  (R_earth**2)*np.radians(LC_lonres)
              * (np.sin(np.radians(LC_lats_2d+LC_latres))-np.sin(np.radians(LC_lats_2d)))
              )

LC_points = [ Point(point) for point in zip(LC_lons_2d.ravel(), LC_lats_2d.ravel()) ]
Sabah_Region_prep = [ prep(polygon) for polygon in Sabah_Region.geometry ]

Region_mask = []
for region_polygon in Sabah_Region_prep:
    Region_mask.append([region_polygon.contains(point) for point in LC_points])

Region_mask = np.sum(Region_mask,axis=0).astype(bool).reshape(LC_lons_2d.shape)

LC_subset = inf_LC.variables[in_LCvar][LC_yrange[0]:LC_yrange[1],LC_xrange[0]:LC_xrange[1]].data
LC_subset = LC_subset.astype(int)

LC_subset = np.ma.masked_array(LC_subset,mask=np.logical_not(Region_mask),fill_value=0)
LC_subset.data[LC_subset.mask]=LC_subset.fill_value


unique_LCs = np.unique(LC_subset.data[LC_subset.mask==False])
unique_LCs = np.append(unique_LCs,unique_LCs[-1])
#unique_LCs = np.unique(LC_subset).data
n_uni_LC = len(unique_LCs)

LC_total_areas = [ np.sum(LC_area_2d[(LC_subset.data==unique_LC)&(LC_subset.mask==False)])
                for unique_LC in unique_LCs ]

LC_index,LC_labels,LC_cmap,LC_norm=CMAPs.ESA_CCI_cmap_reduced(unique_LCs)

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

fig.savefig('Sabah_Landcover.png',bbox_inches='tight')











#region = next(Regions)
#for att in region.attributes.keys(): 
#    print(att, region.attributes[att])
#
#print(region.attributes['admin'])
