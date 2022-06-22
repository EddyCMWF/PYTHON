#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:29:37 2019

@author: edwcom
"""
import os   #,sys
import numpy as np
import gdal
import cartopy.crs as ccrs
from cartopy.io import shapereader
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from shapely.geometry import Point#, shape, mapping
from shapely.prepared import prep
#from cartopy.feature import ShapelyFeature
from PlotTools import CMAPs

lat_range = [4,7.5]
latdel=0.5
lon_range = [115,119.5]
londel=0.5

CiFOR_projection = ccrs.TransverseMercator(false_easting=500000,central_longitude=117.,
                                     scale_factor=0.9996)
plot_projection = ccrs.PlateCarree()
border_projection = ccrs.PlateCarree()

plot_corners = [ (lon_range[0],lat_range[0]), (lon_range[1],lat_range[0]), 
                 (lon_range[0],lat_range[1]), (lon_range[1],lat_range[1]) ]

plot_corners_CiFOR = [ CiFOR_projection.transform_point( corner[0], corner[1], plot_projection )
                        for corner in plot_corners ] 
plot_Xcorners_CiFor = [ corner[0] for corner in plot_corners_CiFOR ]
plot_Ycorners_CiFor = [ corner[1] for corner in plot_corners_CiFOR ]

plot_extent = ( lon_range[0], lon_range[1], lat_range[0], lat_range[1] )

# Use max/min corener to ensure encapsulating all:
plot_extent_CiFOR = ( min(plot_Xcorners_CiFor), max(plot_Xcorners_CiFor),
                      min(plot_Ycorners_CiFor), max(plot_Ycorners_CiFor) )

YEARS = ['1970s', '1990', '1995', '2000', '2005', '2010', '2015', '2020']
#YEARS = ['1990', '2000', '2010', '2015', '2020']
LC_types = [ '', 'Uncertain', 
             'Scrub', 'Non Forest since 1973', 'Other Non Forest',
             'Regrowth', 'Intact Forest', 'Logged Forest',
             'IOPP', 'ITP' ]
LC_colours = [ 'lightgrey', 'dimgray',
               'orange', 'gold', 'darkgoldenrod',
               'yellowgreen','olivedrab','y',
               'blueviolet', 'plum' ]



############################################################################################
# Read in the Sabah and surreounding region border shapes:
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
############################################################################################

############################################################################################
# Read in the CiFOR TIF data of 1973 Forest Extent:
Forest1973_dir = '/users/eow/edwcom/LOMBOK/CiFOR/data/Forest_area_in_1973/'
Forest_TIFfile = Forest1973_dir+'REGIONBorneo_ForestArea_1973_CIFOR.tif'
Forest_numpy_datafile = Forest1973_dir+'npy_files/Forest1973_data.npy'
Forest_numpy_maskfile = Forest1973_dir+'npy_files/Forest1973_mask.npy'
Forest_numpy_Xvalsfile = Forest1973_dir+'npy_files/Forest1973_Xvals.npy'
Forest_numpy_Yvalsfile = Forest1973_dir+'npy_files/Forest1973_Yvals.npy'

if os.path.isfile(Forest_numpy_datafile):
    Forest1973_data = np.load(Forest_numpy_datafile)
    Forest1973_mask = np.load(Forest_numpy_maskfile)
    #Forest1973_data_masked = np.ma.masked_array(Forest1973_data, mask=Forest1973_mask)
    
    Forest1973_Xvals = np.load(Forest_numpy_Xvalsfile)
    Forest1973_Yvals = np.load(Forest_numpy_Yvalsfile)
    Forest1973_Xvals_2d, Forest1973_Yvals_2d = np.meshgrid(Forest1973_Xvals,Forest1973_Yvals)
else:
    ForestTIF = gdal.Open(Forest_TIFfile)
    For_trans = ForestTIF.GetGeoTransform()
    ForNx = ForestTIF.RasterXSize
    ForXstart = For_trans[0]
    ForXres = For_trans[1]
    ForNy = ForestTIF.RasterYSize
    ForYres = For_trans[5]
    ForYstart = For_trans[3]
    #Full_extent = (For_trans[0], For_trans[0]+ForNx*For_trans[1],
    #               For_trans[3]+ForNy*For_trans[5], For_trans[3] )
    # offset is the smaller of the plot extents minus y TIF extent,
    #  Or zero if less thane zero
    Xoff = int( max( 0, min((plot_extent_CiFOR[0]-ForXstart)/ForXres,
                            (plot_extent_CiFOR[1]-ForXstart)/ForXres) ) )
    Yoff = int( max( 0, min((plot_extent_CiFOR[2]-ForYstart)/ForYres,
                            (plot_extent_CiFOR[3]-ForYstart)/ForYres) ) )
    # Size is the larger of plot extent minus TIF extent, minus Offset
    #  Or Nx-Xoff if greater than that
    Xsize = min(ForNx-Xoff, int( max((plot_extent_CiFOR[0]-ForXstart)/ForXres,
                                     (plot_extent_CiFOR[1]-ForXstart)/ForXres) ) - Xoff )
    Ysize = min(ForNy-Yoff, int( max((plot_extent_CiFOR[2]-ForYstart)/ForYres,
                                     (plot_extent_CiFOR[3]-ForYstart)/ForYres) ) - Yoff )
    #FullForest_data = ForestTIF.ReadAsArray(xoff=0,yoff=0,xsize=ForNx,ysize=ForNy)
    Forest1973_data = ForestTIF.ReadAsArray(xoff=Xoff,yoff=Yoff,xsize=Xsize,ysize=Ysize)
    Forest1973_data_oceanmask = Forest1973_data==0.
    # Dimensions = [ y, x ]
    Forest1973_Xvals = (np.arange(Xoff,Xoff+Xsize)*ForXres)+ForXstart 
    Forest1973_Yvals = (np.arange(Yoff,Yoff+Ysize)*ForYres)+ForYstart
    np.save(Forest_numpy_Xvalsfile, Forest1973_Xvals)
    np.save(Forest_numpy_Yvalsfile, Forest1973_Yvals)
    
    Forest1973_Xvals_2d, Forest1973_Yvals_2d = np.meshgrid(Forest1973_Xvals,Forest1973_Yvals)
    Forest1973_points = plot_projection.transform_points( CiFOR_projection,
                                                         Forest1973_Xvals_2d, Forest1973_Yvals_2d)
    
    SR_ll_corner = CiFOR_projection.transform_point( Sabah_Region.bounds[0],Sabah_Region.bounds[1],
                                                    plot_projection )
    SR_ur_corner = CiFOR_projection.transform_point( Sabah_Region.bounds[2],Sabah_Region.bounds[3],
                                                    plot_projection)
    
    ll_xindex = np.argmin( np.abs(Forest1973_Xvals-SR_ll_corner[0]) )
    ur_xindex = np.argmin( np.abs(Forest1973_Xvals-SR_ur_corner[0]) )
    ll_yindex = np.argmin( np.abs(Forest1973_Yvals-SR_ll_corner[1]) )
    ur_yindex = np.argmin( np.abs(Forest1973_Yvals-SR_ur_corner[1]) )
    
    Forest1973_mask = np.zeros_like(Forest1973_data,dtype=bool)
    
    F1973_lndpts = np.where((Forest1973_data!=0.)&
                            (Forest1973_Xvals_2d>SR_ll_corner[0])&
                            (Forest1973_Yvals_2d>SR_ll_corner[1])&
                            (Forest1973_Xvals_2d<SR_ur_corner[0])&
                            (Forest1973_Yvals_2d<SR_ur_corner[1])
                           )
    
    for j,i in zip(F1973_lndpts[0],F1973_lndpts[1]):
        Forest1973_mask[j,i] = any([region_polygon.contains(Point(Forest1973_points[j,i,:]))
                                    for region_polygon in Sabah_Region_prep])
    # Invert the mask:
    Forest1973_mask = ~Forest1973_mask
    
    # Save numpy files for next time:
    np.save(Forest_numpy_datafile, Forest1973_data)
    np.save(Forest_numpy_maskfile, Forest1973_mask)
    
Forest1973_data_masked = np.ma.masked_array(Forest1973_data, mask=Forest1973_mask)
    
############################################################################################


############################################################################################
#  Read in the CiFOR OP time-series data.
CiFOR_indir = '/users/eow/edwcom/LOMBOK/CiFOR/data/REGBorneo_1973to2016_CIFOR/'
CiFOR_shpfile = CiFOR_indir+\
            'REGBorneo_OriginOfLandConvertedToITPAndIOPPComplexTrajectory_1973to2016_CIFOR.shp'
CiFOR_outdir = CiFOR_indir + 'plots/'
shpReader = shapereader.Reader(CiFOR_shpfile)

LC_label = 'klas'
Records = shpReader.records()
geometry_dictionaries = { year:{ LC:[] for LC in LC_types } for year in YEARS}
for record in Records:
    if ( (str(record.attributes['year']) in YEARS) &
         (str(record.attributes[LC_label]) in LC_types) ):
        geometry_dictionaries[str(record.attributes['year'])][str(record.attributes[LC_label])
                            ].append(record.geometry)
############################################################################################


############################################################################################
# Plot section:
For_index,For_labels,For_cmap,For_norm= CMAPs.CiFOR_1973_cmap()

for year in YEARS:
    print(year)
    #fig, ax = plt.subplots(figsize=(12,6),subplot_kw={'projection': plot_projection})
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111,projection=plot_projection)
    ax.set_extent(lon_range+lat_range)
    
    img=ax.pcolormesh(Forest1973_Xvals_2d, Forest1973_Yvals_2d, Forest1973_data_masked,
                      cmap=For_cmap,norm=For_norm, transform=CiFOR_projection)
    
    handles = []
    labels = []
    for i,LC_type in enumerate(LC_types):
        if len(geometry_dictionaries[year][LC_type])>0:
            print(LC_type,len(geometry_dictionaries[year][LC_type]),LC_colours[i])
            labels.append(LC_type)
            handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=LC_colours[i]))
            ax.add_geometries(geometry_dictionaries[year][LC_type], CiFOR_projection, 
                              #transform = plot_projection,
                              edgecolor=LC_colours[i], facecolor=LC_colours[i]) 
    
    gl = ax.gridlines(crs=plot_projection,color='k',linestyle='--',draw_labels=True)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_range[0],lon_range[1]+1,londel))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_range[0],lat_range[1]+1,latdel))
    
    ax.add_geometries(MalayIndo_Regions, border_projection, edgecolor='k',facecolor='none')
    fig.suptitle('CiFOR Land Cover - '+str(year),fontsize=20)
    fig.legend(handles,labels)
    fig.savefig('CiFOR_LandCover_'+year+'.png',bbox_inches='tight')
    plt.close()


fig.savefig('test.png',bbox_inches='tight')
############################################################################################

