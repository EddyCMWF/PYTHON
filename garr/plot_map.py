#
# Python module to selected datasets
#
# Garry Hayman
# Centre for Ecology and Hydrology
# December 2011
#
# Based on plot_maps.py of Simon Dadson
#
from matplotlib import pylab, rc, rcParams, cm
from mpl_toolkits.basemap import Basemap, cm, shiftgrid, addcyclic
import numpy as np
import copy
import matplotlib.colors as col
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
#
def plot_map(DATA_IMAGE,START_LONG,END_LONG,DLONG, \
	START_LAT,END_LAT,DLAT,DATA_MAX,DATA_MIN, \
	PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG):
#
	FIG     = plt.figure(figsize=(12.0,8.0))
#
	FONT    = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 10 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	LONG    = np.arange(START_LONG,END_LONG+DLONG,DLONG)
	LAT     = np.arange(START_LAT, END_LAT+DLAT, DLAT )
#
# Create figure
#
# Create polar stereographic Basemap instance.
#
	M        = Basemap(projection='cyl', \
			llcrnrlat=START_LAT,urcrnrlat=END_LAT, \
			llcrnrlon =START_LONG,urcrnrlon=END_LONG,fix_aspect=False)
	AX       = plt.gca()
#
# Set up a colormap:
#
#	PALETTE  = cm.gist_rainbow
#	PALETTE  = cm.spectral
	PALETTE  = cm.jet
#
	PALETTE.set_under(SET_UNDER, 1.0)
	PALETTE.set_bad('black', 1.0)
	PALETTE.set_over('black', 1.0)
#
	IMAGE    = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest', \
		norm=Normalize(vmin=DATA_MIN,vmax=DATA_MAX))
#
	print(PLOT_TITLE)
	plt.title(PLOT_TITLE,fontsize=14)
#
# draw coastlines.
#
	M.drawrivers(linewidth=0.1,color='grey')
	M.drawcountries(linewidth =0.2)
	M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
	M.drawmapboundary()
#
# draw parallels and meridians.
#
	M.drawparallels(LAT, labels=[1,0,0,0], linewidth = 0.2)
	M.drawmeridians(LONG,labels=[0,0,0,1], linewidth = 0.2)
#
#	plt.xlabel('Longitude',fontsize=12)
#	plt.ylabel('Latitude',fontsize=12)
#
# New axis for colorbar.
#
#	AXIS     = inset_axes(AX,
#		width="5%", # width = 10% of parent_bbox width
#		height="50%", # height : 50%
#		loc=3)
#
	DIVIDER  = make_axes_locatable(AX)
#
	AXIS     = DIVIDER.append_axes('bottom', size='4%', pad=0.75)
#	COLORBAR = plt.colorbar(IMAGE,cax=AXIS,extend='min',orientation='horizontal')
	COLORBAR = plt.colorbar(IMAGE,cax=AXIS,orientation='horizontal')
	COLORBAR.set_label(PLOT_LABEL)
#
#	plt.tight_layout()
#
	if iDISPLAY=='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map2(DATA_IMAGE,START_LONG,END_LONG,LONG,LONG_MAP, \
	START_LAT,END_LAT,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,\
	PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,SET_OVER,DEBUG):
#
# Create figure
#
	FIG      = plt.figure(figsize=(12.0,8.0))
#
	FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 10 }
	matplotlib.rc('font', **FONT) 
	rc('text',usetex=True)
#
# Create polar stereographic Basemap instance.
#
	M         = Basemap(projection='cyl', \
			llcrnrlat =START_LAT,urcrnrlat=END_LAT, \
			llcrnrlon =START_LONG,urcrnrlon=END_LONG,fix_aspect=False)
	AX        = plt.gca()
#
# Set up a colormap:
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Contour0' or MAP_TYPE.find('JULES'):
#
		NCOLOURS  = len(CLEVELS)-1
#
# Set up colour map as used in GRADS
#
#		COLOURS   = [ '#a000c8', '#6e00dc', '#1e3cff', '#00a0ff', '#00c8c8',
#			      '#00d28c', '#00dc00', '#a0e632', '#e6dc32', '#e6af2d',
#			      '#f08228', '#fa3c3c', '#f00082' ]
#
		if len(COLOURS) == 0:
			COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]
#
		COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
		cm.register_cmap(cmap=COL_MAP)
		PALETTE   = COL_MAP
#
		NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
#
	elif MAP_TYPE == 'Contour1':
		PALETTE  = cm.gist_rainbow
	elif MAP_TYPE == 'Contour2':
		PALETTE  = cm.spectral
	elif MAP_TYPE == 'Contour3':
		PALETTE  = cm.jet
#
	PALETTE.set_under(SET_UNDER, 1.0)
	PALETTE.set_over(SET_OVER,  1.0)
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Map_JULES':
		IMAGE    = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)
	elif MAP_TYPE == 'Contour0' or MAP_TYPE == 'Contour0_JULES':
		IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE,norm=NORM)
	else:
		IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE)
#
	plt.title(PLOT_TITLE,fontsize =14)
#
# draw coastlines.
#
#	M.drawrivers(linewidth =0.1,color='grey')
	M.drawcountries(linewidth =0.2)
	M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
	M.drawmapboundary()
#
# draw parallels and meridians.
#
	M.drawparallels(LAT, labels =[1,0,0,0], linewidth = 0.2)
	M.drawmeridians(LONG,labels =[0,0,0,1], linewidth = 0.2)
#
	DIVIDER   = make_axes_locatable(AX)
#
	AXIS      = DIVIDER.append_axes('bottom', size='4%', pad=0.75)
	COLORBAR  = plt.colorbar(IMAGE,cax=AXIS,extend='both',orientation='horizontal',ticks=CLEVELS)
	COLORBAR.set_label(PLOT_LABEL)
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map3(DATA_IMAGE,START_LONG,END_LONG,LONG,LONG_MAP, \
	START_LAT,END_LAT,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH,HEIGHT,ASPECT, \
	RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES,SET_UNDER,SET_OVER,DEBUG):
#
# Create figure
#
	plt.figure(figsize=(WIDTH,HEIGHT))
#
	FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : FONTSIZES[0] }
	matplotlib.rc('font', **FONT) 
	rc('text',usetex=True)
#
# Create polar stereographic Basemap instance.
#
	M         = Basemap(projection='cyl', \
			llcrnrlat =START_LAT,urcrnrlat=END_LAT, \
			llcrnrlon =START_LONG,urcrnrlon=END_LONG, \
			resolution=RESOLUTION,fix_aspect=ASPECT)
	AX        = plt.gca()
#
# Set up a colormap:
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Contour0' or MAP_TYPE == 'Mesh' or MAP_TYPE.find('JULES'):
#
		NCOLOURS  = len(CLEVELS)-1
#
# Set up colour map as used in GRADS
#
#		COLOURS   = [ '#a000c8', '#6e00dc', '#1e3cff', '#00a0ff', '#00c8c8',
#			      '#00d28c', '#00dc00', '#a0e632', '#e6dc32', '#e6af2d',
#			      '#f08228', '#fa3c3c', '#f00082' ]
#
		if len(COLOURS) == 0:
			COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]
#
		COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
		cm.register_cmap(cmap=COL_MAP)
		PALETTE   = COL_MAP
#
		NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
#
	elif MAP_TYPE == 'Contour1':
		PALETTE  = cm.gist_rainbow
	elif MAP_TYPE == 'Contour2':
		PALETTE  = cm.spectral
	elif MAP_TYPE == 'Contour3':
		PALETTE  = cm.jet
#
	PALETTE.set_under(SET_UNDER, 1.0)
	PALETTE.set_over(SET_OVER,  1.0)
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Map_JULES':
		IMAGE    = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)
	elif MAP_TYPE == 'Mesh':
		IMAGE    = M.pcolormesh(LONG_MAP,LAT_MAP,DATA_IMAGE,cmap=PALETTE)
	elif MAP_TYPE == 'Contour0' or MAP_TYPE == 'Contour0_JULES':
		IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE,norm=NORM)
	else:
		IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE)
#
	plt.title(PLOT_TITLE,fontsize = FONTSIZES[3])
#
# draw coastlines.
#
#	M.drawrivers(linewidth =0.1,color='grey')
	M.drawcountries(linewidth =0.2)
	M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
	M.drawmapboundary()
#
# draw parallels and meridians.
#
	M.drawparallels(LAT, labels =[1,0,0,0], linewidth = 0.2)
	M.drawmeridians(LONG,labels =[0,0,0,1], linewidth = 0.2)
#
	DIVIDER   = make_axes_locatable(AX)
#
	AXIS      = DIVIDER.append_axes('bottom', size='4%', pad=0.75)
#	if NCOLOURS > 15:
#		COL_LABS  = []
#		for iCOL in range(NCOLOURS):
#			if iCOL == iCOL/2:
#				COL_LABS.append(CLEVELS[iCOL])
#			else:
#				COL_LABS.append(float('nan'))
#		COL_LABS  = np.array(COL_LABS)
#			
#	COLORBAR  = plt.colorbar(IMAGE,cax=AXIS,extend='both',orientation='horizontal',ticks=COL_LABS)
	COLORBAR  = plt.colorbar(IMAGE,cax=AXIS,extend='both',orientation='horizontal',ticks=CLEVELS)
	COLORBAR.set_label(PLOT_LABEL)
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map3_sites(START_LONG,END_LONG,LONG,LONG_MAP, \
	START_LAT,END_LAT,LAT,LAT_MAP,WIDTH,HEIGHT,ASPECT,RESOLUTION,XOFFSET,YOFFSET, \
	PLOT_LABELS,PLOT_TITLE,PLOT_CODE,FILE_PLOT,iDISPLAY,FONTSIZES,DEBUG):
#
# Create figure
#
	plt.figure(figsize=(WIDTH,HEIGHT))
#
	FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : FONTSIZES[0] }
	matplotlib.rc('font', **FONT) 
	rc('text',usetex=True)
#
# Create polar stereographic Basemap instance.
#
	M         = Basemap(projection='cyl', \
			llcrnrlat =START_LAT,urcrnrlat=END_LAT, \
			llcrnrlon =START_LONG,urcrnrlon=END_LONG, \
			resolution=RESOLUTION,fix_aspect=ASPECT)
#
	for LABEL,X,Y in zip(PLOT_LABELS,LONG_MAP,LAT_MAP):
		plt.text(X+XOFFSET,Y+YOFFSET,LABEL.upper(),fontsize = FONTSIZES[1])

	plt.plot(LONG_MAP,LAT_MAP,PLOT_CODE,markersize=4)
	plt.title(PLOT_TITLE,fontsize = FONTSIZES[3])
#
# draw coastlines
#
#	M.drawrivers(linewidth =0.1,color='grey')
	M.drawcountries(linewidth =0.2)
	M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
	M.drawmapboundary()
#
# draw parallels and meridians.
#
	M.drawparallels(LAT, labels =[1,0,0,0], linewidth = 0.2)
	M.drawmeridians(LONG,labels =[0,0,0,1], linewidth = 0.2)
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map3_GE(DATA_IMAGE,LONG_MAP,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,SIZE, \
	FILE_PLOT,iDISPLAY,FONTSIZES,SET_UNDER,SET_OVER,DEBUG):
#
# Create figure
#
	plt.figure(figsize=(SIZE,round(SIZE*float(DATA_IMAGE.shape[0])/float(DATA_IMAGE.shape[1]),2)))
	plt.axes([0,0,1,1]) 
	plt.axis('off')
#
# Set up a colormap:
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Contour0' or MAP_TYPE.find('JULES'):
#
		NCOLOURS  = len(CLEVELS)-1
#
# Set up colour map as used in GRADS
#
#		COLOURS   = [ '#a000c8', '#6e00dc', '#1e3cff', '#00a0ff', '#00c8c8',
#			      '#00d28c', '#00dc00', '#a0e632', '#e6dc32', '#e6af2d',
#			      '#f08228', '#fa3c3c', '#f00082' ]
#
		if len(COLOURS) == 0:
			COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]
#
		COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
		cm.register_cmap(cmap=COL_MAP)
		PALETTE   = COL_MAP
#
		NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
#
	elif MAP_TYPE == 'Contour1':
		PALETTE  = cm.gist_rainbow
	elif MAP_TYPE == 'Contour2':
		PALETTE  = cm.spectral
	elif MAP_TYPE == 'Contour3':
		PALETTE  = cm.jet
#
#	PALETTE.set_under(SET_UNDER, 1.0)
#	PALETTE.set_over(SET_OVER,  1.0)
#
	if MAP_TYPE == 'Map' or MAP_TYPE == 'Map_JULES':
		DATA_IMAGE=pylab.flipud(DATA_IMAGE)
		plt.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)
	elif MAP_TYPE == 'Contour0' or MAP_TYPE == 'Contour0_JULES':
		plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE,norm=NORM)
	else:
		plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE)
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT,dpi=300)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map3_multi(NROWS,NCOLUMNS,DATA_IMAGE_ALL,START_LONG,END_LONG,LONG,LONG_MAP, \
	START_LAT,END_LAT,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPE,WIDTH,HEIGHT,ASPECT, \
	RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,iDISPLAY,FONTSIZES,SET_UNDER,SET_OVER,DEBUG):
#
# Create figure
#
	plt.figure(figsize=(WIDTH,HEIGHT))
#
	FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('text',usetex=True)
	matplotlib.rc('font', **FONT) 
#
	NPLOTS      = NROWS*NCOLUMNS
#
	for iPLOT in range(NPLOTS):
#
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
# Create polar stereographic Basemap instance.
#
		M         = Basemap(projection='cyl', \
				llcrnrlat =START_LAT,urcrnrlat=END_LAT, \
				llcrnrlon =START_LONG,urcrnrlon=END_LONG, \
#				rsphere=6371007.181,ellps='WGS84', \
				resolution=RESOLUTION,fix_aspect=ASPECT)
		AX        = plt.gca()
#
		COLOURS     = COLOURS_ALL[iPLOT]
		CLEVELS     = CLEVELS_ALL[iPLOT]
		DATA_IMAGE  = DATA_IMAGE_ALL[iPLOT]
		PLOT_LABEL  = PLOT_LABELS[iPLOT]
		SUB_TITLE   = SUB_TITLES[iPLOT]
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
#
# Set up a colormap:
#
		if MAP_TYPE[iPLOT] == 'Map' or MAP_TYPE[iPLOT] == 'Contour0' or MAP_TYPE[iPLOT].find('JULES'):
#
			NCOLOURS  = len(CLEVELS)-1
#
			if len(COLOURS) == 0:
				COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]
#
			COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
			cm.register_cmap(cmap=COL_MAP)
			PALETTE   = COL_MAP
#
			NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
#
		elif MAP_TYPE[iPLOT] == 'Contour1':
			PALETTE  = cm.gist_rainbow
		elif MAP_TYPE[iPLOT] == 'Contour2':
			PALETTE  = cm.spectral
		elif MAP_TYPE[iPLOT] == 'Contour3':
			PALETTE  = cm.jet
#
		PALETTE.set_under(SET_UNDER[iPLOT], 1.0)
		PALETTE.set_over(SET_OVER[iPLOT],  1.0)
#
		if MAP_TYPE[iPLOT] == 'Map' or MAP_TYPE[iPLOT] == 'Map_JULES':
			IMAGE    = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)
		elif MAP_TYPE[iPLOT] == 'Contour0' or MAP_TYPE[iPLOT] == 'Contour0_JULES':
			IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE,norm=NORM)
		else:
			IMAGE    = plt.contourf(LONG_MAP,LAT_MAP,DATA_IMAGE,CLEVELS,cmap=PALETTE)
#
# draw coastlines.
#
#		M.drawrivers(linewidth =0.1,color='grey')
		M.drawcountries(linewidth =0.2)
		M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
		M.drawmapboundary()
#
# draw parallels and meridians.
#
		M.drawparallels(LAT, labels =[1,0,0,0], linewidth = 0.2)
		M.drawmeridians(LONG,labels =[0,0,0,1], linewidth = 0.2)
#
		AX.set_xticklabels(AX.get_xticks(),FONT)
		AX.set_yticklabels(AX.get_yticks(),FONT)
#
		DIVIDER   = make_axes_locatable(AX)
		AXIS      = DIVIDER.append_axes('bottom', size='4%', pad=0.25)
		COLORBAR  = plt.colorbar(IMAGE,cax=AXIS,extend='both',orientation='horizontal',ticks=CLEVELS)
		COLORBAR.set_label(PLOT_LABEL)
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def plot_map4_multi(NROWS,NCOLUMNS,DATA_IMAGE_ALL,START_LONG,END_LONG,LONG,LONG_MAP, \
	START_LAT,END_LAT,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPE,WIDTH,HEIGHT,ASPECT, \
	RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,iDISPLAY,FONTSIZES,SET_UNDER,SET_OVER,DEBUG):
#
# Create figure
#
	plt.figure(figsize=(WIDTH,HEIGHT))
#
	FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('text',usetex=True)
	matplotlib.rc('font', **FONT) 
#
	NPLOTS      = NROWS*NCOLUMNS

	for iPLOT in range(NPLOTS):
#
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
# Create polar stereographic Basemap instance.
#
		M         = Basemap(projection='cyl', \
				llcrnrlat =START_LAT,urcrnrlat=END_LAT, \
				llcrnrlon =START_LONG,urcrnrlon=END_LONG, \
#				rsphere=6371007.181,ellps='WGS84', \
				resolution=RESOLUTION,fix_aspect=ASPECT)
		AX        = plt.gca()
#
		COLOURS     = COLOURS_ALL[iPLOT]
		CLEVELS     = CLEVELS_ALL[iPLOT]
		DATA_IMAGE  = DATA_IMAGE_ALL[iPLOT]
		PLOT_LABEL  = PLOT_LABELS[iPLOT]
		SUB_TITLE   = SUB_TITLES[iPLOT]
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
#
# Set up a colormap:
#
		if MAP_TYPE[iPLOT] == 'Map' or MAP_TYPE[iPLOT] == 'Contour0' or MAP_TYPE[iPLOT].find('JULES'):
#
			NCOLOURS  = len(CLEVELS)-1
#
			if len(COLOURS) == 0:
				COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]
#
			COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
			cm.register_cmap(cmap=COL_MAP)
			PALETTE   = COL_MAP
#
			NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
#
		elif MAP_TYPE[iPLOT] == 'Contour1':
			PALETTE  = cm.gist_rainbow
		elif MAP_TYPE[iPLOT] == 'Contour2':
			PALETTE  = cm.spectral
		elif MAP_TYPE[iPLOT] == 'Contour3':
			PALETTE  = cm.jet
#
		PALETTE.set_under(SET_UNDER[iPLOT], 1.0)
		PALETTE.set_over(SET_OVER[iPLOT],  1.0)
#
		if MAP_TYPE[iPLOT] == 'Map' or MAP_TYPE[iPLOT] == 'Map_JULES':
			IMAGE    = M.imshow(DATA_IMAGE,cmap=PALETTE,interpolation='nearest',norm=NORM)
		elif MAP_TYPE[iPLOT] == 'Contour0' or MAP_TYPE[iPLOT] == 'Contour0_JULES':
			IMAGE    = plt.contourf(LONG_MAP[iPLOT],LAT_MAP[iPLOT],DATA_IMAGE,CLEVELS,cmap=PALETTE,norm=NORM)
		else:
			IMAGE    = plt.contourf(LONG_MAP[iPLOT],LAT_MAP[iPLOT],DATA_IMAGE,CLEVELS,cmap=PALETTE)
#
# draw coastlines.
#
#		M.drawrivers(linewidth =0.1,color='grey')
		M.drawcountries(linewidth =0.2)
		M.drawcoastlines(linewidth =0.2)
#
# Draw a line around the map region.
#
		M.drawmapboundary()
#
# draw parallels and meridians.
#
		M.drawparallels(LAT, labels =[1,0,0,0], linewidth = 0.2)
		M.drawmeridians(LONG,labels =[0,0,0,1], linewidth = 0.2)
#
		AX.set_xticklabels(AX.get_xticks(),FONT)
		AX.set_yticklabels(AX.get_yticks(),FONT)
#
		DIVIDER   = make_axes_locatable(AX)
		AXIS      = DIVIDER.append_axes('bottom', size='4%', pad=0.25)
		COLORBAR  = plt.colorbar(IMAGE,cax=AXIS,extend='both',orientation='horizontal',ticks=CLEVELS)
		COLORBAR.set_label(PLOT_LABEL)
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if iDISPLAY =='Y':
#
# Display onscreen.
		plt.show() # display onscreen.
	else:
#
# Write to file
		plt.savefig(FILE_PLOT)
#
	plt.close()
#
# Return to calling routine
#
	return
#
def switch_long_axis(DATA_IN):
#
	NLONG     = DATA_IN.shape[0]
	NLONG2    = NLONG/2
	DATA_OUT  = np.zeros(NLONG)
	print NLONG,NLONG2,DATA_IN.shape,DATA_OUT.shape
#
	for iLONG in range(NLONG2):
		DATA_OUT[iLONG]          = DATA_IN[iLONG+NLONG2]
		DATA_OUT[iLONG+NLONG2]   = DATA_IN[iLONG]
#
# Return to calling routine
#
	return DATA_OUT
#
def switch_long(DATA_IN,DEBUG='N'):
#
	NLONG     = DATA_IN.shape[1]
	NLAT      = DATA_IN.shape[0]
#
	NLONG2    = NLONG/2
	DATA_OUT  = np.zeros((NLAT,NLONG))
	if DEBUG == 'Y': print NLONG,NLONG2,DATA_IN.shape,DATA_OUT.shape
#
	for iLONG in range(NLONG2):
		DATA_OUT[:,iLONG]        = DATA_IN[:,iLONG+NLONG2]
		DATA_OUT[:,iLONG+NLONG2] = DATA_IN[:,iLONG]
#
# Return to calling routine
#
	return DATA_OUT
#
def switch_long_time(DATA_IN):
#
	NLONG     = DATA_IN.shape[2]
	NLAT      = DATA_IN.shape[1]
	NTIMES    = DATA_IN.shape[0]
#
	NLONG2    = NLONG/2
	DATA_OUT  = np.zeros((NTIMES,NLAT,NLONG))
#
	for iLONG in range(NLONG2):
		DATA_OUT[:,:,iLONG]        = DATA_IN[:,:,iLONG+NLONG2]
		DATA_OUT[:,:,iLONG+NLONG2] = DATA_IN[:,:,iLONG]
#
# Return to calling routine
#
	return DATA_OUT
#
#
# *****************************************************************************************
#
# End of Program
