#!/usr/bin/python
#
# Python module to aggregate output from JULES runs
#
# Garry Hayman
# Centre for Ecology and Hydrology
# September 2013
#
# Contains
#
import os
import sys
import numpy as np
#
import data_info
import data_netCDF
import gdal_functions
import plot_map
import plot_contour
import plot_functions
#
DEBUG          = sys.argv[1]
INTER          = sys.argv[2]
MAP_TYPE       = sys.argv[3]
iDISPLAY       = sys.argv[4]
PLOT_OPT       = sys.argv[5]
#
EPSG_IN        =  4326
EPSG_OUT       = 27700
PROJ_IN,PROJ_OUT,TX=gdal_functions.define_projections(EPSG_IN,EPSG_OUT)
#
DAYS_MONTH     = [  31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31  ]

NETCDF_DIR     = '/prj/ALANIS/deposition/'
MISS_DATA      = -999.9
#
SET_UNDER      = 'white'
SET_OVER       = '#800000' # maroon
FONTSIZES      = [ 12,12,14,16 ]
RESOLUTION     = 'i'
ASPECT         = 'False'
#
COLOURS        = [ '#dcdcdc', '#9400d3', '#8a2be2', '#0000cd', '#0000ff', '#4169e1', '#6495ed', \
		   '#4682b4', '#228B22', '#00ff00', '#98fb98', '#adff2f', '#ffff00', \
		   '#ffeb00', '#ffd700', '#ffa500', '#ff8c00', '#ff0000', '#ff1493' ]
#
CHANGE_VAR     = 1
CHANGE_YEAR    = 1
CONTINUE       = 1
#
LONG_START     =    -15.0
LONG_END       =     15.0
DLONG          =      5.0
LAT_START      =     45.0
LAT_END        =     65.0
DLAT           =      5.0
#
NORTH_START    =  -4000.0
NORTH_END      =  14000.0
DNORTH         =   4000.0
EAST_START     =  -4000.0
EAST_END       =  10000.0
DEAST          =   4000.0
#
LONG_GRID      = np.arange(LONG_START, LONG_END+DLONG,  DLONG)
LAT_GRID       = np.arange(LAT_START,  LAT_END+DLAT,    DLAT )
NORTH_GRID     = np.arange(NORTH_START,NORTH_END+DNORTH,DNORTH)
EAST_GRID      = np.arange(EAST_START, EAST_END+DEAST,  DEAST )
#
# Get shapefile
#
FILE_SHAPE     = '/users/eow/garr/Work/CODE/MAP_DATA/UK_Ireland/UK&IRL_osgb.shp'
SCALE          = 100.0
PATCHES_ALL    = gdal_functions.get_shapefile(FILE_SHAPE,SCALE,DEBUG)
#
# Input data
#
while CONTINUE == 1:
#
	if CHANGE_YEAR == 1:
		print '\nInput date (YYYYMMDD): '
		SDATE          = input()
		SYEAR          = SDATE[0:4]
		YEAR           = int(SYEAR)
		MONTH          = int(SDATE[4:6])
		DAY            = int(SDATE[6:8])
#
	FILE_NAME      = 'EMEP4UK_Base_emep_4.3_'+SYEAR+'_day'
	FILE_NAME_CDF  = NETCDF_DIR+FILE_NAME+'.nc'
#
# (1) Get variable
#
	if CHANGE_VAR == 1:
		VAR_NAMES      = data_netCDF.data_netCDF_getVARNAMES(FILE_NAME_CDF) 
		print '\nSelect variable'
		for iVAR in range(len(VAR_NAMES)):
			TEXT         = ('%04d %s' % (iVAR+1,VAR_NAMES[iVAR]))
			print(TEXT)
#
# (2) Select variable
#
		print '\nInput index of variable to plot: '
		iVAR           = int(input())-1
		VAR_NAME       = VAR_NAMES[iVAR]
#
# Get data
#
		DIMS,LATS      = data_netCDF.data_netCDF_array_var(FILE_NAME_CDF,'lat')
		DIMS,LONGS     = data_netCDF.data_netCDF_array_var(FILE_NAME_CDF,'lon')
		DIMS,DATA      = data_netCDF.data_netCDF_array_var(FILE_NAME_CDF,VAR_NAME)
#
# Convert LATS and LONGS to new projection
#
		OS_EASTS,OS_NORTHS=gdal_functions.convert_points(PROJ_IN,PROJ_OUT,TX,LONGS,LATS,DEBUG)
#
	VAR_NAME_PLT = VAR_NAME.replace('_',' ')
	PLOT_TITLE   = VAR_NAME_PLT+' '+SDATE
	PLOT_LABEL   = VAR_NAME_PLT
#
	iJULIAN      = DAY
	for iMONTH in range(1,MONTH):
		iJULIAN      = iJULIAN+DAYS_MONTH[iMONTH-1]
	print YEAR,MONTH,DAY,iJULIAN
#
	print '\nMinimum and maximum values: '
	print DATA[iJULIAN,:,:].min(),DATA[iJULIAN,:,:].max()
	print 
	print 'Input minimum, maximum and increment: '
	INPUT          = input().split(':')
	MAP_MIN        = float(INPUT[0])
	MAP_MAX        = float(INPUT[1])
	MAP_INC        = float(INPUT[2])
#
	NLEVELS        = int((MAP_MAX-MAP_MIN)/MAP_INC)+1
	CLEVELS        = MAP_MIN+MAP_INC*np.arange(NLEVELS+1)
#

	if PLOT_OPT[0] == 'Y':
		WIDTH          = 12.0
		HEIGHT         =  8.0
		FILE_PLOT      = '/users/eow/edwcom/test.png'  #NETCDF_DIR+'/PLOTS/'+FILE_NAME+'_'+VAR_NAME+'.png'
		print FILE_PLOT
#
		plot_map.plot_map3(DATA[iJULIAN-1,:,:],LONG_START,LONG_END,LONG_GRID,LONGS, \
			LAT_START,LAT_END,LAT_GRID,LATS,CLEVELS,COLOURS,MAP_TYPE,WIDTH,HEIGHT,ASPECT, \
			RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES,SET_UNDER,SET_OVER,DEBUG) 
#
	if PLOT_OPT[1] == 'Y':
		FILE_PLOT    = NETCDF_DIR+'/PLOTS/'+FILE_NAME+'_'+VAR_NAME+'_regrid.png'
		print FILE_PLOT
#
		NROWS          = 1
		NCOLUMNS       = 2
		NPLOTS         = NROWS*NCOLUMNS
		NCONTOURS      = 1
		WIDTH          = 12.0
		HEIGHT         = 12.0
		LAYOUT_DATA    = [[4,8],[1,4,0,4],[0,4,5,8]]
#
		NLEVELS_ALL    = []
		CLEVELS_ALL    = []
		COLOURS_ALL    = []
		MAP_TYPE_ALL   = []
		PLOT_LABELS    = []
		SET_UNDERS     = []
		SET_OVERS      = []
#
		DATA_ALL       = [DATA[iJULIAN-1,:,:],DATA[iJULIAN-1,:,:]]
		SUB_TITLES     = ['Original','Regridded']
#
		START_X        = [LONG_START, EAST_START   ]
		END_X          = [LONG_END,   EAST_END     ]
		X_ALL          = [LONGS,      OS_EASTS/100.0 ]
		XTICKS         = [LONG_GRID,  EAST_GRID    ]
		XTICKS_LABELS  = [LONG_GRID,  EAST_GRID    ]
		XLABELS        = ['',         'OS Easting' ]
#
		START_Y        = [LAT_START,  NORTH_START  ]
		END_Y          = [LAT_END,    NORTH_END    ]
		Y_ALL          = [LATS,       OS_NORTHS/100.0]
		YTICKS         = [LAT_GRID,   NORTH_GRID   ]
		YTICKS_LABELS  = [LAT_GRID,   NORTH_GRID   ]
		YLABELS        = ['',         'OS Northing']
#
		PLOT_LABELS    = [PLOT_LABEL,PLOT_LABEL]
		MAP_DATA       = ['i','N']
		PATCH_DATA     = [1,PATCHES_ALL]
#
		for iPLOT in range(NPLOTS):
			NLEVELS_ALL.append(NLEVELS)
			CLEVELS_ALL.append(CLEVELS)
			COLOURS_ALL.append(COLOURS)
			MAP_TYPE_ALL.append(MAP_TYPE)
			SET_UNDERS.append(SET_UNDER)
			SET_OVERS.append(SET_OVER)
#
		plot_contour.plot_contour_EMEP_multi(NPLOTS,NROWS,NCOLUMNS,NCONTOURS,DATA_ALL, \
			START_X,END_X,X_ALL,XTICKS,XTICKS_LABELS,XLABELS, \
			START_Y,END_Y,Y_ALL,YTICKS,YTICKS_LABELS,YLABELS, \
			CLEVELS_ALL,COLOURS_ALL,MAP_TYPE_ALL,WIDTH,HEIGHT, \
			PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,iDISPLAY, \
			FONTSIZES,SET_UNDERS,SET_OVERS,DEBUG,MAP_DATA,LAYOUT_DATA,PATCH_DATA)
#
	print 'Change variable (1/0): '
	CHANGE_VAR    = int(input())
#
	print 'Continue (1/0): '
	CONTINUE      = int(input())
#
# End of Program
