# Python module to selected datasets
#
# Garry Hayman
# Centre for Ecology and Hydrology
# December 2011
#
# Contains
#
import sys
import numpy as np
from numpy import arange,dtype
import data_netCDF
import plot_map
#
iINTER        = sys.argv[1]
DEBUG         = sys.argv[2]
SDATE         = sys.argv[3]
START_YEAR    = int(sys.argv[4])
RESOL         = float(sys.argv[5])
#
OUT_DIR       = '/prj/ALANIS/z_TEST_PLOTS/'
SET_UNDER     = 'lightgrey'
MISS_DATA     = -999.9
#
# Select output option
#
LAT_START     = -90.0
LAT_END       =  90.0
NLAT          = (LAT_END-LAT_START)/RESOL
DLAT          = RESOL
DEL_LAT       = 30.0
#
LONG_START    =   0.0
LONG_END      = 360.0
NLONG         = (LONG_END-LONG_START)/RESOL
DLONG         = RESOL
DEL_LONG      = 30.0
#
if LONG_END == 360.0:
	LONG_PLOTS   = -180.0
	LONG_PLOTE   =  180.0
#
# Input data
#
print('Input name of directory containing netCDF file : ')
NETCDF_DIR    = input()
#
print('Input name of netCDF file : ')
FILE_IN       = input()
#
print('Input name of variable to plot: ')
DATA_NAME     = input()
#
print('Input name of time variable: ')
TIME_NAME     = input()
#
FILE_CDF_IN  = NETCDF_DIR + '/' + FILE_IN + '.nc'
print(FILE_CDF_IN)
#
DIMS,VAR_DATA_IN,TIME_IN = data_netCDF.data_netCDF_array(FILE_CDF_IN,DATA_NAME,TIME_NAME)
#
print('Minimum and maximum data values are :',VAR_DATA_IN.min(),VAR_DATA_IN.max(),'\n')
print('Input minimum data value: ')
DATA_MIN      = input()
#
print('Input maximum data value: ')
DATA_MAX      = input()
#
print('Input scale factor: ')
SCALE_FAC     = input()
#
print('Input title for display for selected variable: ')
PLOT_TITLE_0  = input()
#
print('Input title for colour bar: ')
PLOT_LABEL    = input()
PLOT_LABEL    = PLOT_LABEL + ' x ' + str(SCALE_FAC)
#
for iTIME in range(VAR_DATA_IN.shape[0]):
#
	iMONTH       = iTIME % 12
	iYEAR        = int(iTIME/12)
	if DEBUG=='Y': print(START_YEAR,iYEAR,iMONTH)
#
	SYEAR        = str(START_YEAR+iYEAR)
	SMONTH       = '%02d' % (iMONTH+1)
	FILE_PLOT    = OUT_DIR + FILE_IN + '_' + SYEAR + SMONTH + '_' + SDATE + '.png'
	PLOT_TITLE   = PLOT_TITLE_0 + ' for ' + SYEAR + SMONTH
	print('%4d %4d %4s %4s %s' % (iMONTH+1,iYEAR,SMONTH,SYEAR,FILE_PLOT))
#
	DATA_PLOT    = VAR_DATA_IN[iTIME,:,:]
	DATA_PLOT    = DATA_PLOT*SCALE_FAC
	if LONG_END == 360.0:
		DATA_PLOT    = plot_map.switch_long(DATA_PLOT)
#
	iDISPLAY     = 'Y'
	plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
		LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
		PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
	iDISPLAY     = 'N'
	plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
		LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
		PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
# End of Program
