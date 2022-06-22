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
import calculate_area
import data_netCDF
import output_netCDF_Sciamachy
import output_netCDF_UM
import plot_map
#
iINTER        = sys.argv[1]
DEBUG         = sys.argv[2]
sRESOL        = sys.argv[3]
SDATE         = sys.argv[4]
START_YEAR    = int(sys.argv[5])
END_YEAR      = int(sys.argv[6])
iDISPLAY      = sys.argv[7]
DATA_MIN      = float(sys.argv[8])
DATA_MAX      = float(sys.argv[9])
NO_SECS_MONTH = 30.0*24.0*60.0*60.0
#
# Assign parameters
#
NYEARS        = END_YEAR-START_YEAR+1
NMONTHS       = 12
NTIMES        = NMONTHS
NETCDF_DIR    = '/prj/ALANIS/UM_Modelling/EMISSIONS/'
OUT_DIR       = '/prj/ALANIS/z_TEST_PLOTS/'
TIME_NAME     = 't'
SET_UNDER     = 'lightgrey'
#
# Select output option
#
PLOT_TITLE_0   = 'CH4 Emissions '
PLOT_LABEL     = 'CH4 Emissions (in Mtonnes CH4 per annum)'
MISS_DATA      = -999.9
DAYS_MONTH     = [  31 , 31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30  ]
#
LAT_START      = -90.0
LAT_END        =  90.0
DEL_LAT        =  30.0
DLAT           = LAT_END-LAT_START
#
LONG_START     =   0.0
LONG_END       = 360.0
DEL_LONG       =  30.0
DLONG          = LONG_END-LONG_START
#
if sRESOL=='UM':
	RESOL_LONG     = 1.875
	RESOL_LAT      = 1.250
	NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)+1
	NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
else:
	RESOL_LONG     = float(sRESOL)
	RESOL_LAT      = float(sRESOL)
	NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
	NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
#
# - Create longitude and latitude values for output (Note these are at centre of grid squares)
#
LONG           = LONG_START + 0.5 + DLONG*arange(NLONG,dtype='float32')
LAT            = LAT_START  + 0.5 + DLAT*arange(NLAT,dtype='float32')
#
if LONG_END == 360.0:
	LONG_PLOTS   = -180.0
	LONG_PLOTE   =  180.0
#
# Calculate area of grid squares
#
NLONG_A,NLAT_A,LONG_A,LAT_A,AREA= \
	calculate_area.calculate_area_var_UM(LONG_START,DLONG,RESOL_LONG,LAT_START,DLAT,RESOL_LAT,'N')
if DEBUG == 'Y': print(NLONG_A,NLAT_A,LONG_A,LAT_A)
#
# Input data
#
# (a) All CH4 sources
#
FILE_CDF_IN  = NETCDF_DIR + 'StdTrop_AR5_surfems_2000.nc'
print(FILE_CDF_IN)
#
DATA_NAME    = 'longitude'
DIMS,LONG,TIME_IN = data_netCDF.data_netCDF_array(FILE_CDF_IN,DATA_NAME,TIME_NAME)
if DEBUG == 'Y': print(LONG.min(),LONG.max())
#
DATA_NAME    = 'latitude'
DIMS,LAT,TIME_IN = data_netCDF.data_netCDF_array(FILE_CDF_IN,DATA_NAME,TIME_NAME)
if DEBUG == 'Y': print(LAT.min(),LAT.max())
#
DATA_NAME    = 'ch4_surf_emiss'
DIMS,CH4_EMISS_ALL,TIME_IN = data_netCDF.data_netCDF_array(FILE_CDF_IN,DATA_NAME,TIME_NAME)
if DEBUG == 'Y': print(CH4_EMISS_ALL.shape,CH4_EMISS_ALL.min(),CH4_EMISS_ALL.max())
#
# (b) Without CH4 wetland emissions
#
FILE_CDF_IN  = NETCDF_DIR + 'StdTrop_AR5_surfems_nowetl_2000.nc'
print(FILE_CDF_IN)
#
DATA_NAME    = 'ch4_surf_emiss'
DIMS,CH4_EMISS_NOWET,TIME_IN = data_netCDF.data_netCDF_array(FILE_CDF_IN,DATA_NAME,TIME_NAME)
if DEBUG == 'Y': print(CH4_EMISS_NOWET.shape,CH4_EMISS_NOWET.min(),CH4_EMISS_NOWET.max())
#
CH4_EMISS_WET= CH4_EMISS_ALL-CH4_EMISS_NOWET
#
for iTIME in range(CH4_EMISS_ALL.shape[0]):
#
	iMONTH       = iTIME % 12
	iYEAR        = int(iTIME/12)
	if DEBUG=='Y': print(START_YEAR,iYEAR,iMONTH)
#
	SYEAR        = str(START_YEAR+iYEAR)
	SMONTH       = '%02d' % (iMONTH+1)
#
	FILE_PLOT    = OUT_DIR + 'CH4_Emissions_All_' + SYEAR + SMONTH + '_' + SDATE + '.png'
	PLOT_TITLE   = PLOT_TITLE_0 + ' (all sources) for ' + SYEAR + SMONTH
	print('%4d %s' % (iTIME+1,FILE_PLOT))
#
	DATA_PLOT    = CH4_EMISS_ALL[iTIME,0,:,:]*AREA[:,:]*NO_SECS_MONTH/1.0E+09
	print(DATA_PLOT.shape,DATA_PLOT.min(),DATA_PLOT.max())
	if LONG_END == 360.0:
		DATA_PLOT    = plot_map.switch_long(DATA_PLOT)
#
	plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
		LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
		PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
	FILE_PLOT    = OUT_DIR + 'CH4_Emissions_nowet_' + SYEAR + SMONTH + '_' + SDATE + '.png'
	PLOT_TITLE   = PLOT_TITLE_0 + ' (all sources except wetlands) for ' + SYEAR + SMONTH
	print('%4d %s' % (iTIME+1,FILE_PLOT))
#
	DATA_PLOT    = CH4_EMISS_NOWET[iTIME,0,:,:]*AREA[:,:]*NO_SECS_MONTH/1.0E+09
	print(DATA_PLOT.shape,DATA_PLOT.min(),DATA_PLOT.max())
	if LONG_END == 360.0:
		DATA_PLOT    = plot_map.switch_long(DATA_PLOT)
#
	plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
		LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
		PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
	FILE_PLOT    = OUT_DIR + 'CH4_Emissions_wet_' + SYEAR + SMONTH + '_' + SDATE + '.png'
	PLOT_TITLE   = PLOT_TITLE_0 + ' (wetlands) for ' + SYEAR + SMONTH
	print('%4d %s' % (iTIME+1,FILE_PLOT))
#
	DATA_PLOT    = CH4_EMISS_WET[iTIME,0,:,:]*AREA[:,:]*NO_SECS_MONTH/1.0E+09
	print(DATA_PLOT.shape,DATA_PLOT.min(),DATA_PLOT.max())
	if LONG_END == 360.0:
		DATA_PLOT    = plot_map.switch_long(DATA_PLOT)
#
	plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
		LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
		PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
# Output as netCDF
#
# - Call routine to output Column CH4 as netCDF
#
#	FILE_CDF_OUT = NETCDF_DIR + JOB_ID + '_' + SYEAR + '_column_ch4_' + SDATE + '.nc'
#	VAR_INFO     = ['XCH4','float32']
#
#	output_netCDF_Sciamachy.output_netCDF_Sciamachy(FILE_CDF_OUT,VAR_INFO,VAR_DATA_OUT, \
#        	TIME_CDF,LONG,LAT,NLONG,NLAT,NTIMES,START_YEAR,MISS_DATA,DEBUG)
#
# End of Program
