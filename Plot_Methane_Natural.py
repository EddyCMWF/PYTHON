#!/usr/bin/python
#
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
import numpy.ma as ma
from numpy import arange,dtype
import data_netCDF
import plot_map
import output_netCDF_Sciamachy
#
iINTER        = sys.argv[1]
DEBUG         = sys.argv[2]
SPLOT         = sys.argv[3]
SDATE         = sys.argv[4]
#
HOME_DIR      = '/prj/ALANIS/UM_Modelling/EMISSIONS/'
OUT_DIR       = HOME_DIR+'a_NATURAL/'
NETCDF_DIR    = HOME_DIR+'a_NATURAL/'
FILE_PARTS    = [ 'nvoc_land_1x1','nvoc_ocean_1x1' ]
SET_UNDER     = 'lightgrey'
MISS_DATA     = -999.9
SRESOL        =   '1.0'
RESOL         = float(SRESOL)
#
# Select output option
#
LAT_START     =  -90.0
LAT_END       =   90.0
NLAT          = int((LAT_END-LAT_START)/RESOL)
DLAT          = RESOL
DEL_LAT       =   30.0
LAT           = LAT_START + DLAT/2.0 + DLAT*np.arange(NLAT,dtype='float32')
#
LONG_START    = -180.0
LONG_END      =  180.0
NLONG         = int((LONG_END-LONG_START)/RESOL)
DLONG         = RESOL
DEL_LONG      =   30.0
LONG          = LONG_START + DLONG/2.0 + DLONG*np.arange(NLONG,dtype='float32')
#
LONG_PLOTS    = -180.0
LONG_PLOTE    =  180.0
#
NFILES        =    2
NPOINTS_ROW   =    6
START_YEAR    = 1980
END_YEAR      = 2010
NYEARS        = END_YEAR-START_YEAR+1
NMONTHS       =   12
NTIMES        = NYEARS*NMONTHS
ELAPSED_DAYS  =  -16
NSECS_DAY     = 24.0*60.0*60.0
DAYS_MONTH    = [  31 , 31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30  ]
#
NVOC          = np.zeros((NFILES,NTIMES,NLAT,NLONG))
#
# Get area of grid squares
#
FILE_CDF_A    = HOME_DIR + '/GridCell_Area_'+SRESOL+'.nc'
print(FILE_CDF_A)
DATA_NAME    = 'area'
DIMS,AREA    = data_netCDF.data_netCDF_array2D(FILE_CDF_A,DATA_NAME)
AREA = np.squeeze(AREA)
if DEBUG == 'Y':
	print(AREA.min(),AREA.max(),AREA.sum())
#
# Input data
#
for iFILE in range(NFILES):
	for iMONTH in range(NMONTHS):
		sMONTH        = '%02d' % (iMONTH+1)
		FILE_CDF_IN   = NETCDF_DIR + FILE_PARTS[iFILE] + sMONTH + '.dat'
		print(FILE_CDF_IN)
		NSECS_MONTH   = NSECS_DAY*DAYS_MONTH[iMONTH]
#
		iLONG         = -NPOINTS_ROW
		iLAT          =  0
		iLINE         =  1
#
		for LINE in open(FILE_CDF_IN):
			INPUT         = LINE.replace('\n','').split()
#
# Assign to master array - six datapoints per line
#
			if iLONG >= 360-NPOINTS_ROW:
				iLONG         = 0
				iLAT         += 1
			else:
				iLONG        += NPOINTS_ROW
#
			if DEBUG == 'Y': print iLINE,iLAT,iLONG,INPUT
			for iPOINT in range(NPOINTS_ROW):
				NVOC[iFILE,iMONTH,iLAT,iLONG+iPOINT] = float(INPUT[iPOINT])
			if DEBUG == 'Y': print NVOC[iFILE,iMONTH,iLAT,iLONG:iLONG+NPOINTS_ROW]
#
			iLINE        += 1
#
# Correct from kg per grid square per month to kg m-2 s-1
#
		NVOC[iFILE,iMONTH,:,:] = NVOC[iFILE,iMONTH,:,:]/(NSECS_MONTH*AREA)
#
# Create TIME array
#
TIME_CDF      = np.zeros(NMONTHS*NYEARS)
iTIME         =   0
#
for iYEAR in range(NYEARS):
#
	EMISS_LAND    = 0.0
	EMISS_OCEAN   = 0.0

# Check if leap year
#
	YEAR          = START_YEAR+iYEAR
	SYEAR         = str(YEAR)
#
	if YEAR ==4*int(YEAR/4):
		DAYS_MONTH[2] = 29
	else:
		DAYS_MONTH[2] = 28
	if DEBUG == 'Y': print(iYEAR,YEAR,DAYS_MONTH)
#
	for iMONTH in range(NMONTHS):
		ELAPSED_DAYS    = ELAPSED_DAYS + DAYS_MONTH[iMONTH]
		TIME_CDF[iTIME] = ELAPSED_DAYS
		CONV_FACTOR     = NSECS_DAY*DAYS_MONTH[iMONTH]/1.00E+09
#
# Replicate annual cycle after first year
#
		if iTIME >= 12:
			for iFILE in range(NFILES):
				NVOC[iFILE,iTIME,:,:] = NVOC[iFILE,iMONTH,:,:]
#
		TEMP_LAND     = NVOC[0,iTIME,:,:]*AREA
		TEMP_OCEAN    = NVOC[1,iTIME,:,:]*AREA
		EMISS_LAND    = EMISS_LAND +TEMP_LAND.sum()*CONV_FACTOR
		EMISS_OCEAN   = EMISS_OCEAN+TEMP_OCEAN.sum()*CONV_FACTOR
#
		iTIME	       += 1
#
	TEXT         = 'Global annual emissions (land, ocean) for '+SYEAR+' = '+ \
		('%12.4f %12.4f %12.4f %12.4f' % (EMISS_LAND,EMISS_LAND, \
			EMISS_OCEAN,EMISS_OCEAN))+' Tg per annum'
	print TEXT
#
if DEBUG == 'Y': print(TIME_CDF)
#
NVOC[NVOC == 0.0] = MISS_DATA
#
if SPLOT[0] == 'Y':
#
	print NVOC[iFILE,:,:,:].min(),NVOC[iFILE,:,:,:].max()
	print 'Input minimum data value: '
	DATA_MIN      = input()
#
	print 'Input maximum data value: '
	DATA_MAX      = input()
#
	print 'Input title for display for selected variable: '
	PLOT_TITLE_0  = input()
#
	print 'Input title for colour bar: '
	PLOT_LABEL    = input()
#
	for iTIME in range(NTIMES):
#
		iMONTH       = iTIME % 12
		iYEAR        = int(iTIME/12)
		if DEBUG=='Y': print(START_YEAR,iYEAR,iMONTH)
#
		SYEAR        = str(START_YEAR+iYEAR)
		SMONTH       = '%02d' % (iMONTH+1)
		FILE_PLOT    = OUT_DIR + FILE_PARTS[iFILE] + '_' + SYEAR + SMONTH + '_' + SDATE + '.png'
		PLOT_TITLE   = PLOT_TITLE_0 + ' for ' + SYEAR + SMONTH
		print('%4d %4d %4s %4s %s' % (iMONTH+1,iYEAR,SMONTH,SYEAR,FILE_PLOT))
#
		DATA_PLOT    = NVOC[1,iTIME,:,:]
#
		if LONG_END == 360.0:
			DATA_PLOT    = plot_map.switch_long(DATA_PLOT)
#
		iDISPLAY     = SPLOT[1]
		plot_map.plot_map(DATA_PLOT,LONG_PLOTS,LONG_PLOTE,DEL_LONG, \
			LAT_START,LAT_END,DEL_LAT,DATA_MAX,DATA_MIN, \
			PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,SET_UNDER,DEBUG)
#
# Call routine to output as netCDF
#
for iFILE in range(NFILES):
#
	SYEARS       = str(START_YEAR)+'_'+str(END_YEAR)
	FILE_CDF_OUT = NETCDF_DIR + FILE_PARTS[iFILE] + '_' + SYEARS + '.nc'
	VAR_INFO     = ['flux','float32']
#
	output_netCDF_Sciamachy.output_netCDF_Sciamachy(FILE_CDF_OUT,VAR_INFO, \
		NVOC[iFILE,:,:,:],TIME_CDF,LONG,LAT,NLONG,NLAT,NTIMES,START_YEAR,MISS_DATA,DEBUG)
#
# End of Program
