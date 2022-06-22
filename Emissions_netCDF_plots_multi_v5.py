#!/usr/bin/python
#
# Python module to derive emissions from Biomass burning
#
# Garry Hayman
# Centre for Ecology and Hydrology
# July 2013
#
# python Emissions_netCDF_plots_multi.py Y Y 2000 2000 NNNNNNNNNNNYNNN 1 Contour0 1 ~/Work/CODE/GOOGLE_EARTH/ YY
# python Emissions_netCDF_plots_multi.py Y N 2000 2000 NNNNNNNNNNNYNNN 1 Map 1 ~/Work/CODE/GOOGLE_EARTH/ YY
#
# ./Emissions_netCDF_plots_multi.py N N 1997 2009 NNNNNNNNNNNNNNN 1 Map 1 ~/Work/CODE/GOOGLE_EARTH/ NNN Y Y 0 4 > /prj/wetlands_africa/Sciamachy/a_PLOTS/Wetland_Emissions_Sudd_JULES_20140327.out
#
# ./Emissions_netCDF_plots_multi.py N N 1997 2009 NNNNNNNNNNNNNNN 1 Map 1 ~/Work/CODE/GOOGLE_EARTH/ NNN Y Y 0 5 > /prj/wetlands_africa/Sciamachy/a_PLOTS/Wetland_Emissions_Sudd_JULES_GIEMS_20140327.out
#
# PLOTS[0]:  Annual/monthly emission maps
# PLOTS[1]:  Emission time series
# PLOTS[2]:  Emission annual cycle
# PLOTS[3]:  Latitudinal  (zonal) plot (single/multi)
# PLOTS[4]:  Longitudinal (meridional) plot (single/multi)
# PLOTS[5]:  Use wetland fraction
# PLOTS[6]:  Emission time series (multi)
# PLOTS[7]:  Emission annual cycle (multi)
# PLOTS[9]:  Maps of wetland fraction
# PLOTS[10]: Emission anomaly time series (multi)
# PLOTS[11]: Images and kmz files for Google Earth
# PLOTS[12]: Maps of wetland fraction and emissions
# PLOTS[13]: CarbonSAT plots
# PLOTS[14]: Site-specific plots
#
# Contains
#
import os
import sys
import numpy as np
from numpy import array,arange,dtype
#
import data_info
import data_netCDF
import data_regrid_new
import merge_Data
import plot_map
import plot_functions
import Emissions_Wetlands_Zonal
#
iINTER         = sys.argv[1]
DEBUG          = sys.argv[2]
START_YEAR     = int(sys.argv[3])
END_YEAR       = int(sys.argv[4])
PLOTS          = sys.argv[5]
PLOT_OPT       = sys.argv[6]
MAP_TYPE       = sys.argv[7]
iMASS          = sys.argv[8]
KMZ_DIR        = sys.argv[9]
KMZ_CODES      = sys.argv[10]
iFIG_PAPER     = sys.argv[11]
iUSE_STORED    = sys.argv[12]
iLEAP          = sys.argv[13]
#
NETCDF_DIR     = '/prj/ALANIS/UM_Modelling/EMISSIONS/'
TIME_NAME      = 'time'
#
DAYS_MONTH     = [  31 , 31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31  ]
SMONTHS        = ['January','February','March','April','May','June', \
		  'July','August','September','October','November','December' ]
#
ORIENTATION    = 'Portrait'
FONTSIZES      = [ 12,12,14,16 ]
NROWS_0        = 4
NCOLUMNS_0     = 3
#
PLOT_TEXT      = ''        # Not used
SET_UNDER      = '#ffffff' # white
SET_OVER       = '#800000' # maroon
SET_UNDER_DIFF = '#191970'
PLOT_CODES_GR  = ['ko-','ro-','bo-','go-','mo-','yo-','co-' ]
PLOT_CODES_AN  = ['k-' ,'r-' ,'b-' ,'g-' ,'m-' ,'y-' ,'c-' ]
PLOT_CODES_ZL  = ['r:' ,'r-' ,'b-' ,'g-' ,'k-' ,'m-' ,'y-'  ,'c-' ]
MAX_TREND_0    = 100
MAX_CLIM_0     =  50
MAX_SITE_0     = 400
MAX_ANOM_0     =  50
MIN_EMISS      = 1.00E-20
PLOT_SCALE     = [1.00,0.50,0.20,0.10,0.05,0.02,0.01]
PLOT_MAX       = [0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0,5000.0]
DATA_NAMES     = ['ch4_surf_emiss','fch4wetl','fch4wetl_wet','FCH4WETL','flux', \
	'agriwaste','fossil','bbur','wetlands','others','soils','total','Prior']
JULES_DATES    = ['20120203','20120224','20120405','20120522','20121003','20130319']
LEGENDS_MULTI  = []
CLEVELS_C      = [  0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 ]
COLOURS_C      = [ '#1e3cff', '#00a0ff', '#00c8c8', '#00d28c', '#00dc00', '#a0e632', '#e6dc32', '#e6af2d', \
		   '#f08228', '#fa3c3c', '#f00082' ]
#
COLOURS_D      = [ '#1e3cff', '#4169e1', '#00bfff', '#87cefa', '#ffc0cb', '#fa8072', '#fa3c3c', '#ff0000' ]
CLEVELS_D      = [ -1.00,-0.75,-0.50,-0.25,0.0,0.25,0.50,0.75,1.00 ]
#
SITE_DATA      = [ \
#		    [ 'Salmisuo'  ,  '31.00',  '63.00'], \
		    [ 'Salmisuo'  ,  '30.93',  '62.75'], \
#		    [ 'Degero'    ,  '19.55',  '64.18'], \
		    [ 'Degero'    ,  '19.55',  '64.30'], \
#		    [ 'BOREAS-NSA', '-99.00',  '56.00'], \
		    [ 'BOREAS-NSA', '-99.00',  '56.00'], \
#		    [ 'Minnesota' , '-93.47',  '47.53'], \
		    [ 'Minnesota' , '-93.47',  '47.53'], \
#		    [ 'Michigan'  , '-84.02',  '42.45'], \
		    [ 'Michigan'  , '-84.02',  '42.45'], \
#		    [ 'Panama'    , '-79.63',   '9.31']  \
		    [ 'Panama'    , '-80.00',   '9.00']  \
		 ]
#
if PLOT_OPT == '1':
	PLOT_MAP       = 'Y'
else:
	PLOT_MAP       = 'N'
#
NSITES         = len(SITE_DATA)
NYEARS         = END_YEAR-START_YEAR+1
NMONTHS        = 12
NTIMES_1       =  1
NTIMES         = NMONTHS*NYEARS
NTRANS         = 12
iNAME          =  0
#
LAT_ZONAL_ALL  = []
LON_MERID_ALL  = []
DATA_MONTH_ALL = []
DATA_MONTH_ALL2= []
DATA_ANN_ALL   = []
EMISS_FACTORS  = []
SOURCE_ALL     = []
SOURCE_FILES   = []
SOURCE_KMZ     = []
SRESOL_KMZ     = []
COLOURS_KMZ    = []
CLEVELS_KMZ    = []
EMISS_KMZ      = []
NAMES_KMZ      = []
NDIMS_KMZ      = []
LATS_KMZ       = []
LONGS_KMZ      = []
LATS_ALL       = []
LONGS_ALL      = []
RESOL_LAT_ALL  = []
RESOL_LONG_ALL = []
EMISS_SITE     = []
AREAS_ALL      = []
UNITS_KMZ      = []
SDATE_TIMES_KMZ= []
EDATE_TIMES_KMZ= []
WET_FRACT_KMZ  = []
MISS_DATA_KMZ  = []
MON_EMISS_ALL  = []
NTIMES_ALL     = []
#
if PLOTS[9] == 'Y':
	print('\nMap wetland fraction (Y/N) : ')
	MAP_FRACT      = input()
else:
	MAP_FRACT      = 'N'
#
DOMAIN         = data_info.data_DOMAIN_info()
#
for i in range(len(DOMAIN)):
	print DOMAIN[i][0]
iDOMAIN        = int(input())
#
SDOMAIN        = DOMAIN[iDOMAIN][1]
LONG_DOMS      = float(DOMAIN[iDOMAIN][2])
LONG_DOME      = float(DOMAIN[iDOMAIN][3])
LAT_DOMS       = float(DOMAIN[iDOMAIN][5])
LAT_DOME       = float(DOMAIN[iDOMAIN][6])
LAT_START      = -90.0
LAT_END        =  90.0
DEL_LAT        = int(DOMAIN[iDOMAIN][7])
DEL_LONG       = int(DOMAIN[iDOMAIN][4])
DLAT           = LAT_END-LAT_START
WIDTH0         = float(DOMAIN[iDOMAIN][8])
HEIGHT0        = float(DOMAIN[iDOMAIN][9])
RESOLUTION     = DOMAIN[iDOMAIN][10]
ASPECT         = DOMAIN[iDOMAIN][11]
#
SRESOL_OUT     = '5.0'
O_RESOL_LONG   = float(SRESOL_OUT)
O_RESOL_LAT    = float(SRESOL_OUT)
NLONG_OUT      = int(360.0/O_RESOL_LONG)
NLAT_OUT       = int(180.0/O_RESOL_LAT)
NLONG_AREA     = int((LONG_DOME-LONG_DOMS)/O_RESOL_LONG)
NLAT_AREA      = int((LAT_DOME-LAT_DOMS)/O_RESOL_LAT)
#
if NYEARS > 14:
	XINC           =  5
elif NYEARS > 6:
	XINC           =  2
elif NYEARS < 6:
	XINC           =  1
#
if iUSE_STORED == 'Y':
	DATA_SOURCES,DATA_STORED = data_info.data_DATA_SOURCES(iUSE_STORED)
else:
	DATA_SOURCES             = data_info.data_DATA_SOURCES(iUSE_STORED)
#
FIG_PLOT_SCALES= data_info.data_PLOT_SCALES_info()
#
print('Input number of files: ')
NFILES         = int(input())
#
XDATA_MULTI    = np.zeros((NFILES,NTRANS,NYEARS))
YDATA_MULTI    = np.zeros((NFILES,NTRANS,NYEARS))
XCLIM_MULTI    = np.zeros((NFILES,NTRANS,NMONTHS))
YCLIM_MULTI    = np.zeros((NFILES,NTRANS,NMONTHS))
EMISS_ANNUAL   = np.zeros((NYEARS,NFILES))
EMISS_DOMAIN   = np.zeros((NYEARS,NFILES))
EMISS_DOMAIN_M = np.zeros((NMONTHS*NYEARS,NFILES))
TRANS_EMISS    = np.zeros((NYEARS,NTRANS))
TRANS_EMISS_CL = np.zeros((NMONTHS,NTRANS))
TRANS_EMISS_VL = np.zeros((NMONTHS,NTRANS))
YSCALE         = np.ones((NTRANS))
#
TRANS_EMISS_MN = np.zeros((NFILES,NTRANS,NMONTHS*NYEARS))
TRANS_EMISS_RM = np.zeros((NFILES,NTRANS,NMONTHS*NYEARS))
TRANS_EMISS_AN = np.zeros((NFILES,NTRANS,NMONTHS*NYEARS))
#
for iFILE in range(NFILES):
#
	TRANS_EMISS[:,:]    = 0.0
	TRANS_EMISS_CL[:,:] = 0.0
	TRANS_EMISS_VL[:,:] = 0.0
#
	SDATE_TIMES_KMZ.append([])
	EDATE_TIMES_KMZ.append([])
#
	print
	print 'Data Source'
	for i in range(len(DATA_SOURCES)):
		print DATA_SOURCES[i][0]
	TEXT           = 'Input dependent on value of iUSE_STORED (X or X.YY) for file %3d:' % (iFILE+1)
	print TEXT
#
	if iUSE_STORED == 'N':
		iSOURCE        = int(input())
	else:
		dSOURCE        = str(input())
		iSOURCE        = int((dSOURCE.split('.'))[0])
#
	DATA_SOURCE    = DATA_SOURCES[iSOURCE][1]
	DATA_LABEL     = DATA_SOURCES[iSOURCE][2]
	SRESOL         = DATA_SOURCES[iSOURCE][3]
	START_DATA     = int(DATA_SOURCES[iSOURCE][4])
	END_DATA       = int(DATA_SOURCES[iSOURCE][5])
	MISS_DATA      = float(DATA_SOURCES[iSOURCE][6])
	LONG_END       = DATA_SOURCES[iSOURCE][7]
	GRIDDED        = DATA_SOURCES[iSOURCE][8]
#
	if iUSE_STORED == 'N':
		SOURCE_FILES.append(DATA_LABEL)
		print 'Provide run details (e.g., RCP scenario): '
		DETAILS        = str(input())
	else:
		SOURCE_FILES.append(DATA_STORED[dSOURCE][0])
		DETAILS        = DATA_STORED[dSOURCE][0]
#
	print DATA_LABEL
	if 'XXXX' in DATA_LABEL:
		DATA_LABEL     = DETAILS
	print DETAILS
	print DATA_LABEL
#
	SOURCE_ALL.append(DATA_LABEL)
	SOURCE_KMZ.append(DATA_SOURCE)
	SRESOL_KMZ.append(SRESOL)
	MISS_DATA_KMZ.append(MISS_DATA)
#
	if LONG_END == 360.0:
		LONG_START =    0.0
		LONG_PLOTS = -180.0
		LONG_PLOTE =  180.0
	else:
		LONG_START = -180.0
		LONG_PLOTS = -180.0
		LONG_PLOTE =  180.0
#
	DLONG          = LONG_END-LONG_START
#
# Get area of grid squares
#
	FILE_CDF_A   = NETCDF_DIR + '/GridCell_Area_'+SRESOL+'.nc'
	print(FILE_CDF_A)
	DATA_NAME    = 'area'
	DIMS,AREA    = data_netCDF.data_netCDF_array2D(FILE_CDF_A,DATA_NAME)
	AREA         = np.squeeze(AREA)
	if DEBUG == 'Y':
		print(AREA.min(),AREA.max(),AREA.sum())
#
# Get TRANSCOM regions
#
	DATA_NAME    = 'transcom_regions'
	FILE_TRANS   = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_'+SRESOL+'.nc'
	print FILE_TRANS
	DIMS,TRANSCOM= data_netCDF.data_netCDF_array2D(FILE_TRANS,DATA_NAME)
	TRANSCOM     = np.squeeze(TRANSCOM)
#
# Switch E-W hemisphere
#
	TRANSCOM    = plot_map.switch_long(TRANSCOM)
#
# Correct AREA and TRANSCOM if SREOLS is UM or MACC
#
	if SRESOL == 'UM' or SRESOL == 'MACC':
		AREA         = AREA[:-1,:]
		TRANSCOM     = TRANSCOM[:-1,:]
#
	AREAS_ALL.append(AREA)
#
	TRANS_REGS   = data_info.data_TRANSCOM_info()
#
# Start for wetland fraction
#
	wTIME        = (START_YEAR-START_DATA)*12
#
	if 'JULES_Wetlands' in DATA_SOURCE:
#
		START_DATA     = START_YEAR
		END_DATA       = END_YEAR
#
		MODIFIER      = [ \
				['0: original',''], \
				['1: masked by EO inundation','+Masked'], \
				['2: driven with EO inundation','+GIEMS'], \
				]
#
		print("\nJULES Wetland option")
		for i in range(len(MODIFIER)):
			print(MODIFIER[i][0])
#
		JULES_OPT,OPT_CORR=data_info.data_WETLAND_EMISS_info()
		FILE_CORR     = ''
		iOPT_CORR     = -1
#
		if iUSE_STORED == 'N':
#
# Input from keyboard
#
			print('\nInput option      : ')
			jFILE         = input()
#
			print("\nJULES Run Dates")
			for i in range(len(JULES_DATES)):
				print(i,JULES_DATES[i])
			print('\nInput index for date of JULES run (0-'+str(len(JULES_DATES)-1)+'): ')
			iDATA         = int(input())
#
			print('\nInput JULES run (e.g., m43, m45): ')
			JULES_RUN     = input()
#
			print('\nMean annual wetland emission (Tg per annum): ')
			WET_TOTAL     = float(input())

			if iDATA < 4 and JULES_RUN != 'm43':
				print 'STOP - incompatible Jules run and data option'
				quit()
#
			if iINTER == 'Y':
				print("\nJULES Parameterisation")
				for i in range(len(JULES_OPT)):
					print(JULES_OPT[i][0])
#
				print('\nInput version     : ')
				iNPP	  = input()
			else:
				iNPP	  = int(sys.argv[14])
#
			if iDATA >= 4:
#
				print("Correction option")
				for i in range(len(OPT_CORR)):
					print(OPT_CORR[i][0])
#
				if iINTER == 'Y':
					print('\nInput correction option: ')
					iOPT_CORR     = input()
				else:
					iOPT_CORR     = int(sys.argv[15])
#
				FILE_CORR     = OPT_CORR[iOPT_CORR][1]
#
		else:
#
# Get parameters from stored values
#
			jFILE         = int(DATA_STORED[dSOURCE][9])
			iDATA         = int(DATA_STORED[dSOURCE][10])
			JULES_RUN     = DATA_STORED[dSOURCE][11]
			WET_TOTAL     = float(DATA_STORED[dSOURCE][12])
#
			if iDATA < 4 and JULES_RUN != 'm43':
				print 'STOP - incompatible Jules run and data option'
				quit()
#
			iNPP	  = int(DATA_STORED[dSOURCE][13])
#
			if iDATA >= 4:
				iOPT_CORR     = int(DATA_STORED[dSOURCE][14])
				FILE_CORR     = OPT_CORR[iOPT_CORR][1]
#
# Correct iDATA for iNPP option
#
		FILE_PART     = JULES_OPT[iNPP][1]
#
# Filenames
#
		TEXT_NAMES	      = [ 'CH4 Wetlands' ]
#
		FILE_DATA	       = [ \
					    [ 'M','CH4','fch4wetl', \
					      '/Wetland_Emissions_CH4_'+FILE_PART+'_'+JULES_RUN, \
					      '20120203','20120224','20120405','20120522','20121003','20130319','20130402'], \
					    [ 'M','CH4','fch4wetl', \
					      '/Wetland_Emissions_CH4_'+FILE_PART+'_'+JULES_RUN+'_Masked', \
					      '20120203','20120224','20120405','20120522','20121003','20130319','20130402'], \
					    [ 'M','CH4','fch4wetl', \
					      '/Wetland_Emissions_CH4_'+FILE_PART+'_'+JULES_RUN+'_EO', \
					      '20120203','20120224','20120405','20120522','20121003','20130319','20130402'] \
					  ]
#
# Get data - wetland fraction
#
		if PLOTS[5] == 'Y' or PLOTS[11] == 'Y' or PLOTS[14] == 'Y' or DATA_SOURCE == 'JULES_Wetlands_Area':
#
			SET_OVER     = 'k' # black
			SDATE        = FILE_DATA[jFILE][4+iDATA]
			FILENAME     = NETCDF_DIR+'a_WETLANDS_'+SDATE+FILE_DATA[jFILE][3]+'_1993_2009'+FILE_CORR+'.nc'
			FILENAME     = FILENAME.replace('Wetland_Emissions_CH4','Wetland_Fraction')
			print(' ')
			print(FILENAME)
#
			DATA_NAME    = 'fwetl'
			DIMS,WET_FRACT = data_netCDF.data_netCDF_array_var(FILENAME,DATA_NAME)
			WET_FRACT    = np.squeeze(WET_FRACT)
#
		if DATA_SOURCE == 'JULES_Wetlands_Area':
			EMISS       = WET_FRACT.copy()
		else:
#
			EMISS       = []
			iTIME       = 0
			iREAD       = 0
			sCODE       = FILE_DATA[jFILE][0]
			FSPECIES    = FILE_DATA[jFILE][1]
			SDATE       = FILE_DATA[jFILE][4+iDATA]
#
# Get data - wetland emissions
#
			for iYEAR in range(NYEARS):
#
				YEAR        = START_YEAR+iYEAR
				SYEAR       = '%4d' % (YEAR)
				FILENAME    = NETCDF_DIR+'a_WETLANDS_'+SDATE+FILE_DATA[jFILE][3]+'_'+SYEAR+FILE_CORR+'.nc'
				FILE_PLOT_PART= FILENAME[:-3].replace('.','_')
				print(' ')
#
				DATA_NAME   = FILE_DATA[jFILE][2]
				TEXT_NAME   = TEXT_NAMES[iNAME]
#
				if DATA_SOURCE == 'JULES_Wetlands' and os.path.exists(FILENAME):
					print(FILENAME)
					DIMS,EMISS_IN= data_netCDF.data_netCDF_array_var(FILENAME,DATA_NAME)
					EMISS_IN     = np.squeeze(EMISS_IN)
#
# On first pass, need to define multi-year emission array
#
					print YEAR,iREAD,iTIME
					if iREAD == 0:
						iREAD       = 1
						NTIMES_JUL  = (min(END_YEAR,END_DATA)-max(START_YEAR,START_DATA)+1) \
							*NMONTHS
						print NTIMES_JUL
						EMISS       = np.zeros((NTIMES_JUL,EMISS_IN.shape[1],EMISS_IN.shape[2]))
					if PLOTS[5] == 'Y':
						DATA_NAME    = DATA_NAME+'_wet'
						WET_FRACT_ANN = WET_FRACT[wTIME:wTIME+12,:,:]
						print wTIME,WET_FRACT_ANN.shape,EMISS.shape
						INDICES      = WET_FRACT_ANN[:,:,:] > 0
						print len(INDICES[INDICES]),WET_FRACT_ANN.size
#
						EMISS_WET    = np.zeros(EMISS_IN.shape)
						EMISS_WET[:,:,:] = 0.0
						EMISS_WET[INDICES] = EMISS_IN[INDICES]/ \
							WET_FRACT_ANN[INDICES]
						WET_FRACT_ANN[WET_FRACT_ANN <= 0.0] = MISS_DATA
#
						if DEBUG == 'Y':
							print oLAT,oLONG
							print WET_FRACT_ANN[:,oLAT,oLONG]
							print EMISS_IN[:,oLAT,oLONG]
							print EMISS_WET[:,oLAT,oLONG]
#
					if LONG_END == 360.0:
						EMISS_IN    = plot_map.switch_long_time(EMISS_IN)
						if PLOTS[5] == 'Y':
							EMISS_WET     = plot_map.switch_long_time(EMISS_WET)
							WET_FRACT_ANN = plot_map.switch_long_time(WET_FRACT_ANN)
#
					EMISS[iTIME:iTIME+NMONTHS,:,:] = EMISS_IN
#
				iTIME        += NMONTHS
#
		if len(EMISS) != 0:
			EMISS[EMISS < 0] = 0.0
		else:
			EMISS            = np.zeros((NTIMES,AREA.shape[0],AREA.shape[1]))
#
		PLOT_DIR       = NETCDF_DIR+'a_WETLANDS_SUMMARY/'
#
		LEGEND_JULES   = (SOURCE_ALL[iFILE]+MODIFIER[jFILE][1]+OPT_CORR[iOPT_CORR][2]).replace('_',' ')
		LEGENDS_MULTI.append(LEGEND_JULES)
#
	else:
		PLOTS          = list(PLOTS)
		PLOTS[5]       = 'N'
		PLOTS          = ''.join(PLOTS)
#
		LEGENDS_MULTI.append(SOURCE_ALL[iFILE])
#
		if iUSE_STORED == 'N':
			print('Input number of netCDF files: ')
			NSUB_FILES     = int(input())
		else:
			NSUB_FILES     = int(DATA_STORED[dSOURCE][1])
#
		FILE_CDF       = []
#
		for iSUB_FILE in range(NSUB_FILES):
			if iUSE_STORED == 'N':
				print('Input netCDF filename: '+NETCDF_DIR)
				FILENAME       = NETCDF_DIR+input()
			else:
				FILENAME       = NETCDF_DIR+DATA_STORED[dSOURCE][2][iSUB_FILE]
#
			print(FILENAME)
			FILE_CDF.append(FILENAME)
#
		PLOT_DIR       = NETCDF_DIR
#
# Input data
#
# Get Variable Names
#
		if GRIDDED == 'LAND_MON':
			SDATE         = str(START_DATA)+'01'
			FILE_VAR      = FILENAME.replace('XXXXXX',SDATE)
		else:
			FILE_VAR      = FILENAME
#
		print FILE_VAR
		VAR_NAMES=data_netCDF.data_netCDF_getVARNAMES(FILE_VAR)
#
		if iUSE_STORED == 'N':
			for iVAR in range(len(VAR_NAMES)):
				TEXT         = ('%4d %s' % (iVAR,VAR_NAMES[iVAR]))
				print(TEXT)
# Select variable
			print('Input index for variable: ')
			sVAR         = str(input()).split(':')
		else:
			sVAR         = DATA_STORED[dSOURCE][3].split(':')
#
	if SRESOL=='UM':
		RESOL_LONG     = 1.875
		RESOL_LAT      = 1.250
		NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
		NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
		DATA_IN        = np.zeros((NTIMES_1,NLAT,NLONG))
		LONG_ALL       = LONG_PLOTS+RESOL_LONG*arange(NLONG)
		LAT_ALL        = LAT_START +RESOL_LAT/2.0 +RESOL_LAT *arange(NLAT)
	elif SRESOL=='MACC':
		RESOL_LONG     = 3.75
		RESOL_LAT      = 2.50
		NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
		NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
		DATA_IN        = np.zeros((NTIMES_1,NLAT,NLONG))
		LONG_ALL       = LONG_PLOTS+RESOL_LONG*arange(NLONG)
		LAT_ALL        = LAT_START +RESOL_LAT/2.0 +RESOL_LAT *arange(NLAT)
	else:
		RESOL_LONG     = float(SRESOL)
		RESOL_LAT      = float(SRESOL)
		NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
		NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
		DATA_IN        = np.zeros((NTIMES_1,NLAT,NLONG))
		LONG_ALL       = LONG_PLOTS+RESOL_LONG/2.0+RESOL_LONG*arange(NLONG)
		LAT_ALL        = LAT_START +RESOL_LAT/2.0 +RESOL_LAT *arange(NLAT)
#
	iDOM_LONGS  = int((LONG_DOMS-LONG_PLOTS)/RESOL_LONG) 
	iDOM_LONGE  = int((LONG_DOME-LONG_PLOTS)/RESOL_LONG)
	iDOM_LATS   = int((LAT_DOMS-LAT_START)/RESOL_LAT) 
	iDOM_LATE   = int((LAT_DOME-LAT_START)/RESOL_LAT)
#
	iDOM_LONGSR = int((LONG_DOMS-LONG_PLOTS)/O_RESOL_LONG) 
	iDOM_LONGER = int((LONG_DOME-LONG_PLOTS)/O_RESOL_LONG)
	iDOM_LATSR  = int((LAT_DOMS-LAT_START)/O_RESOL_LAT) 
	iDOM_LATER  = int((LAT_DOME-LAT_START)/O_RESOL_LAT)
#
	print iDOM_LONGS,iDOM_LONGE,iDOM_LATS,iDOM_LATE, \
		iDOM_LONGSR,iDOM_LONGER,iDOM_LATSR,iDOM_LATER
#
	LONG       = np.arange(LONG_DOMS,LONG_DOME+DEL_LONG,DEL_LONG)
	LAT        = np.arange(LAT_DOMS,LAT_DOME+DEL_LAT,DEL_LAT )
#
	LONG_CSAT  = LONG_DOMS+O_RESOL_LONG*(0.5+arange(NLONG_AREA))
	LAT_CSAT   = LAT_DOMS +O_RESOL_LAT*(0.5+arange(NLAT_AREA))
#
	LONG_MAP   = LONG_ALL[(LONG_ALL >= LONG_DOMS ) & (LONG_ALL <= LONG_DOME)]
	LAT_MAP    = LAT_ALL[(LAT_ALL >= LAT_DOMS)     & (LAT_ALL <= LAT_DOME)]
#
	if DEBUG == 'Y':
		print LONG_DOMS,LONG_DOME,LAT_DOMS,LAT_DOME
		print LONG
		print LAT
		print LONG_ALL
		print LAT_ALL
		print LONG_MAP
		print LAT_MAP
		print LONG_CSAT
		print LAT_CSAT
#
	if PLOTS[11] == 'Y':
#
		MLAT           = iDOM_LATE -iDOM_LATS
		MLONG          = iDOM_LONGE-iDOM_LONGS
		MTIMES         = NYEARS*NMONTHS
#
		NDIMS_KMZ.append([MTIMES,MLAT,MLONG])
		LATS_KMZ.append([LAT_DOMS, LAT_DOME ])
		LONGS_KMZ.append([LONG_DOMS,LONG_DOME])
#
# Emissions will be from 180W to 180 E, hence use LONG_PLOTS
#
	if PLOTS[16] == 'Y':
		LATS_ALL.append( LAT_START +RESOL_LAT/2.0 +RESOL_LAT *arange(NLAT))
		LONGS_ALL.append(LONG_PLOTS+RESOL_LONG/2.0+RESOL_LONG*arange(NLONG))
#
	ANN_EMISS_FILE= np.zeros((NYEARS,NLAT,NLONG))
	MON_EMISS_FILE= np.zeros((NTIMES,NLAT,NLONG))
	MON_EMISS     = np.zeros((NLAT,NLONG))
	TOTAL_EMISS   = np.zeros((NLAT,NLONG))
	TOTAL_EMISS_1 = np.zeros((NLAT,NLONG))
	TOTAL_EMISS_2 = np.zeros((NLAT,NLONG))
	DATA_FRACT    = np.zeros((NLAT,NLONG))
	NUM_IN        = np.ones((NTIMES_1,NLAT,NLONG))
#
	SOPTIONS      = [ 'Annual in months',  \
			  'Annual in seconds', \
			  'Monthly in months',  \
			  'Monthly in seconds',
			  'Monthly in hours',
			  'Monthly in days' ]
#
	if iUSE_STORED == 'N':
		print('Input scaling factor: ')
		FACTOR        = float(input())
#
		print('Input N for numpy array or M for masked array: ')
		iMASK         = input()
#
		for iOPT in range(len(SOPTIONS)):
			TEXT         = ('%4d %s' % (iOPT,SOPTIONS[iOPT]))
			print(TEXT)
#
		print('Input timebase for variable: ')
		iOPT          = int(input())
#
	else:
		FACTOR        = float(DATA_STORED[dSOURCE][4])
		iMASK         = DATA_STORED[dSOURCE][5]
		iOPT          = int(DATA_STORED[dSOURCE][6])
#
	if PLOTS[3] == 'Y':
#
		if iFILE == 0:
#
			NLAT_INT2      = 0
			while NLAT_INT2 == 0:
				print('Input latitude interval: ')
				dLAT           = float(input())
				NLAT_INT       = int(180/dLAT)
				NLAT_INT2      = int(dLAT/RESOL_LAT)
#
			LAT_ZONAL      = LAT_START + (dLAT/2.0) + dLAT*arange(NLAT_INT,dtype='float32')
			LAT_TICKS      = int(LAT_START)+10*arange(19)
#
			EMISS_LAT      = np.zeros((NFILES,NYEARS,NLAT_INT))
			XZONAL         = np.zeros((NFILES,NLAT_INT))
			YZONAL         = np.zeros((NFILES,NLAT_INT))
#
		else:
			NLAT_INT2      = int(dLAT/RESOL_LAT)
#
	if PLOTS[4] == 'Y':
#
		if iFILE == 0:
#
			NLON_INT2      = 0
			while NLON_INT2 == 0:
				print('Input longitude interval: ')
				dLON           = float(input())
				NLON_INT       = int(360/dLON)
				NLON_INT2      = int(dLON/RESOL_LONG)
#
			LON_MERID      = LONG_PLOTS + (dLON/2.0) + dLON*arange(NLON_INT,dtype='float32')
			LON_TICKS      = int(LONG_PLOTS)+30*arange(13)
#
			EMISS_LON      = np.zeros((NFILES,NYEARS,NLON_INT))
			XMERID         = np.zeros((NFILES,NLON_INT))
			YMERID         = np.zeros((NFILES,NLON_INT))
#
		else:
			NLON_INT2      = int(dLON/RESOL_LONG)
#
# Input lat,long for debug
#
	if DEBUG == 'Y' and PLOTS[5] == 'Y':
#
		print('Input grid square latitude  index: ')
		oLAT         = int(input())
#
		print('Input grid square longitude index: ')
		oLONG        = int(input())
#
# Need to have second longitude
#
		if LONG_END == 360.0:
			if oLONG < NLONG/2:
				oLONG2       = oLONG+NLONG/2
			else:
				oLONG2       = oLONG-NLONG/2
#
# Select threshold if PLOTS[13] == 'Y'
#
	if PLOTS[13] == 'Y':
		print('Input threshold: ')
		THRESHOLD    = float(input())
#
# Get data
#
	if not 'JULES_Wetlands' in DATA_SOURCE:
#
		for iSUB_FILE in range(NSUB_FILES):
#
			iVAR          = int(sVAR[iSUB_FILE])
			DATA_NAME     = VAR_NAMES[iVAR]
#
			FILENAME      = FILE_CDF[iSUB_FILE]
			FILE_PLOT_PART= FILENAME[:-3].replace('.','_')
			print GRIDDED
			if GRIDDED == 'GRID':
				DIMS,EMISS_IN = data_netCDF.data_netCDF_array_var(FILENAME,DATA_NAME)
			elif GRIDDED == 'LAND_MON':
				SDATE         = str(START_YEAR)+'_'+str(END_YEAR)
				FILE_PLOT_PART= FILE_PLOT_PART.replace('XXXXXX',SDATE)
				LAT_NAME      = 'latitude'
				LONG_NAME     = 'longitude'
				DIMS,EMISS_IN = data_netCDF.data_netCDF_array_land_var \
					(FILENAME,DATA_NAME,LAT_NAME,LONG_NAME, \
					 START_DATA,END_DATA,NLONG,NLAT,LONG_START,LAT_START,MISS_DATA,GRIDDED,DEBUG)
			EMISS_IN      = np.squeeze(EMISS_IN)
#
# Need to invert both latitude and longitude for Bloom inventory
#
			if 'Bloom' in FILENAME:
				EMISS_IN             = EMISS_IN[:,-1::-1,-1::-1]
				print EMISS_IN.shape
			if 'fch4.ref' in FILENAME:
				SOURCE_ALL[iFILE]    = SOURCE_ALL[iFILE]+' Fung'
				LEGENDS_MULTI[iFILE] = LEGENDS_MULTI[iFILE]+' Fung'
				FILE_PLOT_PART       = FILE_PLOT_PART+'_'+DATA_NAME
			elif 'fch4.kaplan' in FILENAME:
				SOURCE_ALL[iFILE]    = SOURCE_ALL[iFILE]+' Kaplan'
				LEGENDS_MULTI[iFILE] = LEGENDS_MULTI[iFILE]+' Kaplan'
				FILE_PLOT_PART       = FILE_PLOT_PART+'_'+DATA_NAME
			elif 'fung' in FILENAME:
				if 'rice' in FILENAME:
					SOURCE_ALL[iFILE]    = SOURCE_ALL[iFILE]+'+Rice'
					LEGENDS_MULTI[iFILE] = LEGENDS_MULTI[iFILE]+'+Rice'
			elif 'A06' in FILENAME:
				SOURCE_ALL[iFILE]    = SOURCE_ALL[iFILE]+' JULES'
				LEGENDS_MULTI[iFILE] = LEGENDS_MULTI[iFILE]+' JULES'
			elif 'C06' in FILENAME:
				SOURCE_ALL[iFILE]    = SOURCE_ALL[iFILE]+' JULES+GIEMS'
				LEGENDS_MULTI[iFILE] = LEGENDS_MULTI[iFILE]+' JULES+GIEMS'
			print SOURCE_ALL[iFILE]
#
			if 'M' in iMASK:
				import numpy.ma as ma
				if iMASK == 'M':
					MASK_FACTOR          = 1.0
				else:
					MASK_FACTOR          = float(iMASK.split('_')[1])
				TEMP     = EMISS_IN
				EMISS_IN = MASK_FACTOR*ma.getdata(TEMP)
				EMISS_IN[EMISS_IN == MISS_DATA] = 0.0
				EMISS_IN[EMISS_IN != MISS_DATA] = EMISS_IN[EMISS_IN != MISS_DATA]/MASK_FACTOR
				print TEMP.min(),TEMP.max(),TEMP.sum(),EMISS_IN.min(),EMISS_IN.max(),EMISS_IN.sum()
			else:
				print EMISS_IN.min(),EMISS_IN.max(),EMISS_IN.sum()
#
# Set missing data to zero
#
			EMISS_IN[(EMISS_IN == MISS_DATA) | (EMISS_IN < 0)] = 0.0
#
# Need to invert latitude dimension if DATA_SOURCE = 'Inverse_N2O'
#
			if DATA_SOURCE == 'Inverse_N2O':
				EMISS_IN = EMISS_IN[:,::-1,:]
#
# Need to correct latitude array if SRESOL is UM or MACC
#
			if SRESOL == 'UM' or SRESOL == 'MACC':
				EMISS_IN      = EMISS_IN[:,:-1,:]
#
			if iSUB_FILE == 0:
				if len(EMISS_IN.shape) == 2:
					EMISS            = np.zeros((NYEARS,EMISS_IN.shape[0],EMISS_IN.shape[1]))
					for iYEAR in range(NYEARS):
						EMISS[iYEAR,:,:] = EMISS_IN
				else: 
					EMISS            = EMISS_IN
			else:
				if len(EMISS_IN.shape) == 2:
					for iYEAR in range(NYEARS):
						EMISS[iYEAR,:,:] = EMISS[iYEAR,:,:] + EMISS_IN
				else:
					EMISS            = EMISS + EMISS_IN
#
		if LONG_END == 360.0:
			EMISS         = plot_map.switch_long_time(EMISS)
			print 'Emission dataset: E-W hemispheres switched'
#
#
	if 'JULES' in DATA_SOURCE and (PLOTS[11] == 'Y' or PLOTS[14] == 'Y'):
#
# Switch E-W hemisphere as needed
#
		if LONG_END == 360.0:
			WET_FRACT = plot_map.switch_long_time(WET_FRACT)
#
		WET_FRACT_KMZ.append(WET_FRACT[wTIME:wTIME+NYEARS*NMONTHS,iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE])
#
# Calculate annual emissions
#
	if   iOPT      == 0:
		NSECS       = 1.0
#		NSECS       = 24.0*60.0*60.0*365.25/12.0
		NSTEP       = 1
	elif iOPT      == 1:
		NSECS       = 24.0*60.0*60.0*365.25
		NSTEP       = 1
	elif iOPT      == 2:
		NSECS       = 1.0
		NSTEP       = 12
	elif iOPT      == 3 or iOPT      == 6:
		NSECS       = 24.0*60.0*60.0
		NSTEP       = 12
	elif iOPT      == 4:
		NSECS       = 24.0
		NSTEP       = 12
	elif iOPT      == 5:
		NSECS       =  1.0
		NSTEP       = 12
#
	print(NSECS)
#
# Assign to master emission array
#
	EMISS_ALL  = np.zeros((NTIMES,EMISS.shape[1],EMISS.shape[2]))
#
	iTIME      = 0
	for iYEAR in range(NYEARS):
#
		YEAR               = START_YEAR+iYEAR
#
		if EMISS.shape[0] == 12:
#
# If single year, replicate this for other years
#
			EMISS_ALL[iTIME:iTIME+NMONTHS] = EMISS[:,:,:]
		else:
#
# Identify array elements and assign annual emissions
#
			iSTART_ALL = iTIME
			iEND_ALL   = iSTART_ALL+NSTEP
#
			if YEAR >= START_DATA and YEAR <= END_DATA:
#
				iSTART_DAT = (YEAR-START_DATA)*NSTEP
				iEND_DAT   = iSTART_DAT+NSTEP
				print iYEAR,YEAR,iTIME,iSTART_ALL,iEND_ALL,iSTART_DAT,iEND_DAT
#
				EMISS_ALL[iSTART_ALL:iEND_ALL,:,:] = \
					EMISS[iSTART_DAT:iEND_DAT,:,:]
#
		iTIME += NSTEP
#			
	iTIME      = 0
	jTIME      = 0
	iSTART     = START_YEAR-START_DATA
	iEND       = iSTART+(END_YEAR-START_YEAR)
	dTIME      = iSTART
#
	for iYEAR in range(NYEARS):
#
		YEAR               = START_YEAR+iYEAR
		SYEAR              = '%4d' % (YEAR)
#
# Check for leap year
#
		if ('Bloom' in DATA_SOURCE or iLEAP == 'Y') and YEAR == 4*int(YEAR/4):
			DAYS_MONTH[2]      = 29
		else:
			DAYS_MONTH[2]      = 28
#
		print YEAR,DAYS_MONTH[2]
#
		TOTAL_EMISS[:,:]   = 0.0
		TOTAL_EMISS_1[:,:] = 0.0
		TOTAL_EMISS_2[:,:] = 0.0
#
		MON_EMISS_2        = np.zeros((NMONTHS,NLAT,NLONG))
		print MON_EMISS_2.shape
		print iYEAR,NYEARS,iTIME,jTIME,dTIME
#
		if iOPT == 0 or iOPT == 1: # Annual dataset
			CONV_FACTOR  = FACTOR*NSECS/1.0E+09
		elif DATA_SOURCE == 'JULES_Wetlands_Area':
			CONV_FACTOR  = FACTOR
#
		if iOPT == 0 or iOPT == 1 or DATA_SOURCE == 'JULES_Wetlands_Area':
#
# Convert from kg m-2 s-1 to Tg per annum
#
			ANN_EMISS_FILE[iYEAR,:,:] = EMISS_ALL[iYEAR,:,:]*AREA[:,:]*CONV_FACTOR
			TOTAL_EMISS[:,:]          = TOTAL_EMISS[:,:] + ANN_EMISS_FILE[iYEAR,:,:]
#
			TOTAL_EMISS_1[:,:]  = TOTAL_EMISS_1[:,:] + \
				EMISS_ALL[iYEAR,:,:]*1.00E+12*CONV_FACTOR
#
			TOTAL_EMISS_2[:,:]  = TOTAL_EMISS_2[:,:] + \
				EMISS_ALL[iYEAR,:,:]*1.00E+12*CONV_FACTOR
#
		if iOPT >= 2: # Loop over months
#
# Reset iTIME if only single year present in datafile - replicate this annual cycle
#
			if EMISS.shape[0] == 12:
				iTIME        = 0
#
			for iMONTH in range(NMONTHS):
#
				if iOPT == 2:   CONV_FACTOR  = FACTOR*NSECS/1.0E+09
				elif iOPT >= 3: CONV_FACTOR  = FACTOR*NSECS*DAYS_MONTH[iMONTH+1]/1.0E+09
#				elif iOPT == 4: CONV_FACTOR  = FACTOR*NSECS*30.0/1.0E+09
				elif DATA_SOURCE == 'JULES_Wetlands_Area':
					CONV_FACTOR  = FACTOR
#
# Convert from kg m-2 s-1 to Tg per annum
#
				if iOPT == 6: # Daily data
					print iMONTH+1,YEAR,DAYS_MONTH[iMONTH+1],dTIME
					MON_EMISS[:,:]      = 0.0
					for iDAY in range(DAYS_MONTH[iMONTH+1]):
						MON_EMISS[:,:]      = MON_EMISS[:,:] \
							+EMISS_ALL[dTIME,:,:]*AREA[:,:] \
							*CONV_FACTOR/DAYS_MONTH[iMONTH+1]
						dTIME += 1
				else:
					MON_EMISS[:,:]      = EMISS_ALL[iTIME,:,:]*AREA[:,:]*CONV_FACTOR
#
				MON_EMISS_FILE[jTIME,:,:] = MON_EMISS
				TOTAL_EMISS               = TOTAL_EMISS + MON_EMISS
#
				if PLOTS[5] != 'Y':
					TOTAL_EMISS_1[:,:]  = TOTAL_EMISS_1[:,:] + \
						EMISS_ALL[iTIME,:,:]*CONV_FACTOR*1.00E+12
				else:
					TOTAL_EMISS_1[:,:]  = TOTAL_EMISS_1[:,:] + \
						EMISS_WET[iTIME,:,:]*CONV_FACTOR*1.00E+12	
#
					if DEBUG == 'Y':
						print iTIME,oLONG2,oLAT,EMISS_WET[iTIME,oLAT,oLONG2],CONV_FACTOR, \
							EMISS_WET[iTIME,oLAT,oLONG2]*CONV_FACTOR*1.00E+12, \
							TOTAL_EMISS_1[oLAT,oLONG2]
#
				if PLOTS[0] == 'Y' or PLOTS[13] == 'Y':
					print EMISS_ALL.shape,MON_EMISS_2.shape
					MON_EMISS_2[iMONTH,:,:] = EMISS_ALL[iTIME,:,:]*CONV_FACTOR*1.00E+15/(DAYS_MONTH[iMONTH+1])
					TOTAL_EMISS_2[:,:]  = TOTAL_EMISS_2[:,:] + MON_EMISS_2[iMONTH,:,:]/12.0
#					TOTAL_EMISS_2[:,:]  = TOTAL_EMISS_2[:,:] + \
#						EMISS_ALL[iTIME,:,:]*CONV_FACTOR*1.00E+15/(12*30.0)
#
				if DEBUG == 'Y':
					print iMONTH,CONV_FACTOR,EMISS_ALL[iTIME,:,:].min(),EMISS_ALL[iTIME,:,:].max(),EMISS_ALL[iTIME,:,:].sum(),TOTAL_EMISS.sum()
#
# Data for climatology of annual emissions
#
				if (PLOTS[1] == 'Y' or PLOTS[7] == 'Y' or PLOTS[10] == 'Y') and iOPT >= 2:
#
					for iTRANS in range(NTRANS):
#
						if iTRANS == 0:
							TRANS_EMISS_MN[iFILE,iTRANS,jTIME] = MON_EMISS.sum()
						else:
							TRANS_EMISS_MN[iFILE,iTRANS,jTIME] = MON_EMISS[(TRANSCOM == iTRANS)].sum()
#
						if TRANS_EMISS_MN[iFILE,iTRANS,jTIME] > 0:
							TRANS_EMISS_CL[iMONTH,iTRANS] = TRANS_EMISS_CL[iMONTH,iTRANS] + \
                                                       		TRANS_EMISS_MN[iFILE,iTRANS,jTIME]
							TRANS_EMISS_VL[iMONTH,iTRANS] = TRANS_EMISS_VL[iMONTH,iTRANS] + 1
#
				if PLOTS[9] == 'Y' or PLOTS[12] == 'Y':
					if MAP_FRACT == 'Y':
						DATA_MONTH_ALL.append(WET_FRACT_ANN[iTIME,iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE])
						DATA_MONTH_ALL.append(EMISS_WET[iTIME,iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE]*CONV_FACTOR*1.00E+12)
					else:
						DATA_MONTH_ALL.append(EMISS_ALL[iTIME,iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE]*CONV_FACTOR*1.00E+12)
#
# Google Earth - PLOTS[11]
#
				if PLOTS[11] == 'Y' and iOPT >= 2:
#
					SDATE_TIME_BEG = '%4d-%02d-%02dT%02d:%02d:%02d' % (iYEAR+START_YEAR,iMONTH+1,1,0,0,0)
					SDATE_TIME_END = '%4d-%02d-%02dT%02d:%02d:%02d' % (iYEAR+START_YEAR,iMONTH+1,DAYS_MONTH[iMONTH+1],23,59,59)
#
					SDATE_TIMES_KMZ[iFILE].append(SDATE_TIME_BEG)
					EDATE_TIMES_KMZ[iFILE].append(SDATE_TIME_END)
#
					print iTIME,iDOM_LATS,iDOM_LATE,iDOM_LONGS,iDOM_LONGE
					EMISS_KMZ.append(EMISS_ALL[iTIME,iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE]*CONV_FACTOR*1.00E+15/DAYS_MONTH[iMONTH+1])
#
				if PLOTS[15] == 'Y':
					EMISS_DOMAIN_M[jTIME,iFILE] = MON_EMISS_FILE[jTIME,iDOM_LATS:iDOM_LATE+1,iDOM_LONGS:iDOM_LONGE+1].sum()
#
				iTIME        += 1
				jTIME        += 1
#
# Data for time series of annual emissions
#
		EMISS_ANNUAL[iYEAR,iFILE] = TOTAL_EMISS[:,:].sum()
		EMISS_DOMAIN[iYEAR,iFILE] = TOTAL_EMISS[iDOM_LATS:iDOM_LATE+1,iDOM_LONGS:iDOM_LONGE+1].sum()
#
		for iTRANS in range(NTRANS):

			if iTRANS == 0:
				TRANS_EMISS[iYEAR,iTRANS] = TOTAL_EMISS[:,:].sum()
			else:
				TRANS_EMISS[iYEAR,iTRANS] = TOTAL_EMISS[(TRANSCOM == iTRANS)].sum()
#
	FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME)
#
# Convert zero to Nan or missing data
#

	if iOPT >= 2:
		MON_EMISS_ALL.append(MON_EMISS_FILE)
		NTIMES_ALL.append(MON_EMISS_FILE.shape[0])
	else:
		MON_EMISS_ALL.append(ANN_EMISS_FILE)
		NTIMES_ALL.append(NYEARS)
#
	if 'JULES_Wetlands' in DATA_SOURCE:
		NAMES_KMZ.append(FILE_DATA[jFILE][3][1:]+FILE_CORR+SDOMAIN+str(START_YEAR)+'_'+str(END_YEAR))
	else:
		NAMES_KMZ.append(SOURCE_ALL[iFILE].replace('Wetlands','Wetland_Emissions')+SDOMAIN+str(START_YEAR)+'_'+str(END_YEAR))
#
	COLOURS_KMZ.append(COLOURS)
	CLEVELS_KMZ.append(CLEVELS)
#
	if PLOTS[11] == 'Y' and iOPT >= 2:
		UNITS_KMZ.append('CH$_{4}$ emission flux (mgCH$_{4}$ m$^{-2}$ day$^{-1}$)')
	else:
		UNITS_KMZ.append(UNITS)
#
# Output time series of annual emissions
#
	for iYEAR in range(NYEARS):
#
		iSTART2      = iYEAR*NMONTHS
		iEND2        = iSTART2+NMONTHS
#
		SYEAR        = str(START_YEAR+iYEAR)
		TEXT         = 'Global annual emissions for '+SYEAR+' = '+ \
			('%12.4f' % (EMISS_ANNUAL[iYEAR,iFILE]))+' Tg per annum'
		print(TEXT)
		TEXT         = 'Domain annual emissions for '+SYEAR+' = '+ \
			('%12.4f' % (EMISS_DOMAIN[iYEAR,iFILE]))+' Tg per annum'
		print(TEXT)
		print TRANS_EMISS[iYEAR,:],TRANS_EMISS[iYEAR,1:].sum()
#
		if PLOTS[15] == 'Y':
			TEXT         = 'Domain monthly emissions for '+SYEAR+' = '
			for iMONTH in range(NMONTHS):
				TEXT         = TEXT+('%12.4f' % (EMISS_DOMAIN_M[iSTART2+iMONTH,iFILE]))
			TEXT         = TEXT+('%12.4f' % (EMISS_DOMAIN_M[iSTART2:iEND2,iFILE].sum()))+' Tg'
			print TEXT
#
# Plots
# 
		if PLOTS[0] == 'Y' or PLOTS[8] == 'Y' or PLOTS[9] == 'Y' or PLOTS[12] == 'Y' or PLOTS[13] == 'Y':
#
			if iMASS == '1':
				if FSPECIES == 'NOx':
					SCALE_FACTOR = 14.0/RMM
					UNITS2       = 'gN m$^{-2}$ yr$^{-1}$' 
				elif FSPECIES == 'N2O':
					SCALE_FACTOR = 1.00E-03/12.0
					UNITS2       = 'gN m$^{-2}$ yr$^{-1}$' 
				else:
					SCALE_FACTOR = 12.0/RMM
					UNITS2       = 'gC m$^{-2}$ yr$^{-1}$'
			elif iMASS == '2' and FSPECIES == 'N2O':
				SCALE_FACTOR = 1.00E-03*10.0*RMM/(12.0*14.0)
				UNITS2       = 'kg ha$^{-1}$ yr$^{-1}$'
				DATA_NAME    = 'n2oflux_ha' 
				FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME)
			else:
				SCALE_FACTOR = 1.0
				UNITS2       = UNITS
#
			print SCALE_FACTOR
			print UNITS2
#
# PLOTS[0]: Annual Maps
#
			DATA_MAP         = TOTAL_EMISS_1[iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE]*SCALE_FACTOR
			DATA_MAP[DATA_MAP <=0] = MISS_DATA
			DATA_ANN_ALL.append(DATA_MAP)
			print SYEAR,DATA_MAP.min(),DATA_MAP.max()
#
			if DATA_SOURCE == 'TRANSCOM_Fung':
				FILE_PLOT  = FILENAME[0:50]+'ch4_actm_wetlands_'+FSPECIES+SDOMAIN+SYEAR+'.png'
			else:
				FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'.png'
			print FILE_PLOT
			PLOT_LABEL = 'Emissions of '+CSPECIES+' ('+UNITS2+')'
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Annual emission of '+CSPECIES+' for '+SYEAR
#
			print DATA_MAP.shape,LAT_MAP.shape,LONG_MAP.shape,LAT.shape,LONG.shape
#
			plot_map.plot_map3(DATA_MAP,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
				LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER,SET_OVER,DEBUG)
#
		if PLOTS[8] == 'Y' and MAP_FRACT == 'Y':
#
			HEIGHT       = HEIGHT0*2
			WIDTH        = WIDTH0
			NROWS        = 2
			NCOLUMNS     = 1
#
			SET_UNDER_ALL= [SET_UNDER,SET_UNDER]
			SET_OVER_ALL = [SET_OVER,SET_OVER]
			CLEVELS_ALL  = [CLEVELS_C,CLEVELS]
			COLOURS_ALL  = [COLOURS_C,COLOURS]
			MAP_TYPES    = [MAP_TYPE,MAP_TYPE]
			DATA_MAP_ALL = []
#
			print 'Max wetland fraction'
			for iLAT in range(NLAT):
				for iLONG in range(NLONG):
					DATA_FRACT[iLAT,iLONG] = np.max(WET_FRACT_ANN[:,iLAT,iLONG])
			print 'Max wetland fraction'
#
			DATA_MAP_ALL.append(DATA_FRACT[iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE])
			DATA_MAP_ALL.append(DATA_MAP)
#
			PLOT_LABELS = [ 'Fraction of grid-square', \
					'Emissions of '+CSPECIES+' ('+UNITS2+')' ]
					
			SUB_TITLES  = [ 'Maximum grid-square wetland fraction', \
					'Global annual emission' ]
#
			PLOT_TITLE  = 'Wetland fraction and emissions of '+CSPECIES+' for '+SYEAR
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+'_wetlands'+SDOMAIN+SYEAR+'.png'
			print FILE_PLOT
#
			plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_MAP_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
				LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
				RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
# PLOTS[0]: Single set of monthly maps
#
		if PLOTS[0] == 'Y' or PLOTS[13] == 'Y':
#
			if DATA_NAME in DATA_NAMES:
				DATA_NAME_R = 'ch4_surf_emiss_base'
				FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME_R)
#
			SRESOL_ALL = [ SRESOL,SRESOL_OUT ]
			INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT= \
				data_regrid_new.data_regrid_factor(DEBUG,2,SRESOL_ALL,LONG_START,LAT_START,LONG_END,LAT_END)
#
			DATA_MAP         = TOTAL_EMISS_2[iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE]
			DATA_MAP[DATA_MAP <=0] = MISS_DATA
			print SYEAR,DATA_MAP.min(),DATA_MAP.max()
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_base.png'
			print FILE_PLOT
			UNITS2     = 'mgCH$_{4}$ m$^{-2}$ d$^{-1}$'
			PLOT_LABEL = 'Emissions of '+CSPECIES+' ('+UNITS2+')'
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			plot_map.plot_map3(DATA_MAP,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
				LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER,SET_OVER,DEBUG)
#
			DATA_ALL     = []
			DATA_R_ALL   = []
			DATA_T_ALL   = []
			CLEVELS_ALL  = []
			COLOURS_ALL  = []
			MAP_TYPES    = []
			PLOT_LABELS  = []
			SUB_TITLES   = SMONTHS
			SET_UNDER_ALL= []
			SET_OVER_ALL = []
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_monthly_base.png'
			print FILE_PLOT
			PLOT_LABEL = 'Emissions of '+CSPECIES+' ('+UNITS2+')'
#
			for iMONTH in range(NMONTHS):
#
				SDATE        = ('%4d%02d' % (iYEAR+START_YEAR,iMONTH+1))
#
				DATA_MON         = MON_EMISS_2[iMONTH,:,:]
				DATA_MON[DATA_MON <=0] = MISS_DATA
				print DATA_MON.max(),DATA_MON[DATA_MON >=0].min()
#
#				if SRESOL == 'UM' or SRESOL == 'MACC':
#					NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
#					DATA_IN[0,:,:] = MON_EMISS_2[iMONTH,:-1,:]
#				else:
#					DATA_IN[0,:,:] = MON_EMISS_2[iMONTH,:,:]
#
				DATA_IN[0,:,:] = MON_EMISS_2[iMONTH,:,:]
#
				DATA_ALL.append(DATA_MON[iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE])
				DATA_MONTH_ALL2.append(DATA_MON[iDOM_LATS:iDOM_LATE,iDOM_LONGS:iDOM_LONGE])
#
				if PLOTS[13] == 'Y':
					DATA_MON_R=data_regrid_new.data_regrid(DEBUG,DATA_IN,NUM_IN,NTIMES_1,MISS_DATA, \
						NLONG,NLAT,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
					DATA_MON_R       = np.squeeze(DATA_MON_R)
					DATA_MON_R[DATA_MON_R <=0] = MISS_DATA
					print DATA_MON_R.max(),DATA_MON_R[DATA_MON_R >=0].min()
					DATA_R_ALL.append(DATA_MON_R[iDOM_LATSR:iDOM_LATER,iDOM_LONGSR:iDOM_LONGER])
#
					DATA_MON_T       = np.copy(DATA_MON_R)
					DATA_MON_T[DATA_MON_R > 0] = DATA_MON_R[DATA_MON_R > 0]/THRESHOLD
					print DATA_MON_T.max(),DATA_MON_T[DATA_MON_T >=0].min()
					DATA_T_ALL.append(DATA_MON_T[iDOM_LATSR:iDOM_LATER,iDOM_LONGSR:iDOM_LONGER])
#
				SET_UNDER_ALL.append(SET_UNDER)
				SET_OVER_ALL.append(SET_OVER)
				CLEVELS_ALL.append(CLEVELS)
				COLOURS_ALL.append(COLOURS)
				PLOT_LABELS.append(PLOT_LABEL)
				MAP_TYPES.append(MAP_TYPE)
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			NROWS        = 6
			NCOLUMNS     = 2
			HEIGHT       = HEIGHT0*NROWS
			WIDTH        = WIDTH0*NCOLUMNS
#
			plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
				LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
				RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
# Regrid to CarbonSAT spatial domain
#
		if PLOTS[13] == 'Y':
#
			if SRESOL == 'UM' or SRESOL == 'MACC':
				NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
				DATA_IN[0,:,:] = TOTAL_EMISS_2[:-1,:]
			else:
				DATA_IN[0,:,:] = TOTAL_EMISS_2[:,:]
#
			DATA_MAP_R=data_regrid_new.data_regrid(DEBUG,DATA_IN,NUM_IN,NTIMES_1,MISS_DATA, \
				NLONG,NLAT,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
			DATA_MAP_R       = np.squeeze(DATA_MAP_R)
			DATA_MAP         = np.copy(DATA_MAP_R)
			DATA_MAP[DATA_MAP <=0] = MISS_DATA
			DATA_AREA        = DATA_MAP[iDOM_LATSR:iDOM_LATER,iDOM_LONGSR:iDOM_LONGER]
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_regrid.png'
			print FILE_PLOT
			PLOT_LABEL = 'Emissions of '+CSPECIES+' ('+UNITS2+')'
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			plot_map.plot_map3(DATA_AREA,LONG_DOMS,LONG_DOME,LONG,LONG_CSAT, \
				LAT_DOMS,LAT_DOME,LAT,LAT_CSAT,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER,SET_OVER,DEBUG)
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_monthly_regrid.png'
			print FILE_PLOT
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_R_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_CSAT, \
				LAT_DOMS,LAT_DOME,LAT,LAT_CSAT,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
				RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
# Repeat but ratioed to CarbonSAT threshold
#
			if DATA_NAME in DATA_NAMES:
				DATA_NAME_R = 'ch4_surf_emiss_thres'
				FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME_R)
#
			DATA_MAP         = np.copy(DATA_MAP_R)
			DATA_MAP[DATA_MAP > 0] = DATA_MAP[DATA_MAP > 0]/THRESHOLD
			DATA_AREA        = DATA_MAP[iDOM_LATSR:iDOM_LATER,iDOM_LONGSR:iDOM_LONGER]
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_thres_'+('%02d' % THRESHOLD)+'.png'
			print FILE_PLOT
			PLOT_LABEL = UNITS+' of '+('%2d' % THRESHOLD)+' mg CH$_{4}$ m$^{-2}$ d$^{-1}$'
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			print DATA_AREA.shape,LONG_MAP.shape,LAT_MAP.shape
			plot_map.plot_map3(DATA_AREA,LONG_DOMS,LONG_DOME,LONG,LONG_CSAT, \
				LAT_DOMS,LAT_DOME,LAT,LAT_CSAT,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER,SET_OVER,DEBUG)
#
# Set of Monthly maps
#
			FILE_PLOT  = FILE_PLOT_PART+'_'+FSPECIES+SDOMAIN+SYEAR+'_monthly_thres_'+('%02d' % THRESHOLD)+'.png'
			print FILE_PLOT
			PLOT_LABEL = UNITS+' of '+('%2d' % THRESHOLD)+' mg CH$_{4}$ m$^{-2}$ d$^{-1}$'
#
			CLEVELS_ALL  = []
			COLOURS_ALL  = []
			MAP_TYPES    = []
			PLOT_LABELS  = []
#
			for iMONTH in range(NMONTHS):
#
				CLEVELS_ALL.append(CLEVELS)
				COLOURS_ALL.append(COLOURS)
				PLOT_LABELS.append(PLOT_LABEL)
				MAP_TYPES.append(MAP_TYPE)
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				PLOT_TITLE = 'Global annual emission of '+CSPECIES
			else:
				PLOT_TITLE = 'Global annual emission of '+CSPECIES+' for '+SYEAR
#
			plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_T_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_CSAT, \
				LAT_DOMS,LAT_DOME,LAT,LAT_CSAT,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
				RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
# PLOTS[3]:
# Data for latitudinal plots
#
		if PLOTS[3] == 'Y':
#
			NDATA_ZONAL    = 1 
			TEMP           = Emissions_Wetlands_Zonal.Emissions_Wetlands_Zonal( \
				NDATA_ZONAL,NMONTHS,NLONG,NLAT_INT,NLAT_INT2,dLAT, \
				TOTAL_EMISS,MON_EMISS_FILE[iSTART2:iEND2,:,:],LAT_ZONAL,FSPECIES,SYEAR,DEBUG)
			EMISS_LAT[iFILE,iYEAR,:]=TEMP[0,:]
#
# PLOTS[4]:
# Data for longitudinal plots
#
		if PLOTS[4] == 'Y':
#
			NDATA_MERID    = 1
			TEMP           = Emissions_Wetlands_Zonal.Emissions_Wetlands_Meridional( \
				NDATA_MERID,NMONTHS,NLAT,NLON_INT,NLON_INT2,dLON,LONG_END,
				TOTAL_EMISS,MON_EMISS_FILE[iSTART2:iEND2,:,:],LON_MERID,FSPECIES,SYEAR,DEBUG)
			EMISS_LON[iFILE,iYEAR,:]=TEMP[0,:]
#
# PLOTS[14]:
# Data for site-specific plots
#
		if PLOTS[14] == 'Y':
#
			print iYEAR,MON_EMISS_FILE[iSTART2:iEND2,:,:].min(),MON_EMISS_FILE[iSTART2:iEND2,:,:].max(),len(EMISS_SITE)
			EMISS_SITE.append(MON_EMISS_FILE[iSTART2:iEND2,:,:])
#
			if iYEAR == 0:
#
				for iSITE in range(NSITES):
#
					SITE_NAME    = SITE_DATA[iSITE][0]
					SITE_LONG    = float(SITE_DATA[iSITE][1])
					SITE_LAT     = float(SITE_DATA[iSITE][2])
#
# Emissions will be from 180W to 180 E
#
					sLONG        = int((SITE_LONG-LONG_PLOTS)/RESOL_LONG)
					sLAT         = int((SITE_LAT -LAT_START )/RESOL_LAT)
					print(iSITE,SITE_NAME,SITE_LONG,SITE_LAT,sLONG,sLAT)
#
					RESOL_LONG_ALL.append(sLONG)
					RESOL_LAT_ALL.append(sLAT)
#
# Derive emission scale factors
#
	EMISS_FACTOR  = 1.0
	if DATA_SOURCE == 'JULES_Wetlands' and WET_TOTAL != -1.0:
		TEMP          = TRANS_EMISS[iSTART:iEND+1,0]
		if len(TEMP[TEMP > 0.0]) > 0:
			MISS_FACTOR  = WET_TOTAL/np.average(TEMP[TEMP > 0.0])
#
	TEMP          = EMISS_DOMAIN[:,iFILE]
#
	print 'Emission Scale Factor: ',EMISS_FACTOR
	print TRANS_EMISS[:,0]
	print TRANS_EMISS[:,0]*EMISS_FACTOR
	print TEMP
	print TEMP*EMISS_FACTOR
#
	if len(TEMP[TEMP>0.0]) > 0:
		print TEMP[TEMP>0.0].min()*EMISS_FACTOR,TEMP[TEMP>0.0].mean()*EMISS_FACTOR, \
			TEMP[TEMP>0.0].max()*EMISS_FACTOR
#
	EMISS_FACTORS.append(EMISS_FACTOR)
#
	TRANS_EMISS[TRANS_EMISS <= 0.0]       = float('NaN')
	TRANS_EMISS_CL[TRANS_EMISS_CL <= 0.0] = float('NaN')
	TRANS_EMISS_VL[TRANS_EMISS_VL <= 0.0] = float('NaN')
	TRANS_EMISS_MN[TRANS_EMISS_MN <= 0.0] = float('NaN')
#
# PLOTS[1]:
# Time series of emissions
	if PLOTS[1] == 'Y' or PLOTS[6] == 'Y' or PLOTS[7] == 'Y':
#
		NDATASETS     = 1
#
		XDATA         = np.zeros((NDATASETS,NYEARS))
		YDATA         = np.zeros((NDATASETS,NYEARS))
		XCLIM         = np.zeros((NDATASETS,NMONTHS))
		YCLIM         = np.zeros((NDATASETS,NMONTHS))
#
		for iTRANS in range(NTRANS):
#
			for iYEAR in range(NYEARS):
				XDATA[0,iYEAR]     = float(iYEAR+START_YEAR)+0.5
				YDATA[0,iYEAR]     = TRANS_EMISS[iYEAR,iTRANS]*EMISS_FACTOR
#
			XDATA_MULTI[iFILE,iTRANS,:] = XDATA[0,:]
			YDATA_MULTI[iFILE,iTRANS,:] = YDATA[0,:]
#
			for iMONTH in range(NMONTHS):
				XCLIM[0,iMONTH]    = float(iMONTH)+0.5
				YCLIM[0,iMONTH]    = TRANS_EMISS_CL[iMONTH,iTRANS]*EMISS_FACTOR \
					/TRANS_EMISS_VL[iMONTH,iTRANS]
#
			XCLIM_MULTI[iFILE,iTRANS,:] = XCLIM[0,:]
			YCLIM_MULTI[iFILE,iTRANS,:] = YCLIM[0,:]
#
			if DEBUG == 'Y':
				print iTRANS,YCLIM_MULTI[iFILE,iTRANS,:],YCLIM_MULTI[iFILE,iTRANS,:].sum(),YDATA_MULTI[iFILE,iTRANS,iYEAR]
#
			DOMAIN	    = TRANS_REGS[iTRANS]
			sTRANS      = '%02d' % iTRANS
#
			if 'JULES_Wetlands' in DATA_SOURCE:
				FILE_PLOT   = PLOT_DIR+FILE_DATA[jFILE][3]+'_'+str(int(WET_TOTAL))+'_timeseries'+FILE_CORR+'_TRANS'+sTRANS+'_'+SDATE+'.png'
				LEGEND      = [LEGEND_JULES]
			else:
				FILE_PLOT   = FILE_PLOT_PART+'_'+FSPECIES+'_'+SYEAR+'_timeseries_TRANS'+sTRANS+'.png'
				LEGEND      = [SOURCE_ALL[iFILE]]
#
			TITLE       = 'Time series of annual '+CSPECIES+' emissions for '+DOMAIN
			XLABEL      = 'Year'
			YLABEL      = 'Annual emissions (Tg '+CSPECIES+' per annum)'
			LEGEND_POS  = 0
			PLOT_CODE   = ['ko-']
#
			XMIN        = START_YEAR
			XMAX        = END_YEAR+1
			XTICKS      = int(XMIN)+XINC*arange(int((XMAX-XMIN)/XINC)+1)
#
			if iTRANS == 0:
				SCALE       = MAX_TREND_0
			else:
				SCALE       = 10.0
#
			YMIN        = 0
			YMAX        = SCALE*int(YDATA[~np.isnan(YDATA)].max()/SCALE+1)
			YINC        = YMAX/10.0
			NYTICKS     = int((YMAX-YMIN)/YINC)+1
			YTICKS      = YMIN + YINC*arange(NYTICKS)
#
			if PLOTS[1] == 'Y' and PLOTS[2] == 'Y':
				plot_functions.Plot_TimeSeries2(NDATASETS,XDATA,YDATA, \
					XTICKS[0],XTICKS[-1],XTICKS,XLABEL,YTICKS[0],YTICKS[-1],YTICKS,YLABEL, \
					TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
			if iOPT >= 2:
#
				if 'JULES_Wetlands' in DATA_SOURCE:
					FILE_PLOT   = PLOT_DIR+FILE_DATA[jFILE][3]+'_'+str(int(WET_TOTAL))+'_climatology'+FILE_CORR+'_TRANS'+sTRANS+'_'+SDATE+'.png'
				else:
					FILE_PLOT   = FILE_PLOT_PART+'_'+FSPECIES+'_'+SYEAR+'_clim_TRANS'+sTRANS+'.png'
#
				TITLE       = 'Climatology of '+CSPECIES+' emissions for '+DOMAIN
				XLABEL      = 'Month'
#
				XMIN        =  0.0
				XMAX        = 12.0
				XTICKS      = 0.5+arange(NMONTHS)
#
				if iTRANS == 0:
					SCALE       = MAX_CLIM_0
				else:
					SCALE       = 10.0
#
				YMIN        = 0
				YMAX        = SCALE*int(YCLIM[~np.isnan(YCLIM)].max()/SCALE+1)
				YINC        = YMAX/10.0
				NYTICKS     = int((YMAX-YMIN)/YINC)+1
				YTICKS      = YMIN + YINC*arange(NYTICKS)
#
				LEGEND_POS  = 0
#
				if PLOTS[2] == 'Y':
					plot_functions.Plot_Climatology(NDATASETS,XCLIM,YCLIM, \
						XMIN,XMAX,XTICKS,XLABEL,YTICKS[0],YTICKS[-1],YTICKS,YLABEL, \
						TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
# PLOTS[3]:
# Data for latitudinal plots
#
	if PLOTS[3] == 'Y':
#
		for iYEAR in range(NYEARS):
#
			SYEAR        = '%4d' % (iYEAR+START_YEAR)
			XZONAL[0,:]  = LAT_ZONAL
			YZONAL[0,:]  = EMISS_LAT[iFILE,iYEAR,:]*EMISS_FACTOR
			LAT_ZONAL_ALL.append(LAT_ZONAL)
#
# Zonal plot (Latitude)
#
			if iYEAR == 0:
				print('Minimum and maximum data values are :',EMISS_LAT.min()*EMISS_FACTOR,EMISS_LAT.max()*EMISS_FACTOR)
				print('Input maximum data value: ')
				EMISS_MAX      = float(input())
				EMISS_TICKS    = (EMISS_MAX/10.0)*arange(11)
#
			Emissions_Wetlands_Zonal.Emissions_Wetlands_Aggregate(FSPECIES,SYEAR,LAT_ZONAL,YZONAL[0,:])
#
			XLABEL      = 'Latitude'
			YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' yr$^{-1}$ '+('%.1f' % dLAT)+' Latitude Band$^{-1}$)'
			TITLE       = 'Latitudinal distribution of CH$_{4}$ emissions for '+SYEAR
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				TITLE      = 'Latitudinal distribution of '+CSPECIES+' emissions'
				LEGEND      = [SOURCE_ALL[iFILE]]
			else:
				TITLE      = 'Latitudinal distribution of '+CSPECIES+' emissions for '+SYEAR
				LEGEND      = [LEGEND_JULES]
#
			PLOT_CODE   = [ 'k-' ]
#
			if DATA_SOURCE == 'JULES_Wetlands':
				FILE_PLOT   = PLOT_DIR+FILE_DATA[jFILE][3]+'_'+str(int(WET_TOTAL))+'_zonal_'+SYEAR+FILE_CORR+'_'+SDATE+'.png'
			else:
				FILE_PLOT   = FILE_PLOT_PART+'_'+FSPECIES+'_'+SYEAR+'_zonal.png'
#
			plot_functions.Plot_Zonal(NDATA_ZONAL,YZONAL,XZONAL,0.0,EMISS_MAX,EMISS_TICKS,YLABEL, \
				LAT_START,LAT_END,LAT_TICKS,XLABEL,TITLE,LEGEND,FONTSIZES, \
				PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
# PLOTS[4]:
# Data for longitudinal plots
#
	if PLOTS[4] == 'Y':
#
		for iYEAR in range(NYEARS):
#
			SYEAR        = '%4d' % (iYEAR+START_YEAR)
			XMERID[0,:]  = LON_MERID
			YMERID[0,:]  = EMISS_LON[iFILE,iYEAR,:]*EMISS_FACTOR
			LON_MERID_ALL.append(LON_MERID)
#
# Zonal plot (Latitude)
#
			if iYEAR == 0:
				print('Minimum and maximum data values are :',EMISS_LON.min()*EMISS_FACTOR,EMISS_LON.max()*EMISS_FACTOR)
				print('Input maximum data value: ')
				EMISS_MAX      = float(input())
				EMISS_TICKS    = (EMISS_MAX/10.0)*arange(11)
#			
			XLABEL      = 'Longitude'
			YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' yr$^{-1}$ '+('%.1f' % dLON)+' Longitude Band$^{-1}$)'
#
			if DATA_SOURCE == 'TRANSCOM_Fung' or DATA_SOURCE == 'Wetland_Fung':
				TITLE      = 'Longitudinal distribution of '+CSPECIES+' emissions'
				LEGEND      = [SOURCE_ALL[iFILE]]
			else:
				TITLE      = 'Longitudinal distribution of '+CSPECIES+' emissions for '+SYEAR
				LEGEND      = [LEGEND_JULES]
#
			PLOT_CODE   = [ 'k-' ]
#
			if DATA_SOURCE == 'JULES_Wetlands':
				FILE_PLOT   = PLOT_DIR+FILE_DATA[jFILE][3]+'_'+str(int(WET_TOTAL))+'_meridional_'+SYEAR+FILE_CORR+'_'+SDATE+'.png'
			else:
				FILE_PLOT   = FILE_PLOT_PART+'_'+FSPECIES+'_'+SYEAR+'_meridional.png'
#
			plot_functions.Plot_Zonal(NDATA_MERID,XMERID,YMERID,LONG_PLOTS,LONG_PLOTE,LON_TICKS,XLABEL, \
				0.0,EMISS_MAX,EMISS_TICKS,YLABEL,TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT, \
				FILE_PLOT,DEBUG)
#
# Multi-file output
#
FILE_PART_MULTI='Test'
if NFILES >= 2 and PLOT_OPT == '2':
	print('Input partial file name for multi-file outputs: ')
	FILE_PART_MULTI= NETCDF_DIR+input()
#
# PLOTS[3]:
# Data for latitudinal plots
#
if PLOTS[3] == 'Y' and NFILES > 1:
#
	LEGEND           = []
	PLOT_CODE        = []
#
	for iYEAR in range(NYEARS):
#
		for iFILE in range(NFILES):
#
			SYEAR            = '%4d' % (iYEAR+START_YEAR)
			LAT_ZONAL        = LAT_ZONAL_ALL[iFILE]
			XZONAL[iFILE,:]  = LAT_ZONAL
			YZONAL[iFILE,:]  = EMISS_LAT[iFILE,iYEAR,:]*EMISS_FACTOR
#
			LEGEND.append(LEGENDS_MULTI[iFILE])
			PLOT_CODE.append(PLOT_CODES_ZL[iFILE])
#
# Zonal plot (Latitude)
#
			if iYEAR == 0 and iFILE == 0:
				print('Minimum and maximum data values are :',EMISS_LAT.min()*EMISS_FACTOR,EMISS_LAT.max()*EMISS_FACTOR)
				print('Input maximum data value: ')
				EMISS_MAX      = float(input())
				EMISS_TICKS    = (EMISS_MAX/10.0)*arange(11)
#
				Emissions_Wetlands_Zonal.Emissions_Wetlands_Aggregate(FSPECIES,SYEAR,LAT_ZONAL,YZONAL[iFILE,:])
#
		XLABEL      = 'Latitude'
		YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' yr$^{-1}$ '+('%.1f' % dLAT)+' Latitude Band$^{-1}$)'
		TITLE       = 'Latitudinal distribution of CH$_{4}$ emissions for '+SYEAR
		FILE_PLOT   = FILE_PART_MULTI+'_'+FSPECIES+'_'+SYEAR+'_zonal_multi.png'
#
		plot_functions.Plot_Zonal(NFILES,YZONAL,XZONAL,0.0,EMISS_MAX,EMISS_TICKS,YLABEL, \
			LAT_START,LAT_END,LAT_TICKS,XLABEL,TITLE,LEGEND,FONTSIZES, \
			PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
# PLOTS[4]:
# Data for longitudinal plots
#
if PLOTS[4] == 'Y' and NFILES > 1:
#
	LEGEND           = []
	PLOT_CODE        = []
#
	for iYEAR in range(NYEARS):
#
		for iFILE in range(NFILES):
#
			SYEAR        = '%4d' % (iYEAR+START_YEAR)
			LON_MERID    = LON_MERID_ALL[iFILE]
			XMERID[iFILE,:]  = LON_MERID
			YMERID[iFILE,:]  = EMISS_LON[iFILE,iYEAR,:]*EMISS_FACTOR
#
			LEGEND.append(LEGENDS_MULTI[iFILE])
			PLOT_CODE.append(PLOT_CODES_ZL[iFILE])
#
# Zonal plot (Latitude)
#
			if iYEAR == 0 and iFILE == 0:
				print('Minimum and maximum data values are :',EMISS_LON.min()*EMISS_FACTOR,EMISS_LON.max()*EMISS_FACTOR)
				print('Input maximum data value: ')
				EMISS_MAX      = float(input())
				EMISS_TICKS    = (EMISS_MAX/10.0)*arange(11)
#			
		XLABEL      = 'Longitude'
		YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' yr$^{-1}$ '+('%.1f' % dLON)+' Longitude Band$^{-1}$)'
		TITLE       = 'Longitudinal distribution of '+CSPECIES+' emissions for '+SYEAR
		FILE_PLOT   = FILE_PART_MULTI+'_'+FSPECIES+'_'+SYEAR+'_meridional_multi.png'
#
		plot_functions.Plot_Zonal(NFILES,XMERID,YMERID,LONG_PLOTS,LONG_PLOTE,LON_TICKS,XLABEL, \
			0.0,EMISS_MAX,EMISS_TICKS,YLABEL,TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT, \
			FILE_PLOT,DEBUG)
#
# DEBUG     = 'Y'
#
if PLOTS[6] == 'Y' and NFILES > 1 and iOPT >= 2:
#
	SDATE       = str(START_YEAR)+'_'+str(END_YEAR)
	FILE_PLOT   = FILE_PART_MULTI+'_'+FSPECIES+'_'+SDATE+'_timeseries_multi.png'
	PLOT_TITLE  = 'Time series of '+CSPECIES+' annual emissions'
	XLABEL      = 'Year'
	if iFIG_PAPER == 'V':
		YLABEL      = 'Emission (Tg '+CSPECIES+' per month)'
	else:
		YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' per month)'
	LEGEND_POS  = 0
#
	LEGEND      = []
	PLOT_CODE   = []
#
	for iFILE in range(NFILES):
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_GR[iFILE])
#
	XMIN        = START_YEAR
	XMAX        = END_YEAR+1
	XTICKS      = int(XMIN)+XINC*arange(int((XMAX-XMIN)/XINC)+1)
	XTICK_LABS  = XTICKS
#
	YMIN        = 0
	YMAX        = MAX_TREND_0*int(YDATA_MULTI[~np.isnan(YDATA_MULTI)].max()/MAX_TREND_0+1)
	YINC        = YMAX/10.0
	NYTICKS     = int((YMAX-YMIN)/YINC)+1
	YTICKS      = YMIN + YINC*arange(NYTICKS)
#
	NROWS       = NROWS_0
	NCOLUMNS    = NCOLUMNS_0
	NPLOTS      = NROWS*NCOLUMNS
	NDATASETS   = NFILES
	WIDTH       = 4.0*NCOLUMNS
	HEIGHT      = 4.0*NROWS
#
	XPLOT       = []
	YPLOT       = []
	PLOT_TRANS  = []
#
	if not iFIG_PAPER == 'N' and not iFIG_PAPER == 'V':
		TEMP        = YDATA_MULTI[:,1:NPLOTS,:]        
		TEMP1       = TEMP[~np.isnan(TEMP)].max()/YDATA_MULTI[~np.isnan(YDATA_MULTI)].max()
		INDEX       = (np.abs(PLOT_SCALE-TEMP1)).argmin()
		TEMP2       = PLOT_SCALE[INDEX]
#
		if TEMP1 > TEMP2:
			TEMP2       = PLOT_SCALE[INDEX+1]
#
		print YDATA_MULTI[:,1:NPLOTS,:].max(),YDATA_MULTI.max(),TEMP1,TEMP2
#
	for iPLOT in range(NPLOTS):
#
		PLOT_TRANS.append(TRANS_REGS[iPLOT])
#
		if iFIG_PAPER == 'N':
			if iPLOT == 0:
				YSCALE[0]     = 1.0
			else:	
				YSCALE[iPLOT] = TEMP2
#
		elif iFIG_PAPER == 'V':
			TEMP        = YDATA_MULTI[:,iPLOT,:]
			TEMP1       = TEMP[~np.isnan(TEMP)].max()
			INDEX       = (np.abs(PLOT_MAX-TEMP1)).argmin()
			TEMP2       = PLOT_MAX[INDEX]
#
			if TEMP1 > TEMP2:
				TEMP2       = PLOT_MAX[INDEX+1]
			if DEBUG == 'Y': print TEMP1,TEMP2
#
			YSCALE[iPLOT] = TEMP2/YMAX
#
		else:
			YSCALE[iPLOT] = FIG_PLOT_SCALES[iFIG_PAPER][1][iPLOT]
#
		for iFILE in range(NFILES):
#
			if DEBUG == 'Y':
				print XDATA_MULTI.min(),XDATA_MULTI.max(),YDATA_MULTI.min(),YDATA_MULTI.max()
			XPLOT.append(XDATA_MULTI[iFILE,iPLOT,:])
			YPLOT.append(YDATA_MULTI[iFILE,iPLOT,:])
#
	print YMAX,YSCALE
	plot_functions.Plot_General_MultiPlot_VarScale3(NROWS,NCOLUMNS,NDATASETS,XPLOT,YPLOT,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
		YMIN,YMAX,YTICKS,YTICKS,YLABEL,YSCALE,WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,PLOT_TRANS,PLOT_TEXT,LEGEND,LEGEND_POS, \
		PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
# Climatology
#
if PLOTS[7] == 'Y' and NFILES > 1 and iOPT >= 2:
#
	SDATE       = str(START_YEAR)+'_'+str(END_YEAR)
	FILE_PLOT   = FILE_PART_MULTI+'_'+FSPECIES+'_'+SDATE+'_climatology_multi.png'
	PLOT_TITLE  = 'Climatology of '+CSPECIES+' emissions'
	XLABEL      = 'Month'
	if iFIG_PAPER == 'V':
		YLABEL      = 'Emission (Tg '+CSPECIES+' per month)'
	else:
		YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' per month)'
	LEGEND_POS  = 0
	XTICK_LABS  = ['J','F','M','A','M','J','J','A','S','O','N','D' ]
#
	LEGEND      = []
	PLOT_CODE   = []
#
	for iFILE in range(NFILES):
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_GR[iFILE])
#
	XMIN        =  0.0
	XMAX        = 12.0
	XTICKS      = 0.5+arange(12)
#
	YMIN        = 0
	YMAX        = MAX_CLIM_0*int(YCLIM_MULTI[~np.isnan(YCLIM_MULTI)].max()/MAX_CLIM_0+1)
	YINC        = YMAX/10.0
	NYTICKS     = int((YMAX-YMIN)/YINC)+1
	YTICKS      = YMIN + YINC*arange(NYTICKS)
#
	NROWS       = NROWS_0
	NCOLUMNS    = NCOLUMNS_0
	NPLOTS      = NROWS*NCOLUMNS
	NDATASETS   = NFILES
	WIDTH       = 4.0*NCOLUMNS
	HEIGHT      = 4.0*NROWS
#
	XPLOT       = []
	YPLOT       = []
	PLOT_TRANS  = []
#
	if not iFIG_PAPER == 'N' and not iFIG_PAPER == 'V':
		TEMP        = YCLIM_MULTI[:,1:NPLOTS,:]        
		TEMP1       = TEMP[~np.isnan(TEMP)].max()/YCLIM_MULTI[~np.isnan(YCLIM_MULTI)].max()
		INDEX       = (np.abs(PLOT_SCALE-TEMP1)).argmin()
		TEMP2       = PLOT_SCALE[INDEX]
#
		if TEMP1 > TEMP2:
			TEMP2       = PLOT_SCALE[INDEX+1]
#
		print YCLIM_MULTI[:,1:NPLOTS,:].max(),YCLIM_MULTI.max(),TEMP1,TEMP2
#
	for iPLOT in range(NPLOTS):
#
		PLOT_TRANS.append(TRANS_REGS[iPLOT])
#
		if iFIG_PAPER == 'N':
			if iPLOT == 0:
				YSCALE[0]     = 1.0
			else:	
				YSCALE[iPLOT] = TEMP2
#
		elif iFIG_PAPER == 'V':
			TEMP        = YCLIM_MULTI[:,iPLOT,:]
			TEMP1       = TEMP[~np.isnan(TEMP)].max()
			INDEX       = (np.abs(PLOT_MAX-TEMP1)).argmin()
			TEMP2       = PLOT_MAX[INDEX]
#
			if TEMP1 > TEMP2:
				TEMP2       = PLOT_MAX[INDEX+1]
			if DEBUG == 'Y': print TEMP1,TEMP2
#
			YSCALE[iPLOT] = TEMP2/YMAX
#
		else:
			YSCALE[iPLOT] = FIG_PLOT_SCALES[iFIG_PAPER][2][iPLOT]
#
		for iFILE in range(NFILES):
#
			if DEBUG == 'Y':
				print XCLIM_MULTI.min(),XCLIM_MULTI.max(),YCLIM_MULTI.min(),YCLIM_MULTI.max()
			XPLOT.append(XCLIM_MULTI[iFILE,iPLOT,:])
			YPLOT.append(YCLIM_MULTI[iFILE,iPLOT,:])
#
	print YMAX,YSCALE
	plot_functions.Plot_General_MultiPlot_VarScale3(NROWS,NCOLUMNS,NDATASETS,XPLOT,YPLOT,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
		YMIN,YMAX,YTICKS,YTICKS,YLABEL,YSCALE,WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,PLOT_TRANS,PLOT_TEXT,LEGEND,LEGEND_POS, \
		PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
if (PLOTS[9] == 'Y' or PLOTS[12] == 'Y') and NFILES > 1:
#
	FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME)
	for i in range(len(CLEVELS)):
		CLEVELS[i] = CLEVELS[i]/10.0
#
	iTIME        = 0
#
	if MAP_FRACT == 'Y':
		HEIGHT       = HEIGHT0*2
		WIDTH        = WIDTH0*NFILES
		NROWS        = 2
		NCOLUMNS     = NFILES
		iSTART       = 2*iFILE*(NYEARS*NMONTHS)
	else:
		if NFILES == 1:
			HEIGHT       = HEIGHT0
			WIDTH        = WIDTH0
			NROWS        = 1
			NCOLUMNS     = 1
		elif NFILES == 2 and PLOTS[12] == 'Y':
			HEIGHT       = HEIGHT0*(NFILES+1)
			WIDTH        = WIDTH0
			NROWS        = NFILES+1
			NCOLUMNS     = 1
			print NFILES,NROWS,NCOLUMNS
		elif NFILES == 2 or NFILES == 3:
			HEIGHT       = HEIGHT0*NFILES
			WIDTH        = WIDTH0
			NROWS        = NFILES
			NCOLUMNS     = 1
		elif NFILES >= 4:
			HEIGHT       = HEIGHT0*2
			WIDTH        = WIDTH0*NFILES/2
			NROWS        = 2
			NCOLUMNS     = NFILES/2
		iSTART       = iFILE*(NYEARS*NMONTHS)
#
	if PLOTS[12] == 'Y':
		if SRESOL_KMZ[0] == SRESOL_KMZ[1]:
			INDEX        = -1
		else:
			INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT= \
				data_regrid_new.data_regrid_factor(DEBUG,NFILES,SRESOL_KMZ,LONG_START,LAT_START,LONG_END,LAT_END)
			NUM_IN_MAP   = np.ones((NLAT_IN,NLONG_IN))
			print INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
#
	for iYEAR in range(NYEARS):
#
		for iMONTH in range(NMONTHS):
#
			SDATE        = ('%4d%02d' % (iYEAR+START_YEAR,iMONTH+1))
			CLEVELS_ALL  = []
			COLOURS_ALL  = []
			MAP_TYPES    = []
			DATA_MAP_ALL = []
			PLOT_LABELS  = []
			SUB_TITLES   = []
			SET_UNDER_ALL= []
			SET_OVER_ALL = []
#
			for iFILE in range(NFILES):
#
				if MAP_FRACT == 'Y':
					iSTART       = 2*iFILE*NYEARS*NMONTHS
					COLOURS_ALL.append(COLOURS_C)
					CLEVELS_ALL.append(CLEVELS_C)
					MAP_TYPES.append(MAP_TYPE)
					PLOT_LABELS.append('Fraction of grid-square')
					SUB_TITLES.append(SOURCE_ALL[iFILE]+': Grid-square wetland fraction')
					SET_UNDER_ALL.append(SET_UNDER)
					SET_OVER_ALL.append(SET_OVER)
#
					DATA_MAP_ALL.append(DATA_MONTH_ALL[iSTART+2*iTIME])
					DATA_MAP     = DATA_MONTH_ALL[iSTART+2*iTIME+1]*EMISS_FACTORS[iFILE]
					DATA_MAP[DATA_MAP <=0] = MISS_DATA
					DATA_MAP_ALL.append(DATA_MAP)
#
					PLOT_TITLE   = 'Wetland fraction and emissions of '+CSPECIES+' for '+SDATE
					FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+'_wetlands'+SDOMAIN+SDATE+'.png'
				else:
					iSTART       = iFILE*NYEARS*NMONTHS
					PLOT_TITLE   = 'Emissions of '+CSPECIES+' for '+SDATE
					FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+SDOMAIN+SDATE+'.png'
					DATA_MAP     = DATA_MONTH_ALL[iSTART+iTIME]
					DATA_MAP[DATA_MAP <=0] = MISS_DATA
					DATA_MAP_ALL.append(DATA_MAP)
#
					if iFILE == NFILES-1 and PLOTS[12] == 'Y':
#
						DATA_0,DATA_1,DATA_D=data_regrid_new.data_map_difference \
							(DEBUG,INDEX,DATA_MAP_ALL,NUM_IN_MAP,MIN_EMISS,MISS_DATA, \
							NLONG,NLAT,NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
#						oLONG       = 239
#						oLAT        = 149
#						print TEMP0[oLAT,oLONG],TEMP1[oLAT,oLONG],TEMPD[oLAT,oLONG]
#						oLONG       = 299
#						oLAT        = 149
#						print TEMP0[oLAT,oLONG],TEMP1[oLAT,oLONG],TEMPD[oLAT,oLONG]
#
						DATA_MAP_ALL = []
						DATA_MAP_ALL.append(DATA_0)
						DATA_MAP_ALL.append(DATA_1)
						DATA_MAP_ALL.append(DATA_D)
#
				CLEVELS[0] = MIN_EMISS
				CLEVELS_ALL.append(CLEVELS)
				COLOURS_ALL.append(COLOURS)
				MAP_TYPES.append(MAP_TYPE)
				PLOT_LABELS.append('Emissions of '+CSPECIES+' ('+UNITS2+')')
				SUB_TITLES.append(SOURCE_ALL[iFILE]+': Monthly emissions')
				SET_UNDER_ALL.append(SET_UNDER)
				SET_OVER_ALL.append(SET_OVER)
#
				if iFILE == NFILES-1 and PLOTS[12] == 'Y':
					CLEVELS_ALL.append(np.array(CLEVELS_D)*max(CLEVELS))
					COLOURS_ALL.append(COLOURS_D)
					MAP_TYPES.append(MAP_TYPE)
					PLOT_LABELS.append('Emissions of '+CSPECIES+' ('+UNITS2+')')
					SUB_TITLES.append('Difference: Monthly emissions')
					SET_UNDER_ALL.append(SET_UNDER_DIFF)
					SET_OVER_ALL.append(SET_OVER)
#
			print FILE_PLOT
#
			plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_MAP_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
				LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
				RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
			iTIME        += 1
#
if PLOTS[10] == 'Y' and NFILES > 1 and iOPT >= 2:
#
	SDATE        = str(START_YEAR)+'_'+str(END_YEAR)
	TRANS_ANOM   = np.zeros(NMONTHS*NYEARS)
	TIME_ANOM    = np.zeros(NMONTHS*NYEARS)
#
	for iFILE in range(NFILES):
		for iTRAN in range(NTRANS):
#
			TRANS_ANOM[:] = float('nan')
#
			for iMONTH in range(NMONTHS*NYEARS):
				if iMONTH > 5 and iMONTH <= NMONTHS*NYEARS-6:
					TRANS_EMISS_RM[iFILE,iTRAN,iMONTH] = TRANS_EMISS_MN[iFILE,iTRAN,iMONTH-6:iMONTH+6].sum()
				else:
					TRANS_EMISS_RM[iFILE,iTRAN,iMONTH] = float('nan')
#
			TRANS_LOCAL  = TRANS_EMISS_RM[iFILE,iTRAN,:]
			TRANS_MEAN   = TRANS_LOCAL[~np.isnan(TRANS_LOCAL)].mean()
			TRANS_ANOM[~np.isnan(TRANS_LOCAL)] = TRANS_LOCAL[~np.isnan(TRANS_LOCAL)]-TRANS_MEAN
#
			TRANS_EMISS_AN[iFILE,iTRAN,:] = TRANS_ANOM
#
			if iTRAN == 0:
				print iFILE,iTRAN
				print TRANS_EMISS_MN[iFILE,0,:]
				print TRANS_EMISS_RM[iFILE,0,:]
				print TRANS_EMISS_AN[iFILE,0,:]
#
				print TRANS_MEAN
#
	for iMONTH in range(NMONTHS*NYEARS):
		TIME_ANOM[iMONTH] = START_YEAR+(iMONTH+0.5)/12
#
	FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+'_'+SDATE+'_timeseries_anom_multi.png'
#
	PLOT_TITLE  = 'Time series of anomalies in deseasonlised '+CSPECIES+' emissions'
	XLABEL      = 'Year'
	YLABEL      = CSPECIES+' Emission (Tg '+CSPECIES+' per annum)'
	LEGEND_POS  = 0
#
	LEGEND      = []
	PLOT_CODE   = []
#
	for iFILE in range(NFILES):
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_AN[iFILE])
#
	XMIN        = START_YEAR
	XMAX        = END_YEAR+1
	XTICKS      = int(XMIN)+XINC*arange(int((XMAX-XMIN)/XINC)+1)
	XTICK_LABS  = XTICKS
#
	MAX_PLOT    = TRANS_EMISS_AN[~np.isnan(TRANS_EMISS_AN)].max()
	YMAX        = MAX_ANOM_0*int(MAX_PLOT/MAX_ANOM_0+1)
	YMIN        = -YMAX
	YINC        = YMAX/5.0
	NYTICKS     = int((YMAX-YMIN)/YINC)+1
	YTICKS      = YMIN + YINC*arange(NYTICKS)
#
	NROWS       = NROWS_0
	NCOLUMNS    = NCOLUMNS_0
	NPLOTS      = NROWS*NCOLUMNS
	NDATASETS   = NFILES
	WIDTH       = 4.0*NCOLUMNS
	HEIGHT      = 4.0*NROWS
#
	XPLOT       = []
	YPLOT       = []
	PLOT_TRANS  = []
#
	for iPLOT in range(NPLOTS):
#
		PLOT_TRANS.append(TRANS_REGS[iPLOT])
#
		if iFIG_PAPER == 'N':
			if iPLOT == 0:
				YSCALE[0]     = 1.0
			else:
				TEMP        = TRANS_EMISS_AN[:,iPLOT,:]
				TEMP1       = TEMP[~np.isnan(TEMP)].max()/MAX_PLOT
				INDEX       = (np.abs(PLOT_SCALE-TEMP1)).argmin()
				TEMP2       = PLOT_SCALE[INDEX]
#
				if TEMP1 > TEMP2:
					TEMP2       = PLOT_SCALE[INDEX-1]
#
				print TEMP[~np.isnan(TEMP)].max(),MAX_PLOT,TEMP1,TEMP2
#
				YSCALE[iPLOT] = TEMP2
#
		elif iFIG_PAPER == 'V':
			TEMP        = TRANS_EMISS_AN[:,iPLOT,:]
			TEMP1       = TEMP[~np.isnan(TEMP)].max()
			INDEX       = (np.abs(PLOT_MAX-TEMP1)).argmin()
			TEMP2       = PLOT_MAX[INDEX]
#
			if TEMP1 > TEMP2:
				TEMP2       = PLOT_MAX[INDEX+1]
			if DEBUG == 'Y': print TEMP1,TEMP2
#
			YSCALE[iPLOT] = TEMP2/YMAX
#
		else:
			YSCALE[iPLOT] = FIG_PLOT_SCALES[iFIG_PAPER][3][iPLOT]
#
		for iFILE in range(NFILES):
#
			if DEBUG == 'Y':
				print XDATA_MULTI.min(),XDATA_MULTI.max(),TRANS_EMISS_AN[~np.isnan(TRANS_EMISS_AN)].min(),TRANS_EMISS_AN[~np.isnan(TRANS_EMISS_AN)].max()
			XPLOT.append(TIME_ANOM)
			YPLOT.append(TRANS_EMISS_AN[iFILE,iPLOT,:])
#
	print YMAX,YSCALE
	plot_functions.Plot_General_MultiPlot_VarScale3(NROWS,NCOLUMNS,NDATASETS,XPLOT,YPLOT,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
		YMIN,YMAX,YTICKS,YTICKS,YLABEL,YSCALE,WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,PLOT_TRANS,PLOT_TEXT,LEGEND,LEGEND_POS, \
		PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
if PLOTS[11] == 'Y' and iOPT >= 2:
#
        import write_GoogleEarth
#
	sANIM          = KMZ_CODES[0]
	sLEGEND        = KMZ_CODES[1] 
	sEXPAND        = KMZ_CODES[2]
#
	for iFILE in range(NFILES):
#
		SDATES         = []
		MTIMES         = NDIMS_KMZ[iFILE][0]
		DATA_MAP_ALL   = np.zeros(NDIMS_KMZ[iFILE])
#
		KMZ_FILE       = NAMES_KMZ[iFILE][:-9].replace(' ','_').replace('+','_')
		KMZ_FILE_END   = '%4d_%4d' % (START_YEAR,END_YEAR)
                PROGNAME       = NAMES_KMZ[iFILE].replace('+',' ')
		LEGEND         = UNITS_KMZ[iFILE]
		MISS_DATA      = MISS_DATA_KMZ[iFILE]
#
		LAT_KMZ        = LATS_KMZ[iFILE]
		LONG_KMZ       = LONGS_KMZ[iFILE]
#
		if sEXPAND == 'Y':
#
			COLOURS        = [ '#9400d3', '#8a2be2', '#0000cd', '#0000ff', '#4169e1', '#6495ed', \
			                   '#4682b4', '#228B22', '#00ff00', '#98fb98', '#adff2f', '#ffff00', \
			                   '#ffeb00', '#ffd700', '#ffa500', '#ff8c00', '#ff0000', '#ff1493' ]
#
			CLEVELS        = [  0.1,    0.2,    0.5,    1.0,     2.0,     3.0,     4.0,     5.0,     6.0,     7.0,     8.0, \
					    9.0,   10.0,   20.0,   50.0,   100.0,   200.0,   500.0  ]
#
		else:
			COLOURS        = COLOURS_KMZ[iFILE]
			CLEVELS        = CLEVELS_KMZ[iFILE]
#
# Create images
#
		for iTIMES in range(MTIMES):
			DATA_MAP_ALL[iTIMES,:,:] = EMISS_KMZ[iTIMES+iFILE*MTIMES]
			SDATES.append(str(SDATE_TIMES_KMZ[iFILE][iTIMES][0:4])+str(SDATE_TIMES_KMZ[iFILE][iTIMES][5:7]))
			print SDATE_TIMES_KMZ[iFILE][iTIMES],EDATE_TIMES_KMZ[iFILE][iTIMES],SDATES[iTIMES]
#
		write_GoogleEarth.control(KMZ_DIR,KMZ_FILE,KMZ_FILE_END,DATA_MAP_ALL, \
			LONG_KMZ,LAT_KMZ,LONG_MAP,LAT_MAP, \
			MAP_TYPE,CLEVELS,COLOURS,PROGNAME,LEGEND,FONTSIZES,SET_UNDER,SET_OVER, \
			SDATES,SDATE_TIMES_KMZ[iFILE],EDATE_TIMES_KMZ[iFILE], \
			WIDTH0,sANIM,sLEGEND,MISS_DATA,DEBUG)
#
# For JULES wetlands, write out wetland fraction
#
		if 'JULES' in KMZ_FILE:
#
			KMZ_FILE       = KMZ_FILE.replace('Emissions','Fraction')
			PROGNAME       = NAMES_KMZ[iFILE].replace('Emissions','Fraction').replace('+',' ')
			LEGEND         = 'Wetland extent (as fraction of grid square)'
#
			DATA_MAP_ALL   = WET_FRACT_KMZ[iFILE][:,:,:]
			print WET_FRACT_KMZ[iFILE].shape
#
			if sEXPAND == 'Y':
#
				COLOURS        = [ '#9400d3', '#8a2be2', '#0000cd', '#0000ff', '#4169e1', '#6495ed', \
				                   '#4682b4', '#228B22', '#00ff00', '#98fb98', '#adff2f', '#ffff00', \
				                   '#ffeb00', '#ffd700', '#ffa500', '#ff8c00', '#ff0000', '#ff1493' ]
#
				CLEVELS        = [ 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09, \
						   0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]
			else:
				COLOURS        = COLOURS_C
				CLEVELS        = CLEVELS_C				
#
		write_GoogleEarth.control(KMZ_DIR,KMZ_FILE,KMZ_FILE_END,DATA_MAP_ALL, \
				LONG_KMZ,LAT_KMZ,LONG_MAP,LAT_MAP, \
				MAP_TYPE,CLEVELS,COLOURS,PROGNAME,LEGEND,FONTSIZES,SET_UNDER,SET_OVER, \
				SDATES,SDATE_TIMES_KMZ[iFILE],EDATE_TIMES_KMZ[iFILE], \
				WIDTH0,sANIM,sLEGEND,MISS_DATA,DEBUG)
#
# Plot annual maps and difference map for each year
#
if PLOTS[0] == 'Y' and NFILES == 2:
#
	FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME)
	NUM_IN_MAP   = np.ones((NLAT_IN,NLONG_IN))
#
	HEIGHT       = HEIGHT0*(NFILES+1)
	WIDTH        = WIDTH0
	NROWS        = NFILES+1
	NCOLUMNS     = 1
	print NFILES,NROWS,NCOLUMNS
#
	if SRESOL_KMZ[0] == SRESOL_KMZ[1]:
		INDEX        = -1
	else:
		INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT= \
			data_regrid_new.data_regrid_factor(DEBUG,NFILES,SRESOL_KMZ,LONG_START,LAT_START,LONG_END,LAT_END)
		NUM_IN_MAP   = np.ones((NLAT_IN,NLONG_IN))
		print INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
#
	for iYEAR in range(NYEARS):
#
		SDATE        = ('%4d' % (iYEAR+START_YEAR))
		CLEVELS_ALL  = []
		COLOURS_ALL  = []
		MAP_TYPES    = []
		DATA_MAP_ALL = []
		PLOT_LABELS  = []
		SUB_TITLES   = []
		SET_UNDER_ALL= []
		SET_OVER_ALL = []
#
		for iFILE in range(NFILES):
#
			iSTART       = iFILE*NYEARS
			PLOT_TITLE   = 'Emissions of '+CSPECIES+' for '+SDATE
			FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+SDOMAIN+SDATE+'_difference_annual.png'
			DATA_MAP     = DATA_ANN_ALL[iSTART+iYEAR]
			DATA_MAP[DATA_MAP <=0] = MISS_DATA
			DATA_MAP_ALL.append(DATA_MAP)
#
			if iFILE == NFILES-1:
#
				DATA_0,DATA_1,DATA_D=data_regrid_new.data_map_difference \
					(DEBUG,INDEX,DATA_MAP_ALL,NUM_IN_MAP,MIN_EMISS,MISS_DATA, \
					NLONG,NLAT,NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
# Reset values for plotting
#
				DATA_0[DATA_0 <= MIN_EMISS] = float('nan')
				DATA_1[DATA_1 <= MIN_EMISS] = float('nan')
				DATA_D[(DATA_D >= -MIN_EMISS) & (DATA_D <= MIN_EMISS)] = float('nan')
				print DATA_0.min(),DATA_0.max(),DATA_1.min(),DATA_1.max(),DATA_D.min(),DATA_D.max()
#
			CLEVELS[0] = MIN_EMISS
			CLEVELS_ALL.append(CLEVELS)
			COLOURS_ALL.append(COLOURS)
			MAP_TYPES.append(MAP_TYPE)
			PLOT_LABELS.append('Emissions of '+CSPECIES+' ('+UNITS2+')')
			SUB_TITLES.append(SOURCE_ALL[iFILE])
			SET_UNDER_ALL.append(SET_UNDER)
			SET_OVER_ALL.append(SET_OVER)
#
			if iFILE == NFILES-1:
				CLEVELS_ALL.append(np.array(CLEVELS_D)*max(CLEVELS))
				COLOURS_ALL.append(COLOURS_D)
				MAP_TYPES.append(MAP_TYPE)
				PLOT_LABELS.append('Emissions of '+CSPECIES+' ('+UNITS2+')')
				SUB_TITLES.append('Difference:'+SOURCE_ALL[0]+'-'+SOURCE_ALL[1])
				SET_UNDER_ALL.append(SET_UNDER_DIFF)
				SET_OVER_ALL.append(SET_OVER)
#
		DATA_MAP_ALL = []
		DATA_MAP_ALL.append(DATA_0)
		DATA_MAP_ALL.append(DATA_1)
		DATA_MAP_ALL.append(DATA_D)
#
		print FILE_PLOT
#
		plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_MAP_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
			LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
			RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
# Plot monthly difference maps for each year
#
if PLOTS[0] == 'Y' and NFILES == 2:
#
	if DATA_NAME in DATA_NAMES:
		DATA_NAME_R = 'ch4_surf_emiss_base'
		FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME_R)
	else:
		FSPECIES,CSPECIES,CLEVELS,UNITS,COLOURS,RMM=data_info.data_map_info(DATA_NAME)
#
	iTIME        = 0
#
	NROWS        = 6
	NCOLUMNS     = 2
	HEIGHT       = HEIGHT0*NROWS
	WIDTH        = WIDTH0*NCOLUMNS
	iSTART       = iFILE*(NYEARS*NMONTHS)
#
	if SRESOL_KMZ[0] == SRESOL_KMZ[1]:
		INDEX        = -1
	else:
		INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT= \
			data_regrid_new.data_regrid_factor(DEBUG,NFILES,SRESOL_KMZ,LONG_START,LAT_START,LONG_END,LAT_END)
		NUM_IN_MAP   = np.ones((NLAT_IN,NLONG_IN))
		print INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
#
	for iYEAR in range(NYEARS):
#
		SDATE        = ('%4d' % (iYEAR+START_YEAR))
		CLEVELS_ALL  = []
		COLOURS_ALL  = []
		MAP_TYPES    = []
		DATA_MAP_ALL = []
		PLOT_LABELS  = []
		SUB_TITLES   = SMONTHS
		SET_UNDER_ALL= []
		SET_OVER_ALL = []
#
		for iMONTH in range(NMONTHS):
#
			DATA_MAP_DIFF = []
#
			for iFILE in range(NFILES):
#
				iSTART       = iFILE*NYEARS*NMONTHS
				PLOT_TITLE   = 'Difference in emissions of '+CSPECIES+' for '+SDATE
				FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+SDOMAIN+SDATE+'_difference_monthly.png'
				DATA_MAP     = DATA_MONTH_ALL2[iSTART+iTIME]
				DATA_MAP[DATA_MAP <=0] = MISS_DATA
				DATA_MAP_DIFF.append(DATA_MAP)
#
				if iFILE == NFILES-1:
#
					DATA_0,DATA_1,DATA_D=data_regrid_new.data_map_difference \
						(DEBUG,INDEX,DATA_MAP_DIFF,NUM_IN_MAP,MIN_EMISS,MISS_DATA, \
						NLONG,NLAT,NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
					DATA_MAP_ALL.append(DATA_D)
#
					CLEVELS_ALL.append(np.array(CLEVELS_D)*max(CLEVELS))
					COLOURS_ALL.append(COLOURS_D)
					MAP_TYPES.append(MAP_TYPE)
					PLOT_LABELS.append('Emissions of '+CSPECIES+' ('+UNITS2+')')
					SET_UNDER_ALL.append(SET_UNDER_DIFF)
					SET_OVER_ALL.append(SET_OVER)
#
			iTIME        += 1
#
		print FILE_PLOT
#
		plot_map.plot_map3_multi(NROWS,NCOLUMNS,DATA_MAP_ALL,LONG_DOMS,LONG_DOME,LONG,LONG_MAP, \
			LAT_DOMS,LAT_DOME,LAT,LAT_MAP,CLEVELS_ALL,COLOURS_ALL,MAP_TYPES,WIDTH,HEIGHT,ASPECT, \
			RESOLUTION,PLOT_TITLE,SUB_TITLES,PLOT_LABELS,FILE_PLOT,PLOT_MAP,FONTSIZES,SET_UNDER_ALL,SET_OVER_ALL,DEBUG)
#
if PLOTS[14] == 'Y':
#
# Annual cycle at sites
#
	NROWS       = 3
	NCOLUMNS    = 2
	HEIGHT      = 4.0*NROWS
	WIDTH       = 4.0*NCOLUMNS
	XLABEL      = 'Month'
	YLABEL      = CSPECIES+' Emission (mg '+CSPECIES+' m$^{-2}$ d$^{-1}$)'
	LEGEND_POS  = 0
	XTICK_LABS  = ['J','F','M','A','M','J','J','A','S','O','N','D' ]
#
	XMIN        =  0.0
	XMAX        = 12.0
	XTICKS      = 0.5+arange(12)
#
	EMISS_DATA  = np.zeros((NMONTHS*NYEARS,NSITES,NFILES))
#
	LEGEND      = []
	PLOT_CODE   = []
#
	wTIME        = (START_YEAR-START_DATA)*12
#
	for iFILE in range(NFILES):
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_GR[iFILE])
#
		for iYEAR in range(NYEARS):
#
# Monthly site-specific fluxes (period means)
#
			iSTART       = iFILE*NYEARS+iYEAR
#
			for iSITE in range(NSITES):
#
				SITE_NAME    = SITE_DATA[iSITE][0]
				SITE_LONG    = float(SITE_DATA[iSITE][1])
				SITE_LAT     = float(SITE_DATA[iSITE][2])
#
				sLONG        = RESOL_LONG_ALL[iFILE*NSITES+iSITE]
				sLAT         = RESOL_LAT_ALL[iFILE*NSITES+iSITE]
#
# Extract emissions for location of site and convert from Tg(CH4) yr-1 to mg(CH4) m-2 day-1
#
				FACTOR       = []
#
				for iMONTH in range(NMONTHS):
#
					WET_FRACT_LOC = WET_FRACT_KMZ[iFILE][iYEAR*12+iMONTH,sLAT,sLONG]
					FACTOR.append(1.00E+15 \
						/(AREAS_ALL[iFILE][sLAT,sLONG]*DAYS_MONTH[iMONTH]))
#
					if np.isnan(WET_FRACT_LOC) or WET_FRACT_LOC <= 0:
						EMISS_DATA[iYEAR*12+iMONTH,iSITE,iFILE] = 0.0
					else:
						EMISS_DATA[iYEAR*12+iMONTH,iSITE,iFILE] = \
							EMISS_SITE[iSTART][iMONTH,sLAT,sLONG]*FACTOR[iMONTH]/WET_FRACT_LOC
#
	YMIN        = 0
	YMAX        = MAX_SITE_0*int(EMISS_DATA.max()/MAX_SITE_0+1)
	YINC        = YMAX/10.0
	NYTICKS     = int((YMAX-YMIN)/YINC)+1
	YTICKS      = YMIN + YINC*arange(NYTICKS)
	print EMISS_DATA.max(),YMIN,YMAX
#
	XPLOT       = []
	YPLOT       = []
	PLOT_SITES  = []
#
	TEMP1       = EMISS_DATA.max()/YMAX
	INDEX       = (np.abs(PLOT_SCALE-TEMP1)).argmin()
	TEMP2       = PLOT_SCALE[INDEX]
#
	if TEMP1 > TEMP2:
		TEMP2       = PLOT_SCALE[INDEX-1]
#
	print EMISS_DATA.max(),EMISS_DATA.max(),TEMP1,TEMP2
#
	for iYEAR in range(NYEARS):
#
# Monthly site-specific fluxes (period means)
#
		SDATE        = ('%4d' % (iYEAR+START_YEAR))
		PLOT_TITLE   = 'Annual cycle of '+CSPECIES+' emissions for '+SDATE
		if NFILES == 1:
			FILE_PLOT    = FILE_PLOT_PART+'_'+FSPECIES+'_'+SDATE+'_site_annual_cycle.png'	
		else:
			FILE_PLOT    = FILE_PART_MULTI+'_'+FSPECIES+'_'+SDATE+'_site_annual_cycle_multi.png'
		PLOT_SUBTITLES = []
#
		for iSITE in range(NSITES):
#
			SITE_NAME     = SITE_DATA[iSITE][0]
			PLOT_SUBTITLES.append(SITE_NAME)
			YSCALE[iSITE] = TEMP2
#
			SITE_LONG    = float(SITE_DATA[iSITE][1])
			SITE_LAT     = float(SITE_DATA[iSITE][2])
			sLONG        = RESOL_LONG_ALL[iFILE*NSITES+iSITE]
			sLAT         = RESOL_LAT_ALL[iFILE*NSITES+iSITE]
			print iSITE,iFILE,iYEAR,SITE_NAME,SITE_LONG,SITE_LAT,sLONG,sLAT
#
			for iFILE in range(NFILES):
#
				if DEBUG == 'Y':
					print EMISS_DATA.min(),EMISS_DATA.max()
				XPLOT.append(XTICKS)
				YPLOT.append(EMISS_DATA[iYEAR*NMONTHS:(iYEAR+1)*NMONTHS,iSITE,iFILE])
#
				iSTART       = iFILE*NYEARS+iYEAR
				print iSTART
				print EMISS_DATA[iYEAR*NMONTHS:(iYEAR+1)*NMONTHS,iSITE,iFILE]
				print WET_FRACT_KMZ[iFILE][iYEAR*12:iYEAR*12+NMONTHS,sLAT,sLONG]
				print EMISS_SITE[iSTART][:,sLAT,sLONG]
				print FACTOR[:]
#
		plot_functions.Plot_General_MultiPlot_VarScale3(NROWS,NCOLUMNS,NFILES,XPLOT,YPLOT,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
			YMIN,YMAX,YTICKS,YTICKS,YLABEL,YSCALE,WIDTH,HEIGHT,FONTSIZES,PLOT_TITLE,PLOT_SUBTITLES,PLOT_TEXT,LEGEND,LEGEND_POS, \
			PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
# PLOTS[16] - Extract relative emission contribution
#
if PLOTS[16] == 'Y':
#
	import data_Surface_CH4
	NUM_SITES,SITE_DATA=data_Surface_CH4.data_Surface_CH4_setup(DEBUG)
#
	NLATS        = np.zeros(NFILES)
	NLONGS       = np.zeros(NFILES)
#
	CONTINUE     = 1
	while CONTINUE == 1:
#
		for iSITE in range(0,NUM_SITES,4):
			TEXT         = '%2d %s %2d %s %2d %s %2d %s' % ( \
				iSITE,SITE_DATA[iSITE][0],iSITE+1,SITE_DATA[iSITE+1][0], \
				iSITE+2,SITE_DATA[iSITE+2][0],iSITE+3,SITE_DATA[iSITE+3][0] \
				)
			print TEXT
		print 'Input site number: '
		iSITE        = int(input())
#
		xLONG        = float(SITE_DATA[iSITE][3])
		xLAT         = float(SITE_DATA[iSITE][2])
#
		print 'Input size of lat/long box centred on site for emission extraction: '
		xDOMAIN_SIZE = float(input())
		print xLAT,xLONG
#
# Identify grid squares to extract
#
		INDICES_ALL  = []
		TEXT_LABEL   = 'Year      '
#
		for iFILE in range(NFILES):
#
			TEXT_LABEL   = TEXT_LABEL+(' %20s' % SOURCE_FILES[iFILE][0:20])
			LONG         = LONGS_ALL[iFILE]
			LAT          = LATS_ALL[iFILE]
#
			NLATS[iFILE] = len(LAT)
			NLONGS[iFILE]= len(LONG)
#
# Emissions will be from 180W to 180 E
#
			LONG_IDX     = (np.where((LONG >= xLONG-xDOMAIN_SIZE) & (LONG <= xLONG+xDOMAIN_SIZE)))[0].tolist()
			LAT_IDX      = (np.where((LAT  >= xLAT -xDOMAIN_SIZE) & (LAT  <= xLAT +xDOMAIN_SIZE)))[0].tolist()
			INDICES_ALL.append([LAT_IDX,LONG_IDX])
#
			if len(LAT_IDX) == 0 or len(LONG_IDX) == 0:
				print iFILE,xLAT,xLONG,xDOMAIN_SIZE
				print 'Latitude'
				print LAT
				print LAT_IDX
				print 'Longitude:'
				print LONG
				print LONG_IDX
#
			print iFILE,len(LAT_IDX),len(LONG_IDX)
			print iFILE,LAT_IDX[0],LAT_IDX[-1],LONG_IDX[0],LONG_IDX[-1], \
				LAT[LAT_IDX[0]],LAT[LAT_IDX[-1]],LONG[LONG_IDX[0]],LONG[LONG_IDX[-1]]
#
# Extract emissions from each dataset
#
		EMIS_EXT     = np.zeros(NFILES)
		TEXT_LABEL   = TEXT_LABEL+'                Total'
		print
		print TEXT_LABEL
#
		for iYEAR in range(NYEARS):
#
			iSTART2      = iYEAR*NMONTHS
			iEND2        = iSTART2+NMONTHS
			SYEAR        = str(START_YEAR+iYEAR)
#
			for iFILE in range(NFILES):
#
				EMISS_TEMP   = np.zeros((NLATS[iFILE],NLONGS[iFILE]))
# Retrieve dataset
#
				if NTIMES_ALL[iFILE] == NYEARS:
					if DEBUG == 'Y': print iFILE,iYEAR,MON_EMISS_ALL[iFILE].shape
					EMISS_TEMP      = MON_EMISS_ALL[iFILE][iYEAR,:,:]
				else:
					for iMONTH in range(NMONTHS):
						EMISS_TEMP      = EMISS_TEMP+MON_EMISS_ALL[iFILE][iSTART2+iMONTH,:,:]
#
# Extract grid squares
				INDICES         = INDICES_ALL[iFILE]
#				print INDICES[0]
#				print INDICES[1]
				EMIS_EXT[iFILE] = EMISS_TEMP[np.ix_(INDICES[0],INDICES[1])].sum()
#
			TEXT_EMISS   = '%-10d' % (iYEAR+START_YEAR)
			for iFILE in range(NFILES):
				TEXT_EMISS   = TEXT_EMISS+ \
					(' %12.4E (%5.1f)' % (EMIS_EXT[iFILE],100.0*EMIS_EXT[iFILE]/EMIS_EXT.sum()))
			TEXT_EMISS   = TEXT_EMISS+ \
				(' %12.4E (100.0)' % (EMIS_EXT.sum()))
			print TEXT_EMISS
#
		print
		print 'Continue (1/0): '
		CONTINUE     = int(input())
#
# End of Program
