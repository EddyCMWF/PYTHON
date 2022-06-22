#!/usr/bin/python
#
# Python code to derive wetland inundation statistics and plots
#
# ./Wetlands_plot.py N Y Linux m45 20130319 1993 2007 0 1 NYYYYY Y
# ./Wetlands_plot.py N N Linux m45 20130319 1993 2007 0 1 NYYYYY N 0 4 1
# ./Wetlands_plot.py N N JASMIN WFDEI 20130319 1993 2014 0 1 NYYYYY N 0 4 1
#
# Garry Hayman
# Centre for Ecology and Hydrology
# September 2013
#
# PLOTS[0]: Monthly maps
# PLOTS[1]: Period max/mean maps  
# PLOTS[2]: Time series
# PLOTS[3]: Anomaly timeseries
# PLOTS[4]: Zonal plots
# PLOTS[5]: Meridional plots
#
import sys
import numpy as np
import numpy.ma as ma
from numpy import arange,dtype
import matplotlib.pylab as plt
import data_info
import data_GRADS
import data_netCDF
import statistics
import plot_map
import plot_functions
import Emissions_Wetlands_Zonal
#
DEBUG         = sys.argv[1]
iINTER        = sys.argv[2]
JULES_SOURCE  = sys.argv[3]
JULES_RUN     = sys.argv[4]
SDATE         = sys.argv[5]
START_YEAR    = int(sys.argv[6])
END_YEAR      = int(sys.argv[7])
jTRANS        = int(sys.argv[8])
PLOT_OPT      = sys.argv[9]
PLOTS         = sys.argv[10]
FROZEN        = sys.argv[11]
#
STAT_OPTIONS  = [ \
		 ['0: Maximum monthly area for grid square','max' ],\
		 ['1: Mean    monthly area for grid square','mean'],\
		 ['2: Minimum monthly area for grid square','min' ] \
		]
#
#if FROZEN == 'N' and JULES_SOURCE != 'G' and JULES_SOURCE != 'JULES_DC' and JULES_SOURCE != 'JULES_WFD':
if FROZEN == 'N' and JULES_SOURCE not in ('G','JULES_DC','JULES_WFD','JASMIN'):
	print 'Incompatible option - switching to Frozen'
	FROZEN        = 'Y'
#
if PLOT_OPT == '1':
	iDISPLAY = 'Y'
else:
	iDISPLAY = 'N'
#
JULES_OPT     = [['0: Original',''],['1: Masked','_Masked'],['2: EO-based wetland fraction','_EO']]
#
if iINTER == 'Y':
	print("\nJULES Wetland Options")
	for i in range(len(JULES_OPT)):
		print(JULES_OPT[i][0])
#
	print('\nInput Wetland option: ')
	iOPT	  = input()
else:
	iOPT	  = int(sys.argv[12])
#
WET_OPT,WET_CORR=data_info.data_WETLAND_EMISS_info()
TRANS_REGS=data_info.data_TRANSCOM_info()
#
if jTRANS == 0:
	DOMAIN     = 'Global'
else:
	DOMAIN     = TRANS_REGS[jTRANS]
#
if iINTER == 'Y':
	print("Correction option")
	for i in range(len(WET_CORR)):
		print(WET_CORR[i][0])
#
	print('\nInput correction option: ')
	iWET_CORR     = input()
else:
	iWET_CORR     = int(sys.argv[13])
#
if iINTER == 'Y':
	print('\nUse array calculation (0/1): ')
	iARRAY    = input()
else:
	iARRAY    = int(sys.argv[14])
#
# Get metadata
#
SET_UNDER      = 'white'
SET_OVER       = '#800000' # maroon
ORIENTATION    = 'Portrait'
FONTSIZES_MAP  = [ 12,12,14,16 ]
FONTSIZES_PLOT = [  8,10,10,12 ]
WIDTH0         = 9.0
HEIGHT0        = 6.0
RESOLUTION     = 'c'
ASPECT         = 'False'
#
PLOT_LABEL     = 'Wetland fraction (\% of grid square)'
DATA_MIN       = 0.0
DATA_MAX       = 1.0
DATA_INC       = 0.1
NLEVELS        = int((DATA_MAX-DATA_MIN)/DATA_INC)
CLEVELS        = DATA_MIN+DATA_INC*arange(NLEVELS+1)
COLOURS        = [ '#1e3cff', '#00a0ff', '#00c8c8', \
		   '#00d28c', '#00dc00', '#a0e632', '#e6dc32', '#e6af2d', \
		   '#f08228', '#fa3c3c', '#f00082' ]
#
NYEARS     = END_YEAR-START_YEAR+1
SYEARS     = '%04d_%04d' % (START_YEAR,END_YEAR)
NMONTHS    = 12
NTIMES     = NYEARS*NMONTHS
SRESOL     = '0.5'
MISS_DATA  = -999.9
NTRANS     = len(TRANS_REGS)+1
START_EO   = 1993
END_EO     = 2007
#
if NYEARS < 20:
	XINC      = 1
else:
	XINC      = 5
#
XTICKS     = arange(START_YEAR,END_YEAR+XINC,XINC)
#
LAT_START  = -90.0
LAT_END    =  90.0
DEL_LAT    =  30.0
DLAT       = LAT_END-LAT_START
#
LONG_START =   0.0
LONG_END   = 360.0
DEL_LONG   =  30.0
DLONG      = LONG_END-LONG_START
#
if LONG_END == 360.0:
	LONG_PLOTS = -180.0
	LONG_PLOTE =  180.0
#
LONG       = np.arange(LONG_PLOTS,LONG_PLOTE+DEL_LONG,DEL_LONG)
LAT        = np.arange(LAT_START,LAT_END+DEL_LAT,DEL_LAT )
#
RESOL_LONG     = float(SRESOL)
RESOL_LAT      = float(SRESOL)
NLAT           = int((LAT_END-LAT_START)/RESOL_LAT)
NLONG          = int((LONG_END-LONG_START)/RESOL_LONG)
LONG_MAP       = LONG_PLOTS+RESOL_LONG*(0.5+arange(NLONG))
LAT_MAP        = LAT_START +RESOL_LAT*(0.5+arange(NLAT))
#
# Get grid call area
#
NETCDF_DIR    = '/prj/ALANIS/UM_Modelling/EMISSIONS/'
DATA_NAME     = 'cell_area'
FILE_CDF_IN   = NETCDF_DIR + 'CMIP5/a_Downloads_EDGAR/area_0.5x0.5.nc'
print FILE_CDF_IN
DIMS,AREA = data_netCDF.data_netCDF_array2D(FILE_CDF_IN,DATA_NAME)
#
# Get EO wetland fraction
#
TIME_NAME     = 'time'
DATA_NAME     = 'fwetl'
FILE_CDF_EO   = '/prj/ALANIS/ALANIS_Products/Regional_Wetland/Global_1993_2007/wetland_new_1993_2007_Global_fraction_0.50.nc'
print FILE_CDF_EO
DIMS,WET_FRAC_EO_INPUT,TIME_DATA = data_netCDF.data_netCDF_array(FILE_CDF_EO,DATA_NAME,TIME_NAME)
WET_FRAC_EO   = ma.getdata(WET_FRAC_EO_INPUT)
#
# Get JULES wetland fraction
#
if JULES_SOURCE == 'G':
	START_JUL     = 1980
	END_JUL       = 2009
	data_GRADS.setup_GRADS('Linux')
#
	if FROZEN == 'Y':
		DATA_NAME     = 'fwetl'
		FILE_ENDING   = JULES_RUN+'_grads_frozen'
	else:
		DATA_NAME     = 'fwetlUnf'
		FILE_ENDING   = JULES_RUN+'_grads_unfrozen'
#
	FILE_CTL      = '/prj/ALANIS/jules_data/jules_output/'+JULES_RUN+'/monthly_nc/m0'+JULES_RUN[1:]+'.monthly.ctl'
	print FILE_CTL
	LONG_GRADS,LAT_GRADS,STIME,WET_FRAC=data_GRADS.data_GRADS_var(FILE_CTL,DATA_NAME,DEBUG)
#
elif JULES_SOURCE == 'JULES_NG':
	GRIDDED       = 'LAND_MON'
	START_JUL     = 1990
	END_JUL       = 2009
	DATA_NAME     = 'fwetl'
	LAT_NAME      = 'latitude'
	LONG_NAME     = 'longitude'
	FILE_CDF_JUL  = '/prj/ALANIS/UM_Modelling/EMISSIONS/a_JULES_Gedney/WFDEI79103h.monthly.XXXXXX.nc'
	FILE_ENDING   = 'WFDEI79103h_'+SDATE
	print FILE_CDF_JUL
	DIMS,WET_FRAC_JUL = data_netCDF.data_netCDF_array_land_var \
		(FILE_CDF_JUL,DATA_NAME,LAT_NAME,LONG_NAME, \
		 START_JUL,END_JUL,NLONG,NLAT,LONG_START,LAT_START,MISS_DATA,GRIDDED,DEBUG)
	SDATE         = 'JULES_Gedney'
#
elif JULES_SOURCE == 'JULES_DC':
	GRIDDED       = 'LAND_ANN'
	START_JUL     = 1980
	END_JUL       = 2009
#
	if FROZEN == 'Y':
		DATA_NAME     = 'fwetl'
		FILE_ENDING   = 'm45_netCDF_frozen'
	else:
		DATA_NAME     = 'fwetlUnf'
		FILE_ENDING   = 'm45_netCDF_unfrozen'
#
	LAT_NAME      = 'lat'
	LONG_NAME     = 'lon'
	FILE_CDF_JUL  = '/prj/ALANIS/jules_data/jules_output/m45/monthly_nc/m045_monthly_XXXX.nc'
	print FILE_CDF_JUL
	DIMS,WET_FRAC_JUL = data_netCDF.data_netCDF_array_land_var \
		(FILE_CDF_JUL,DATA_NAME,LAT_NAME,LONG_NAME, \
		 START_JUL,END_JUL,NLONG,NLAT,LONG_START,LAT_START,MISS_DATA,GRIDDED,DEBUG)
	SDATE         = 'JULES_m45'
#
elif JULES_SOURCE == 'JULES_WFD':
	START_JUL     = 1990
	END_JUL       = 2009
#
	if FROZEN == 'Y':
		DATA_NAME     = 'fwetl'
		FILE_ENDING   = 'wfdei_1990_2009_wetland_fraction_global_20140903'
	else:
		DATA_NAME     = 'fwetlUnf'
		FILE_ENDING   = 'wfdei_1990_2009_wetland_fraction_unfrozen_global_20140903'
#
	FILE_CDF_JUL  = NETCDF_DIR + 'a_JULES_'+SDATE+'/JULES_'+FILE_ENDING+'.nc'
	print FILE_CDF_JUL
	DIMS,WET_FRAC_JUL,TIME_DATA = data_netCDF.data_netCDF_array(FILE_CDF_JUL,DATA_NAME,TIME_NAME)
	WET_FRAC_JUL = plot_map.switch_long_time(WET_FRAC_JUL)
	FILE_ENDING   = FILE_ENDING.replace('wetland_fraction_','')
	SDATE         = 'JULES_Gedney'
#
elif JULES_SOURCE == 'JASMIN':
	GRIDDED       = 'LAND_ALL'
	START_JUL     = 1980
	END_JUL       = 2014
	LAT_NAME      = 'latitude'
	LONG_NAME     = 'longitude'
	if JULES_RUN == 'WFDEI':
		FILE_ENDING   = 'JASMIN_'+SDATE
		SDATE         = 'JASMIN_WFD_EI'
		FILE_CDF_JUL  = \
			'/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/JULES_WFDEI_nti_NG-HWSD.monthly_wetl.nc'
	if JULES_RUN == 'NCEP_CRU':
		FILE_ENDING   = 'JASMIN_'+SDATE
		SDATE         = 'JASMIN_NCEP_CRU'
		FILE_CDF_JUL  = \
			'/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/NCEP_CRU_global/JULES_v42_CRUNCEP_HWSD_DougDiag_newtopo.monthly_wetl.nc'
#
	if FROZEN == 'Y':
		DATA_NAME     = 'fwetl_all'
		FILE_ENDING   = FILE_ENDING+'_wetland_fraction_global'
	else:
		DATA_NAME     = 'fwetl_unf'
		FILE_ENDING   = FILE_ENDING+'_wetland_fraction_unfrozen_global'
#
	print FILE_CDF_JUL
	print DATA_NAME
	DIMS,WET_FRAC_JUL = data_netCDF.data_netCDF_array_land_var \
		(FILE_CDF_JUL,DATA_NAME,LAT_NAME,LONG_NAME, \
		 START_JUL,END_JUL,NLONG,NLAT,LONG_START,LAT_START,MISS_DATA,GRIDDED,DEBUG)
#
else:
	START_JUL     = 1993
	END_JUL       = 2009
	DATA_NAME     = 'fwetl'
	FILE_ENDING   = JULES_RUN+JULES_OPT[iOPT][1]+'_1993_2009'+WET_CORR[iWET_CORR][1]
	FILE_CDF_JUL  = NETCDF_DIR + 'a_WETLANDS_'+SDATE+'/Wetland_Fraction_JULES_'+FILE_ENDING+'.nc'
	print FILE_CDF_JUL
	DIMS,WET_FRAC_JUL,TIME_DATA = data_netCDF.data_netCDF_array(FILE_CDF_JUL,DATA_NAME,TIME_NAME)
	FILE_ENDING   = JULES_RUN+JULES_OPT[iOPT][1]+'_'+SYEARS+WET_CORR[iWET_CORR][1]
#
LEGEND_JULES  = ('JULES: '+JULES_RUN+' '+JULES_OPT[iOPT][1]+' '+WET_CORR[iWET_CORR][1]).replace('_',' ').replace('single','').replace('EO','GIEMS')
print LEGEND_JULES
#
LEGENDS_MULTI = ['GIEMS',LEGEND_JULES]
#
# Get TRANSCOM regions
#
DATA_NAME     = 'transcom_regions'
FILE_TRANS    = '/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_'+SRESOL+'.nc'
print FILE_TRANS
DIMS,TRANSCOM=data_netCDF.data_netCDF_array2D(FILE_TRANS,DATA_NAME)
#
NLAT       = WET_FRAC_EO.shape[1]
NLONG      = WET_FRAC_EO.shape[2]
NFILES     = 2
#
WET_FRAC_ALL = np.zeros((NFILES,NTIMES,NLAT,NLONG))
WET_FRAC_ALL[:,:,:,:] = MISS_DATA
#
WET_AREA_ALL = np.zeros((NFILES,NTIMES,NLAT,NLONG))
WET_AREA_EXT = np.zeros((NFILES,NLAT,NLONG))
WET_TRAN_EO  = np.zeros((NTIMES,NTRANS))
WET_TRAN_JUL = np.zeros((NTIMES,NTRANS))
WET_ANN_EO   = np.zeros((NYEARS,NTRANS))
WET_ANN_JUL  = np.zeros((NYEARS,NTRANS))
WET_SUM_EO   = np.zeros((NYEARS,NTRANS))
WET_SUM_JUL  = np.zeros((NYEARS,NTRANS))
#
if JULES_SOURCE == 'G':
	NLAT_JULES   = WET_FRAC.shape[1]
	iLAT_START   = int((LAT_GRADS[0]+90.0-RESOL_LAT/2.0)/RESOL_LAT)
	WET_FRAC_JUL = np.zeros((WET_FRAC.shape[0],NLAT,NLONG))
	WET_FRAC_JUL[:,:,:] = MISS_DATA
	WET_FRAC_JUL[:,iLAT_START:iLAT_START+NLAT_JULES,0:NLONG]=WET_FRAC[:,0:NLAT_JULES,1:NLONG+1]
#
WET_FRAC_JUL  = ma.getdata(WET_FRAC_JUL)
#
if LONG_END == 360.0:
        WET_FRAC_JUL = plot_map.switch_long_time(WET_FRAC_JUL)
	WET_FRAC_EO  = plot_map.switch_long_time(WET_FRAC_EO)
#
# Derive wetland area
#
iTIME        = 0
#
for iYEAR in range(NYEARS):
#
	YEAR         = iYEAR+START_YEAR
	iSTART       = iYEAR*12
	iEND         = iSTART+NMONTHS
#
# Copy from input to specific year in data array
#
	if YEAR >= START_EO and YEAR <= END_EO:
		iSTART_EO  = (YEAR-START_EO)*NMONTHS
		iEND_EO    = iSTART_EO+NMONTHS
		print iYEAR,YEAR,iSTART_EO,iEND_EO
		WET_FRAC_ALL[0,iSTART:iEND,:,:] = WET_FRAC_EO[iSTART_EO:iEND_EO,:,:]
#
	if YEAR >= START_JUL and YEAR <= END_JUL:
		iSTART_JUL = (YEAR-START_JUL)*NMONTHS
		iEND_JUL   = iSTART_JUL+NMONTHS
		print iYEAR,YEAR,iSTART_JUL,iEND_JUL
		WET_FRAC_ALL[1,iSTART:iEND,:,:] = WET_FRAC_JUL[iSTART_JUL:iEND_JUL,:,:]
#
	for iMONTH in range(NMONTHS):
#
		SYEAR        = '%4d%02d' % (YEAR,iMONTH+1)
#
		if PLOTS[0] == 'Y':
#
			MAP_TYPE   = 'Map'
#
			DATA_MAP    = WET_FRAC_ALL[0,iTIME,:,:]
			DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
			FILE_END   = SYEAR+WET_CORR[iWET_CORR][1]
			if FROZEN == 'Y':
				FILE_END   = FILE_END+'_frozen'
			else:
				FILE_END   = FILE_END+'_unfrozen'
#
			if jTRANS != 0: FILE_END   = FILE_END+'_'+('TRANS%02d' % jTRANS)
#
			FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_EO_'+FILE_END+'.png'
			print FILE_PLOT
			PLOT_TITLE = 'EO wetland fraction for ' + SYEAR
#
			plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
				LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
			DATA_MAP    = WET_FRAC_ALL[1,iTIME,:,:]
			DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
			FILE_END   = JULES_RUN+JULES_OPT[iOPT][1]+'_'+SYEAR+WET_CORR[iWET_CORR][1]
			if FROZEN == 'Y':
				FILE_END   = FILE_END+'_frozen'
			else:
				FILE_END   = FILE_END+'_unfrozen'
#
			if jTRANS != 0: FILE_END   = FILE_END+'_'+('TRANS%02d' % jTRANS)
			FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_JULES_'+FILE_END+'.png'
			print FILE_PLOT
			PLOT_TITLE = 'JULES wetland fraction for ' + SYEAR
#
#			if LONG_END == 360.0:
#			        DATA_MAP    = plot_map.switch_long(DATA_MAP)
#
			plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
				LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
				RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
		WET_AREA_ALL[0,iTIME,:,:] = WET_FRAC_ALL[0,iTIME,:,:]*AREA[:,:]
		WET_AREA_ALL[1,iTIME,:,:] = WET_FRAC_ALL[1,iTIME,:,:]*AREA[:,:]
		WET_AREA_ALL[WET_AREA_ALL < 0] = 0.0
		 
		if iARRAY == 1:
#
			TEMP_EO    = WET_AREA_ALL[0,iTIME,:,:].copy()
			TEMP_JUL   = WET_AREA_ALL[1,iTIME,:,:].copy()
#
			for iTRANS in range(NTRANS):
#
				if iTRANS == 0:
					WET_TRAN_EO[iTIME,iTRANS] = TEMP_EO[TEMP_EO > 0].sum()
				else:
					WET_TRAN_EO[iTIME,iTRANS] = TEMP_EO[(TEMP_EO > 0) & (TRANSCOM == iTRANS) & (~np.isnan(TEMP_EO))].sum()
#
				if iTRANS == 0:
					WET_TRAN_JUL[iTIME,iTRANS] = TEMP_JUL[TEMP_JUL > 0].sum()
				else:
					WET_TRAN_JUL[iTIME,iTRANS] = TEMP_JUL[(TEMP_JUL > 0) & (TRANSCOM == iTRANS) & (~np.isnan(TEMP_JUL))].sum()
#
		else:
			WET_TRAN_EO[iTIME,:]  = 0.0
			WET_TRAN_JUL[iTIME,:] = 0.0
#
			for iLAT in range(NLAT):
				for iLONG in range(NLONG):
#
					if WET_FRAC_ALL[0,iTIME,iLAT,iLONG] > 0.0:
						WET_TRAN_EO[iTIME,0]        = WET_TRAN_EO[iTIME,0]+ \
							WET_FRAC_ALL[0,iTIME,iLAT,iLONG]* \
							AREA[iLAT,iLONG]
#
					if WET_FRAC_ALL[1,iTIME,iLAT,iLONG] > 0.0:
						WET_TRAN_JUL[iTIME,0]       = WET_TRAN_JUL[iTIME,0]+ \
							WET_FRAC_ALL[1,iTIME,iLAT,iLONG]* \
							AREA[iLAT,iLONG]
#
			for iLAT in range(NLAT):
				for iLONG in range(NLONG):
					iTRANS    = TRANSCOM[iLAT,iLONG]
#
					if WET_FRAC_ALL[0,iTIME,iLAT,iLONG] > 0.0:
						WET_TRAN_EO[iTIME,iTRANS]   = WET_TRAN_EO[iTIME,iTRANS]+ \
							WET_FRAC_ALL[0,iTIME,iLAT,iLONG]* \
							AREA[iLAT,iLONG]
#
					if WET_FRAC_ALL[1,iTIME,iLAT,iLONG] > 0.0:
						WET_TRAN_JUL[iTIME,iTRANS]  = WET_TRAN_JUL[iTIME,iTRANS]+ \
							WET_FRAC_ALL[1,iTIME,iLAT,iLONG]* \
							AREA[iLAT,iLONG]
#
			print SYEAR
			print WET_TRAN_EO[iTIME,:]
			print WET_TRAN_JUL[iTIME,:]
#
		iTIME += 1
#
	TEXT_EO    = '%4d EO   : ' % (iYEAR+START_YEAR)
	TEXT_JUL   = '%4d JULES: ' % (iYEAR+START_YEAR)
#
	for iTRANS in range(NTRANS):
		TEMP = WET_TRAN_EO[iSTART:iEND,iTRANS]
		WET_ANN_EO[iYEAR,iTRANS]  = TEMP[~np.isnan(TEMP)].mean() 
		WET_SUM_EO[iYEAR,iTRANS]  = TEMP[~np.isnan(TEMP)].sum()
		TEMP = WET_TRAN_JUL[iSTART:iEND,iTRANS]
		WET_ANN_JUL[iYEAR,iTRANS] = TEMP[~np.isnan(TEMP)].mean() 
		WET_SUM_JUL[iYEAR,iTRANS] = TEMP[~np.isnan(TEMP)].sum() 
#
		TEXT_EO  = TEXT_EO  + ('%8.2e ' % (WET_SUM_EO[iYEAR,iTRANS]))
		TEXT_JUL = TEXT_JUL + ('%8.2e ' % (WET_SUM_JUL[iYEAR,iTRANS]))
#
	TEXT_EO  = TEXT_EO  + ('%8.2e ' % (WET_SUM_EO[iYEAR,1:].sum()))
	TEXT_JUL = TEXT_JUL + ('%8.2e ' % (WET_SUM_JUL[iYEAR,1:].sum()))
#
	print TEXT_EO
	print TEXT_JUL
#
MAX_AREA_EO  = WET_TRAN_EO[~np.isnan(WET_TRAN_EO)].max()
MAX_AREA_JUL = WET_TRAN_JUL[~np.isnan(WET_TRAN_JUL)].max()
MAX_AREA     = max(MAX_AREA_EO,MAX_AREA_JUL)
#
WET_TRAN_EO[WET_TRAN_EO <= 0.0]   = float('nan')
WET_TRAN_JUL[WET_TRAN_JUL <= 0.0] = float('nan')
WET_ANN_EO[WET_ANN_EO <= 0.0]     = float('nan')
WET_ANN_JUL[WET_ANN_JUL <= 0.0]   = float('nan')
#
if PLOTS[1] == 'Y':
#
	MAP_TYPE   = 'Map'
#
	PLOT_EO_MEAN    = np.zeros((NLAT,NLONG))
	PLOT_EO_MAX     = np.zeros((NLAT,NLONG))
	PLOT_JULES_MEAN = np.zeros((NLAT,NLONG))
	PLOT_JULES_MAX  = np.zeros((NLAT,NLONG))
#
	if jTRANS == 0:
		for iLAT in range(NLAT):
			for iLONG in range(NLONG):
#
				DATA_EO                     = WET_FRAC_EO[:,iLAT,iLONG]
				PLOT_EO_MEAN[iLAT,iLONG]    = DATA_EO[~np.isnan(DATA_EO)].mean()
				PLOT_EO_MAX[iLAT,iLONG]     = DATA_EO[~np.isnan(DATA_EO)].max()
#
				DATA_JULES                  = WET_FRAC_JUL[:,iLAT,iLONG]
				PLOT_JULES_MEAN[iLAT,iLONG] = DATA_JULES[~np.isnan(DATA_JULES)].mean()
				PLOT_JULES_MAX[iLAT,iLONG]  = DATA_JULES[~np.isnan(DATA_JULES)].max()
	else:
		for iLAT in range(NLAT):
			for iLONG in range(NLONG):
				if TRANSCOM[iLAT,iLONG] == jTRANS:
#
					DATA_EO                     = WET_FRAC_EO[:,iLAT,iLONG]
					PLOT_EO_MEAN[iLAT,iLONG]    = DATA_EO[~np.isnan(DATA_EO)].mean()
					PLOT_EO_MAX[iLAT,iLONG]     = DATA_EO[~np.isnan(DATA_EO)].max()
#
					DATA_JULES                  = WET_FRAC_JUL[:,iLAT,iLONG]
					PLOT_JULES_MEAN[iLAT,iLONG] = DATA_JULES[~np.isnan(DATA_JULES)].mean()
					PLOT_JULES_MAX[iLAT,iLONG]  = DATA_JULES[~np.isnan(DATA_JULES)].max()
				else:
					PLOT_EO_MEAN[iLAT,iLONG]    = MISS_DATA
					PLOT_EO_MAX[iLAT,iLONG]     = MISS_DATA
					PLOT_JULES_MEAN[iLAT,iLONG] = MISS_DATA
					PLOT_JULES_MAX[iLAT,iLONG]  = MISS_DATA
#
	FYEAR      = str(START_YEAR)+'_'+str(END_YEAR)
	SYEAR      = FYEAR.replace('_','-')
	FILE_END   = FYEAR+WET_CORR[iWET_CORR][1]
	if FROZEN == 'Y':
		FILE_END   = FILE_END+'_frozen'
	else:
		FILE_END   = FILE_END+'_unfrozen'
#
	if jTRANS != 0: FILE_END   = FILE_END+'_'+('TRANS%02d' % jTRANS)
#
# EO - max grid square
#
	DATA_MAP    = PLOT_EO_MAX
	DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_EO_'+FILE_END+'_max.png'
	print FILE_PLOT
	PLOT_TITLE = 'Maximum wetland fraction in grid square for '+SYEAR+' (EO)'
#
#	if LONG_END == 360.0:
#	        DATA_MAP    = plot_map.switch_long(DATA_MAP)
#
	plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
		LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
		RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
# EO - mean grid square
#
	DATA_MAP    = PLOT_EO_MEAN
	DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_EO_'+FILE_END+'_mean.png'
	print FILE_PLOT
	PLOT_TITLE = 'Mean wetland fraction in grid square for '+SYEAR+' (EO)'
#
#	if LONG_END == 360.0:
#	        DATA_MAP    = plot_map.switch_long(DATA_MAP)
#
	plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
		LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
		RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
# JULES - max grid square
#
	FILE_END   = JULES_RUN+JULES_OPT[iOPT][1]+'_'+FYEAR+WET_CORR[iWET_CORR][1]
	if FROZEN == 'Y':
		FILE_END   = FILE_END+'_frozen'
	else:
		FILE_END   = FILE_END+'_unfrozen'
#
	if jTRANS != 0: FILE_END   = FILE_END+'_'+('TRANS%02d' % jTRANS)
#
	DATA_MAP    = PLOT_JULES_MAX
	DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_JULES_'+FILE_END+'_max.png'
	print FILE_PLOT
	PLOT_TITLE = 'Maximum wetland fraction in grid square for '+SYEAR+' (JULES '+JULES_RUN.replace('_',' ')+JULES_OPT[iOPT][1].replace('_','')+')'
#
#	if LONG_END == 360.0:
#	        DATA_MAP    = plot_map.switch_long(DATA_MAP)
#
	plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
		LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
		RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
# JULES - mean grid square
#
	DATA_MAP    = PLOT_JULES_MEAN
	DATA_MAP[DATA_MAP <= 0.0] = MISS_DATA
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Fraction_JULES_'+FILE_END+'_mean.png'
	print FILE_PLOT
	PLOT_TITLE = 'Mean wetland fraction in grid square for '+SYEAR+' (JULES '+JULES_RUN.replace('_',' ')+JULES_OPT[iOPT][1].replace('_','')+')'
#
#	if LONG_END == 360.0:
#	        DATA_MAP    = plot_map.switch_long(DATA_MAP)
#
	plot_map.plot_map3(DATA_MAP,LONG_PLOTS,LONG_PLOTE,LONG,LONG_MAP, \
		LAT_START,LAT_END,LAT,LAT_MAP,CLEVELS,COLOURS,MAP_TYPE,WIDTH0,HEIGHT0,ASPECT, \
		RESOLUTION,PLOT_TITLE,PLOT_LABEL,FILE_PLOT,iDISPLAY,FONTSIZES_MAP,SET_UNDER,SET_OVER,DEBUG)
#
if PLOTS[2] == 'Y':
#
	if jTRANS == 0:
		FILE_ENDING    = FILE_ENDING+'_global'
	else:
		FILE_ENDING    = FILE_ENDING+'_'+('TRANS%02d' % jTRANS)
#
	DATA_MIN   = 0.0
	DATA_MAX   = 1.0
	DATA_INC   = 0.1
	NYTICKS    = int((DATA_MAX-DATA_MIN)/DATA_INC)+1
	YTICKS     = DATA_MIN + DATA_INC*arange(NYTICKS)
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_EO_JULES_'+FILE_ENDING+'.png'
	PLOT_TITLE = 'Normalised EO and modelled wetland area for '+DOMAIN
	YLABEL     = 'Normalised wetland area'
	XLABEL     = 'Year'
	LEGEND     = ['EO (max area):    '+('%.2e km$^{2}$' % (MAX_AREA_EO/1.00E+06)),'JULES (max area): '+('%.2e km$^{2}$' % (MAX_AREA_JUL//1.00E+06))]
	LEGEND_POS = 3
	PLOT_CODE  = ['r-','b-']
	NDATA      = 2
#
	XDATA      = np.zeros((2,NTIMES))
	YDATA      = np.zeros((2,NTIMES))
#
	print YDATA.shape,WET_TRAN_EO.shape,WET_TRAN_JUL.shape
#
	TIME       = float(START_YEAR)+1.0/24+arange(NTIMES)/12.0
	XDATA[0,:] = TIME
	XDATA[1,:] = TIME
	YDATA[0,:] = WET_TRAN_EO[:,jTRANS]/MAX_AREA
	YDATA[1,:] = WET_TRAN_JUL[:,jTRANS]/MAX_AREA
#
	plot_functions.Plot_TimeSeries2(NDATA,XDATA,YDATA,START_YEAR,END_YEAR,XTICKS,XLABEL, \
		DATA_MIN,DATA_MAX,YTICKS,YLABEL,PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_Anomaly_EO_JULES_'+FILE_ENDING+'_annual.png'
	PLOT_TITLE = 'Anomaly in EO and modelled annual mean wetland area for '+DOMAIN
	YLABEL     = 'Normalised percentage anomaly'
	PLOT_CODE  = ['ro-','bo-']
	LEGEND     = ['EO','JULES']
#
	DATA_MIN   = -15
	DATA_MAX   =  15
	DATA_INC   =   5
	NYTICKS    = int((DATA_MAX-DATA_MIN)/DATA_INC)+1
	YTICKS     = DATA_MIN + DATA_INC*arange(NYTICKS)
#
	TEMP_EO    = WET_ANN_EO[:,jTRANS].copy()
	TEMP_JULES = WET_ANN_JUL[:,jTRANS].copy()
#
	MEAN_EO    = TEMP_EO[~np.isnan(TEMP_EO)].mean()
	MEAN_JULES = TEMP_JULES[~np.isnan(TEMP_JULES)].mean()
#
	XDATA      = np.zeros((2,NYEARS))
	YDATA      = np.zeros((2,NYEARS))
#
	XDATA[0,:] = float(START_YEAR)+0.5+arange(NYEARS)
	XDATA[1,:] = XDATA[0,:]
	YDATA[0,:] = 100.0*(TEMP_EO - MEAN_EO)/MEAN_EO
	YDATA[1,:] = 100.0*(TEMP_JULES - MEAN_JULES)/MEAN_JULES
#
	plot_functions.Plot_TimeSeries2(NDATA,XDATA,YDATA,START_YEAR,END_YEAR,XTICKS,XLABEL, \
		DATA_MIN,DATA_MAX,YTICKS,YLABEL,PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
if PLOTS[3] == 'Y':
#
	XDATA      = np.zeros((2,NTIMES))
	YDATA      = np.zeros((2,NTIMES))
#
	TRANS_EO_R  = np.zeros(NTIMES)
	TRANS_JUL_R = np.zeros(NTIMES)
	TRANS_EO_A  = np.zeros(NTIMES)
	TRANS_JUL_A = np.zeros(NTIMES)
	TIME_ANOM   = np.zeros(NTIMES)
#
	TRANS_EO_A[:]  = float('nan')
	TRANS_JUL_A[:] = float('nan')
	TRANS_EO_R[:]  = float('nan')
	TRANS_JUL_R[:] = float('nan')
#
	for iMONTH in range(NTIMES):
		TIME_ANOM[iMONTH] = START_YEAR+(iMONTH+0.5)/12
		if iMONTH > 5 and iMONTH <= NTIMES-6:
			TRANS_EO_R[iMONTH]  = WET_TRAN_EO[iMONTH-6:iMONTH+6,jTRANS].mean()
			TRANS_JUL_R[iMONTH] = WET_TRAN_JUL[iMONTH-6:iMONTH+6,jTRANS].mean()
		else:
			TRANS_EO_R[iMONTH]  = float('nan')
			TRANS_JUL_R[iMONTH] = float('nan')
#
	if DEBUG == 'Y':
		print
		print NTIMES,jTRANS
		print WET_TRAN_EO[:,jTRANS]
		print TRANS_EO_R
		print WET_TRAN_JUL[:,jTRANS]
		print TRANS_JUL_R
#
	print
	print MAX_AREA_EO,MAX_AREA_JUL,MAX_AREA
	print 'Input maximum area for plot: '
	MAX_PLOT    = float(input())
#
	TRANS_EO_M  = TRANS_EO_R[~np.isnan(TRANS_EO_R)].mean()
 	TRANS_JUL_M = TRANS_JUL_R[~np.isnan(TRANS_JUL_R)].mean()
	TRANS_EO_A[~np.isnan(TRANS_EO_R)]   = TRANS_EO_R[~np.isnan(TRANS_EO_R)]-TRANS_EO_M
	TRANS_JUL_A[~np.isnan(TRANS_JUL_R)] = TRANS_JUL_R[~np.isnan(TRANS_JUL_R)]-TRANS_JUL_M
#
	XDATA_ALL  = []
	YDATA_ALL  = []
	YMIN_ALL   = []
	YMAX_ALL   = []
	YTICKS_ALL = []
#
	XDATA[0,:] = TIME_ANOM
	XDATA[1,:] = TIME_ANOM
#
	TEMP_EO                        = WET_TRAN_EO[:,jTRANS].copy()
	TEMP_EO[~np.isnan(TEMP_EO)]    = TEMP_EO[~np.isnan(TEMP_EO)]/MAX_PLOT
	TEMP_EO_M                      = TEMP_EO[~np.isnan(TEMP_EO)].mean()
	TEMP_JUL                       = WET_TRAN_JUL[:,jTRANS].copy()
	TEMP_JUL[~np.isnan(TEMP_JUL)]  = TEMP_JUL[~np.isnan(TEMP_JUL)]/MAX_PLOT
	TEMP_JUL_M                     = TEMP_JUL[~np.isnan(TEMP_JUL)].mean()
#
	YDATA_ALL.append([TEMP_EO,TEMP_JUL])
#
	TEMP_EO1                       = WET_TRAN_EO[:,jTRANS].copy()
	TEMP_EO1[~np.isnan(TEMP_EO)]   = TEMP_EO[~np.isnan(TEMP_EO)]-TEMP_EO_M
	TEMP_JUL1                      = WET_TRAN_JUL[:,jTRANS].copy()
	TEMP_JUL1[~np.isnan(TEMP_JUL)] = TEMP_JUL[~np.isnan(TEMP_JUL)]-TEMP_JUL_M
	YDATA_ALL.append([TEMP_EO1,TEMP_JUL1])
#
	YDATA[:,:] = float('nan')
	YDATA[0,:] = TRANS_EO_A/MAX_PLOT
	YDATA[1,:] = TRANS_JUL_A/MAX_PLOT
	YDATA_ALL.append([YDATA[0,:],YDATA[1,:]])
#
	XDATA_ALL  = [[TIME_ANOM,TIME_ANOM],[TIME_ANOM,TIME_ANOM],[TIME_ANOM,TIME_ANOM]]
#
	print YDATA.shape,WET_TRAN_EO.shape,WET_TRAN_JUL.shape
#
	NYTICKS    = 11
	DATA_MIN   = -0.1
	DATA_MAX   =  0.1
	YTICKS     = plot_functions.get_ticks2(DATA_MAX,DATA_MIN,NYTICKS,DEBUG)
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_Anomaly_EO_JULES_'+FILE_ENDING+'.png'
	PLOT_TITLE = 'Anomaly in EO and modelled wetland area for '+DOMAIN
	YLABEL     = 'Anomaly in wetland area'
	XLABEL     = 'Year'
	LEGEND     = ['EO (max area):    '+('%.2e km$^{2}$' % (MAX_AREA_EO/1.00E+06)),'JULES (max area): '+('%.2e km$^{2}$' % (MAX_AREA_JUL//1.00E+06))]
	LEGEND_POS = 3
	PLOT_CODE  = ['r-','b-']
	NDATA      = 2
#
	plot_functions.Plot_TimeSeries2(NDATA,XDATA,YDATA,START_YEAR,END_YEAR,XTICKS,XLABEL, \
		DATA_MIN,DATA_MAX,YTICKS,YLABEL,PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
	NROWS      = 1
	NCOLUMNS   = 3
	NDATA_ALL  = [NDATA,NDATA,NDATA]
	YLABELS_ALL= ['Wetland Area','Wetland Area',YLABEL]
	WIDTH      = 6.0*NCOLUMNS
	HEIGHT     = 6.0*NROWS
#
	DATA_MIN1  =  0.0
	DATA_MAX1  = 1.0
	DATA_MIN2  = -0.5
	DATA_MAX2  =  0.5
	YMIN_ALL   = [DATA_MIN1,DATA_MIN2,DATA_MIN]
	YMAX_ALL   = [DATA_MAX1,DATA_MAX2,DATA_MAX]
	YTICKS_ALL = [plot_functions.get_ticks2(DATA_MAX1,DATA_MIN1,NYTICKS,DEBUG), \
		      plot_functions.get_ticks2(DATA_MAX2,DATA_MIN2,NYTICKS,DEBUG),YTICKS]
#
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_Anomaly_EO_JULES_'+FILE_ENDING+'_multi.png'
	PLOT_TITLE = 'Anomaly in EO and modelled wetland area for '+DOMAIN
	SUBTITLES  = ['Original','Original-Mean','Anomaly in Running Annual Mean']
#
	plot_functions.Plot_General_MultiPlot_VarScale2a(NROWS,NCOLUMNS,NDATA_ALL,XDATA_ALL,YDATA_ALL, \
	START_YEAR,END_YEAR+1,XTICKS,XTICKS,XLABEL, \
	YMIN_ALL,YMAX_ALL,YTICKS_ALL,YTICKS_ALL,YLABELS_ALL,WIDTH,HEIGHT,FONTSIZES_PLOT, \
	PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG)
#
if PLOTS[4] == 'Y' or PLOTS[5] == 'Y':
#
	PLOT_CODES_ZL = ['k-' ,'r-' ,'g-' ,'b-' ,'m-' ,'y-'  ,'c-' ]
	WIDTH         = 10.0
	HEIGHT        =  5.0
	LEGEND_POS    =  2
#
	print
	print "Option for zonal/meriodional wetland plot"
	for iOPT in range(len(STAT_OPTIONS)):
		print STAT_OPTIONS[iOPT][0]
	print '\nInput option: '
	iOPT          = int(input())
#
	FILE_STAT     = STAT_OPTIONS[iOPT][1]
#
	for iLAT in range(NLAT):
		for iLONG in range(NLONG):
			for iFILE in range(NFILES):
				if iOPT == 0:
					WET_AREA_EXT[iFILE,iLAT,iLONG] = WET_AREA_ALL[iFILE,:,iLAT,iLONG].max()/1.0E+12
				elif iOPT == 1:
					WET_AREA_EXT[iFILE,iLAT,iLONG] = WET_AREA_ALL[iFILE,:,iLAT,iLONG].mean()/1.0E+12
				elif iOPT == 2:
					WET_AREA_EXT[iFILE,iLAT,iLONG] = WET_AREA_ALL[iFILE,:,iLAT,iLONG].min()/1.0E+12
#
# PLOTS[4] - Zonal/Latitudinal plot
#
if PLOTS[4] == 'Y':
#
	LEGEND           = []
	PLOT_CODE        = []
#
	for iFILE in range(NFILES):
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
			XZONAL         = np.zeros((NFILES,NLAT_INT))
			YZONAL         = np.zeros((NFILES,NLAT_INT))
#
		else:
			NLAT_INT2      = int(dLAT/RESOL_LAT)
#
# Data for latitudinal plots
#
	TEMP           = Emissions_Wetlands_Zonal.Area_Wetlands_Zonal( \
		NFILES,NLONG,NLAT_INT,NLAT_INT2,dLAT,WET_AREA_EXT,LAT_ZONAL,DEBUG)
#
	for iFILE in range(NFILES):
#
		XZONAL[iFILE,:]= LAT_ZONAL
		YZONAL[iFILE,:]=TEMP[iFILE,:]
#
# Latitudinal plot
#
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_ZL[iFILE])
#
# Zonal plot (Latitude)
#
		if iFILE == 0:
			print 'Minimum and maximum data values are :',YZONAL.min(),YZONAL.max()
			print 'Input maximum data value: '
			YMAX             = float(input())
			YTICKS           = (YMAX/10.0)*arange(11)
#
	XLABEL      = 'Latitude'
	YLABEL      = 'Wetland Area (10 $^{6}$ km$^{2}$ per '+('%.1f' % dLAT)+' degree Latitude Band)'
	TITLE       = 'Zonal sum of the monthly '+FILE_STAT+' wetland area'
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_'+FILE_STAT+'_EO_JULES_'+FILE_ENDING+'_zonal_multi.png'
	print FILE_PLOT
#
	plot_functions.Plot_Zonal(NFILES,XZONAL,YZONAL,LAT_START,LAT_END,LAT_TICKS,XLABEL, \
		0.0,YMAX,YTICKS,YLABEL,TITLE,LEGEND,FONTSIZES_PLOT,PLOT_CODE,PLOT_OPT, \
		FILE_PLOT,DEBUG,WIDTH,HEIGHT,LEGEND_POS)
#
# PLOTS[5] - longitudinal plots
#
if PLOTS[5] == 'Y':
#
	LEGEND           = []
	PLOT_CODE        = []
#
	for iFILE in range(NFILES):
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
			XMERID         = np.zeros((NFILES,NLON_INT))
			YMERID         = np.zeros((NFILES,NLON_INT))
#
		else:
			NLON_INT2      = int(dLON/RESOL_LONG)
#
# Data for longitudinal plots
#
	NDATA_MERID    = 2
	TEMP           = Emissions_Wetlands_Zonal.Area_Wetlands_Meridional( \
		NDATA_MERID,NLAT,NLON_INT,NLON_INT2,dLON,LONG_END,WET_AREA_EXT,LON_MERID,DEBUG)
#
	for iFILE in range(NFILES):
#
		XMERID[iFILE,:]  = LON_MERID
		YMERID[iFILE,:]  = TEMP[iFILE,:]
#
		LEGEND.append(LEGENDS_MULTI[iFILE])
		PLOT_CODE.append(PLOT_CODES_ZL[iFILE])
#
# Zonal plot (Latitude)
#
		if iFILE == 0:
			print 'Minimum and maximum data values are :',YMERID.min(),YMERID.max()
			print 'Input maximum data value: '
			YMAX             = float(input())
			YTICKS           = (YMAX/10.0)*arange(11)
#		
	XLABEL      = 'Longitude'
	YLABEL      = 'Wetland Area (10 $^{6}$ km$^{2}$ per '+('%.1f' % dLON)+' degree Longitude Band)'
	TITLE       = 'Longitudinal sum of the monthly '+FILE_STAT+' wetland area'
	FILE_PLOT  = '/prj/ALANIS/UM_Modelling/PLOTS/Wetlands_'+SDATE+'/Wetland_Area_'+FILE_STAT+'_EO_JULES_'+FILE_ENDING+'_meridional_multi.png'
	print FILE_PLOT
#
	plot_functions.Plot_Zonal(NFILES,XMERID,YMERID,LONG_PLOTS,LONG_PLOTE,LON_TICKS,XLABEL, \
		0.0,YMAX,YTICKS,YLABEL,TITLE,LEGEND,FONTSIZES_PLOT,PLOT_CODE,PLOT_OPT, \
		FILE_PLOT,DEBUG,WIDTH,HEIGHT,LEGEND_POS)
#
