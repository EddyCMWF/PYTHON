# Python module to derive emissions from Wetlands
#
# Garry Hayman
# Centre for Ecology and Hydrology
# July 2012
#
# Contains
#
import numpy as np
from numpy import arange
#
def Emissions_Wetlands_Zonal(NDATASETS,NMONTHS,NLONG,NLAT_INT,NLAT_INT2,dLAT,TOTAL_EMISS,EMISS,LAT_ZONAL,sSPECIES,sYEAR,DEBUG):
#
	EMISS_LAT      = np.zeros((NDATASETS,NLAT_INT))
	TOTAL_EMISS_LAT= np.zeros((NLAT_INT2,NLONG))
#
# Calculate latitudinal contribution
#
	for iDATA in range(NDATASETS):
#
		EMISS_LAT[iDATA,:]     = 0.0
#
		iLAT_END    =  0
		zLAT_END    = -90.0
#
		for iLAT in range(NLAT_INT):
#
			iLAT_START = iLAT_END
			iLAT_END   = iLAT_START + NLAT_INT2
			zLAT_START = zLAT_END
			zLAT_END   = zLAT_START  + dLAT
			if DEBUG == 'Y': print(iLAT,iLAT_START,iLAT_END,zLAT_START,zLAT_END,NLAT_INT,NLAT_INT2) 
#
			for jLAT in range(iLAT_START,iLAT_END):
				kLAT       = jLAT-iLAT_START
				if DEBUG == 'Y': print(kLAT,jLAT)
				TOTAL_EMISS_LAT[kLAT,:] = TOTAL_EMISS[jLAT,:]
#
				for iMONTH in range(NMONTHS):
					EMISS_LAT[iDATA,iLAT] = EMISS_LAT[iDATA,iLAT]+EMISS[iMONTH,jLAT,:].sum()
#
			if DEBUG == 'Y':
				TEXT         = 'Annual '+sSPECIES+' emissions                  '+' for '+sYEAR+' = '+ \
					('%12.4f' % TOTAL_EMISS_LAT[:,:].sum())+' Tg per annum between '+('%5.1f' % zLAT_START)+ \
					' and '+('%5.1f' % zLAT_END)
				print(TEXT)
#
# Aggregate to larger latitudinal interval
#
		print
		Emissions_Wetlands_Aggregate(sSPECIES,sYEAR,LAT_ZONAL,EMISS_LAT[iDATA,:])
#
		if DEBUG == 'Y':
			print(LAT_ZONAL)
			print(EMISS_LAT[iDATA,:])
#
# Return to calling routine
#
	return EMISS_LAT
#
def Emissions_Wetlands_Meridional(NDATASETS,NMONTHS,NLAT,NLON_INT,NLON_INT2,dLON,LONG_END,TOTAL_EMISS,EMISS,LON_ZONAL,sSPECIES,sYEAR,DEBUG):
#
	EMISS_LON      = np.zeros((NDATASETS,NLON_INT))
	TOTAL_EMISS_LON= np.zeros((NLAT,NLON_INT2))
#
# Calculate latitudinal contribution
#
	for iDATA in range(NDATASETS):
#
		EMISS_LON[iDATA,:]     = 0.0
#
		if LONG_END == 180.0:
			iLON_END    =  0
			zLON_END    = -180.0
		else:
			iLON_END    =  0
			zLON_END    =  0.0
#
		for iLON in range(NLON_INT):
#
			iLON_START = iLON_END
			iLON_END   = iLON_START + NLON_INT2
			zLON_START = zLON_END
			zLON_END   = zLON_START  + dLON
			if DEBUG == 'Y': print(iLON,iLON_START,iLON_END,zLON_START,zLON_END,NLON_INT,NLON_INT2) 
#
			for jLON in range(iLON_START,iLON_END):
				kLON       = jLON-iLON_START
				if DEBUG == 'Y': print(kLON,jLON)
				TOTAL_EMISS_LON[:,kLON] = TOTAL_EMISS[:,jLON]
#
				for iMONTH in range(NMONTHS):
					EMISS_LON[iDATA,iLON] = EMISS_LON[iDATA,iLON]+EMISS[iMONTH,:,jLON].sum()
#
			if DEBUG == 'Y':
				TEXT         = 'Annual '+sSPECIES+' emissions                  '+' for '+sYEAR+' = '+ \
					('%12.4f' % TOTAL_EMISS_LON[:,:].sum())+' Tg per annum between '+('%5.1f' % zLON_START)+ \
					' and '+('%5.1f' % zLON_END)
				print(TEXT)
#
		if DEBUG == 'Y':
			print(LON_ZONAL)
			print(EMISS_LON[iDATA,:])
#
# Return to calling routine
#
	return EMISS_LON
#
def Emissions_Wetlands_Aggregate(sSPECIES,sYEAR,LAT,EMISS):
#
	EMISS_EXTRA_SH = EMISS[(LAT <= -30.0)].sum()
	EMISS_TROPICS  = EMISS[(LAT > -30.0) & (LAT <= 30.0)].sum()
	EMISS_EXTRA_NH = EMISS[(LAT >  30.0) & (LAT <= 50.0)].sum()
	EMISS_BOREAL   = EMISS[(LAT >  50.0) & (LAT <= 90.0)].sum()
	EMISS_TOTAL    = EMISS_EXTRA_SH+EMISS_TROPICS+EMISS_EXTRA_NH+EMISS_BOREAL
#	print EMISS_TOTAL,EMISS_EXTRA_SH,EMISS_TROPICS,EMISS_EXTRA_NH,EMISS_BOREAL
#
	TEXT         = 'Annual '+sSPECIES+' emissions                  '+' for '+sYEAR+' = '+ \
		('%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f' % (EMISS_EXTRA_SH,EMISS_TROPICS,EMISS_EXTRA_NH,EMISS_BOREAL,EMISS_TOTAL,EMISS.sum()))+' Tg per annum'
	print(TEXT)
#
# Return to calling routine
#
	return
#
def Area_Wetlands_Zonal(NDATASETS,NLONG,NLAT_INT,NLAT_INT2,dLAT,WET_AREA,LAT_ZONAL,DEBUG):
#
	WET_AREA_LAT      = np.zeros((NDATASETS,NLAT_INT))
#
# Calculate latitudinal contribution
#
	for iDATA in range(NDATASETS):
#
		WET_AREA_LAT[iDATA,:]     = 0.0
#
		iLAT_END    =  0
		zLAT_END    = -90.0
#
		for iLAT in range(NLAT_INT):
#
			iLAT_START = iLAT_END
			iLAT_END   = iLAT_START + NLAT_INT2
			zLAT_START = zLAT_END
			zLAT_END   = zLAT_START  + dLAT
			if DEBUG == 'Y': print iLAT,iLAT_START,iLAT_END,zLAT_START,zLAT_END,NLAT_INT,NLAT_INT2 
#
			for jLAT in range(iLAT_START,iLAT_END):
				kLAT       = jLAT-iLAT_START
				if DEBUG == 'Y': print iLAT,kLAT,jLAT
				TEMP_AREA = WET_AREA[iDATA,jLAT,:]
				WET_AREA_LAT[iDATA,iLAT] = WET_AREA_LAT[iDATA,iLAT]+TEMP_AREA.sum()
				if DEBUG == 'Y': print TEMP_AREA[TEMP_AREA > 0],WET_AREA_LAT[iDATA,iLAT]
#
		if DEBUG == 'Y':
			print(LAT_ZONAL)
			print(WET_AREA_LAT[iDATA,:])
#
# Return to calling routine
#
	return WET_AREA_LAT
#
def Area_Wetlands_Meridional(NDATASETS,NLAT,NLON_INT,NLON_INT2,dLON,LONG_END,WET_AREA,LON_ZONAL,DEBUG):
#
	WET_AREA_LON      = np.zeros((NDATASETS,NLON_INT))
#
# Calculate latitudinal contribution
#
	for iDATA in range(NDATASETS):
#
		WET_AREA_LON[iDATA,:]     = 0.0
#
		if LONG_END == 180.0:
			iLON_END    =  0
			zLON_END    = -180.0
		else:
			iLON_END    =  0
			zLON_END    =  0.0
#
		for iLON in range(NLON_INT):
#
			iLON_START = iLON_END
			iLON_END   = iLON_START + NLON_INT2
			zLON_START = zLON_END
			zLON_END   = zLON_START  + dLON
			if DEBUG == 'Y': print iLON,iLON_START,iLON_END,zLON_START,zLON_END,NLON_INT,NLON_INT2 
#
			for jLON in range(iLON_START,iLON_END):
				kLON       = jLON-iLON_START
				if DEBUG == 'Y': print iLON,kLON,jLON
				TEMP_AREA = WET_AREA[iDATA,:,jLON]
				WET_AREA_LON[iDATA,iLON] = WET_AREA_LON[iDATA,iLON]+TEMP_AREA.sum()
				if DEBUG == 'Y': print TEMP_AREA[TEMP_AREA > 0],WET_AREA_LON[iDATA,iLON]
#
		if DEBUG == 'Y':
			print(LON_ZONAL)
			print(WET_AREA_LON[iDATA,:])
#
# Return to calling routine
#
	return WET_AREA_LON
#
