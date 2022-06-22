#
# Python module to regird dataset
#
# Garry Hayman
# Centre for Ecology and Hydrology
# December 2011
#
import numpy as np
from fractions import gcd
import statistics
#
def data_regrid_factor(DEBUG,NFILES,SRESOL_ALL,LONG_START,LAT_START,LONG_END,LAT_END):
#
	RESOL_LONG     = []
	RESOL_LAT      = []
	NLAT           = []
	NLONG          = []
#
	for iFILE in range(NFILES):
#
		SRESOL            = SRESOL_ALL[iFILE]
#
		if SRESOL=='UM':
			RESOL_LONG.append(1.875)
			RESOL_LAT.append(1.250)
#			NLAT.append(int((LAT_END-LAT_START)/RESOL_LAT[iFILE])+1)
			NLAT.append(int((LAT_END-LAT_START)/RESOL_LAT[iFILE]))
			NLONG.append(int((LONG_END-LONG_START)/RESOL_LONG[iFILE]))
		elif SRESOL=='MACC' or SRESOL == 'STOCHEM':
			RESOL_LONG.append(3.75)
			RESOL_LAT.append(2.50)
			NLAT.append(int((LAT_END-LAT_START)/RESOL_LAT[iFILE]))
			NLONG.append(int((LONG_END-LONG_START)/RESOL_LONG[iFILE]))
		elif SRESOL=='pTOMCAT_T42':
			RESOL_LONG.append(2.8125)
			RESOL_LAT.append(2.8125)
			NLAT.append(int((LAT_END-LAT_START)/RESOL_LAT[iFILE]))
			NLONG.append(int((LONG_END-LONG_START)/RESOL_LONG[iFILE]))
		else:
			RESOL_LONG.append(float(SRESOL))
			RESOL_LAT.append(float(SRESOL))
			NLAT.append(int((LAT_END-LAT_START)/RESOL_LAT[iFILE]))
			NLONG.append(int((LONG_END-LONG_START)/RESOL_LONG[iFILE]))
#
	RESOL_LONG_HGH_0,INDEX_LON = min((RESOL_LONG_HGH_0,INDEX_LON) for (INDEX_LON,RESOL_LONG_HGH_0) in enumerate(RESOL_LONG))
	RESOL_LONG_HGH    = gcd(RESOL_LONG[0],RESOL_LONG[1])
#
	RESOL_LAT_HGH_0,INDEX_LAT  = min((RESOL_LAT_HGH_0,INDEX_LAT) for (INDEX_LAT,RESOL_LAT_HGH_0) in enumerate(RESOL_LAT))
	RESOL_LAT_HGH     = gcd(RESOL_LAT[0],RESOL_LAT[1])
#
	print RESOL_LONG,RESOL_LONG_HGH_0,RESOL_LONG_HGH
	print RESOL_LAT,RESOL_LAT_HGH_0,RESOL_LAT_HGH
#
	for iFILE in range(NFILES):
		print RESOL_LONG[iFILE]/RESOL_LONG_HGH,RESOL_LAT[iFILE]/RESOL_LAT_HGH
#
	if INDEX_LON == 0:
		NLONG_HGH         = int((LONG_END-LONG_START)/RESOL_LONG_HGH)
		NLONG_IN          = NLONG[INDEX_LON]
		NLONG_OUT         = NLONG[INDEX_LON+1]
	elif RESOL_LONG[0] >= RESOL_LONG[1]:
		NLONG_HGH         = int((LONG_END-LONG_START)/RESOL_LONG_HGH)
		NLONG_IN          = NLONG[INDEX_LON-1]
		NLONG_OUT         = NLONG[INDEX_LON]
	else:
		NLONG_HGH         = int((LONG_END-LONG_START)/RESOL_LONG_HGH)
		NLONG_IN          = NLONG[INDEX_LON]
		NLONG_OUT         = NLONG[INDEX_LON-1]
#
	if INDEX_LAT == 0:
		NLAT_HGH          = int((LAT_END-LAT_START)/RESOL_LAT_HGH)
		NLAT_IN           = NLAT[INDEX_LAT]
		NLAT_OUT          = NLAT[INDEX_LAT+1]
	elif RESOL_LAT[0] >= RESOL_LAT[1]:
		NLAT_HGH          = int((LAT_END-LAT_START)/RESOL_LAT_HGH)
		NLAT_IN           = NLAT[INDEX_LAT-1]
		NLAT_OUT          = NLAT[INDEX_LAT]
	else:
		NLAT_HGH          = int((LAT_END-LAT_START)/RESOL_LAT_HGH)
		NLAT_IN           = NLAT[INDEX_LAT]
		NLAT_OUT          = NLAT[INDEX_LAT-1]
#
	INDEX             = INDEX_LON
#
	print INDEX,INDEX_LON,INDEX_LAT,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
#
	return INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
#
def data_map_difference(DEBUG,INDEX,DATA_MAP,NUM_IN,MIN_EMISS,MISS_DATA, \
	NLONG,NLAT,NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT):
#
	if INDEX == -1:
		DATA_0      = DATA_MAP[0]
		DATA_0[(DATA_0 <= 0.0) | (DATA_0 == float('nan'))] = MIN_EMISS
		DATA_1      = DATA_MAP[1]
		DATA_1[(DATA_1 <= 0.0) | (DATA_1 == float('nan'))] = MIN_EMISS
	elif INDEX == 0:
		print INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
		DATA_0=data_regrid_map(DEBUG,DATA_MAP[0],NUM_IN,MISS_DATA, \
			NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
		DATA_0[(DATA_0 <= 0.0) | (DATA_0 == float('nan'))] = MIN_EMISS
		DATA_1      = DATA_MAP[1]
		DATA_1[(DATA_1 <= 0.0) | (DATA_1 == float('nan'))] = MIN_EMISS
	else:
		print INDEX,NLONG_HGH,NLAT_HGH,NLONG_IN,NLAT_IN,NLONG_OUT,NLAT_OUT
		DATA_1=data_regrid_map(DEBUG,DATA_MAP[1],NUM_IN,MISS_DATA, \
			NLONG,NLAT,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT)
#
		DATA_1[(DATA_1 <= 0.0) | (DATA_1 == float('nan'))] = MIN_EMISS
		DATA_0      = DATA_MAP[0]
		DATA_0[(DATA_0 <= 0.0) | (DATA_0 == float('nan'))] = MIN_EMISS
#
	print INDEX,DATA_0.shape,DATA_1.shape,DATA_MAP[0].shape,DATA_MAP[1].shape
	DATA_D      = DATA_0-DATA_1
	print DATA_0.min(),DATA_0.max(),DATA_1.min(),DATA_1.max(),DATA_D.min(),DATA_D.max()
#
# Reset values for plotting
#
	DATA_0[DATA_0 <= MIN_EMISS] = float('nan')
	DATA_1[DATA_1 <= MIN_EMISS] = float('nan')
	DATA_D[(DATA_D >= -MIN_EMISS) & (DATA_D <= MIN_EMISS)] = float('nan')
	print DATA_0.min(),DATA_0.max(),DATA_1.min(),DATA_1.max(),DATA_D.min(),DATA_D.max()
#
	return DATA_0,DATA_1,DATA_D
#
def data_regrid_map(DEBUG,DATA_IN,NUM_IN,MISS_DATA, \
	NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT):
#
	DATA_INT    = np.zeros((NLAT_HGH,NLONG_HGH))
	NUM_INT     = np.zeros((NLAT_HGH,NLONG_HGH))
#
	DATA_OUT    = np.zeros((NLAT_OUT,NLONG_OUT))
	DATA_OUT[:,:] = MISS_DATA
#
# Regrid to regular output array
# - regrid to intermediate high resolution array with common factor to both input and output array sizes
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_IN)
	NPTS_LAT    = int(NLAT_HGH/NLAT_IN)
#
	print NLONG_IN,NLONG_HGH,NPTS_LONG,NLAT_IN,NLAT_HGH,NPTS_LAT
#
	for iLAT in range(NLAT_IN):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_IN):
#
			iLONG_START = iLONG*NPTS_LONG
			if DEBUG == 'Y': print(iLAT+1,iLONG+1,iLAT_START,iLONG_START)
#
			for jLAT in range(NPTS_LAT):
				for jLONG in range(NPTS_LONG):
					DATA_INT[jLAT+iLAT_START,jLONG+iLONG_START] = DATA_IN[iLAT,iLONG]
					NUM_INT[jLAT+iLAT_START,jLONG+iLONG_START]  = NUM_IN[iLAT,iLONG]
#
# - regrid to output array
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_OUT)
	NPTS_LAT    = int(NLAT_HGH/NLAT_OUT)
#
	if DEBUG == 'Y':
		print(NLONG_OUT,NLONG_HGH,NPTS_LONG,NLAT_OUT,NLAT_HGH,NPTS_LAT)
		print(NUM_IN.min(),NUM_IN.max(),DATA_IN.min(),DATA_IN.max())
		print(NUM_INT.min(),NUM_INT.max(),DATA_INT.min(),DATA_INT.max())
		print(DATA_IN[np.where(DATA_IN >= 0)].min(),DATA_IN[np.where(DATA_IN >= 0)].max())
		print(DATA_INT[np.where(DATA_INT >= 0)].min(),DATA_INT[np.where(DATA_INT >= 0)].max())
#
	for iLAT in range(NLAT_OUT):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_OUT):
#
			iLONG_START = iLONG*NPTS_LONG
#
			COUNT       = 0.0
			SUM         = 0.0
#
			for jLAT in range(NPTS_LAT):
				for jLONG in range(NPTS_LONG):
#
					DATA_LOCAL  = DATA_INT[jLAT+iLAT_START,jLONG+iLONG_START]
					NUM_LOCAL   = NUM_INT[jLAT+iLAT_START,jLONG+iLONG_START]
#
					if DATA_LOCAL >0:
						SUM     = SUM   + DATA_LOCAL*NUM_LOCAL
						COUNT   = COUNT + NUM_LOCAL
#
			if COUNT > 0:
				DATA_OUT[iLAT,iLONG] = float('%.2f' % (SUM/COUNT))
#
	return DATA_OUT
#
def data_regrid(DEBUG,DATA_IN,NUM_IN,NTIMES,MISS_DATA, \
	NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT):
#
	DATA_INT    = np.zeros((NTIMES,NLAT_HGH,NLONG_HGH))
	NUM_INT     = np.zeros((NTIMES,NLAT_HGH,NLONG_HGH))
#
	DATA_OUT    = np.zeros((NTIMES,NLAT_OUT,NLONG_OUT))
	DATA_OUT[:,:,:] = MISS_DATA
#
# Regrid to regular output array
# - regrid to intermediate high resolution array with common factor to both input and output array sizes
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_IN)
	NPTS_LAT    = int(NLAT_HGH/NLAT_IN)
#
	print(NTIMES,NLONG_IN,NLONG_HGH,NPTS_LONG,NLAT_IN,NLAT_HGH,NPTS_LAT)
#
	for iLAT in range(NLAT_IN):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_IN):
#
			iLONG_START = iLONG*NPTS_LONG
			if DEBUG == 'Y': print(iLAT+1,iLONG+1,iLAT_START,iLONG_START)
#
			for jLAT in range(NPTS_LAT):
				for jLONG in range(NPTS_LONG):
					DATA_INT[:,jLAT+iLAT_START,jLONG+iLONG_START] = DATA_IN[:,iLAT,iLONG]
					NUM_INT[:,jLAT+iLAT_START,jLONG+iLONG_START]  = NUM_IN[:,iLAT,iLONG]
#
# - regrid to output array
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_OUT)
	NPTS_LAT    = int(NLAT_HGH/NLAT_OUT)
#
	if DEBUG == 'Y':
		print(NTIMES,NLONG_OUT,NLONG_HGH,NPTS_LONG,NLAT_OUT,NLAT_HGH,NPTS_LAT)
		print(NUM_IN.min(),NUM_IN.max(),DATA_IN.min(),DATA_IN.max())
		print(NUM_INT.min(),NUM_INT.max(),DATA_INT.min(),DATA_INT.max())
		print(DATA_IN[np.where(DATA_IN >= 0)].min(),DATA_IN[np.where(DATA_IN >= 0)].max())
		print(DATA_INT[np.where(DATA_INT >= 0)].min(),DATA_INT[np.where(DATA_INT >= 0)].max())
#
	for iLAT in range(NLAT_OUT):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_OUT):
#
			iLONG_START = iLONG*NPTS_LONG
#
			for iTIME in range(NTIMES):
#
				COUNT       = 0.0
				SUM         = 0.0
#
				for jLAT in range(NPTS_LAT):
					for jLONG in range(NPTS_LONG):
#
						DATA_LOCAL  = DATA_INT[iTIME,jLAT+iLAT_START,jLONG+iLONG_START]
						NUM_LOCAL   = NUM_INT[iTIME,jLAT+iLAT_START,jLONG+iLONG_START]
#
						if DATA_LOCAL >0:
							SUM     = SUM   + DATA_LOCAL*NUM_LOCAL
							COUNT   = COUNT + NUM_LOCAL
#
				if COUNT > 0:
					DATA_OUT[iTIME,iLAT,iLONG] = float('%.2f' % (SUM/COUNT))
#
	return DATA_OUT
#
# Recoding of subroutine aggregate from FORTRAN program Methane_Analysis_process_v10c.f90
#
def data_Sciamachy_regrid(DEBUG,DATA_IN,NUM_IN,NTIMES,MISS_DATA, \
	NLONG_OUT,NLAT_OUT,NLONG_IN,NLAT_IN):
#
# Simply copy if output resolution is the same (or smaller) than the
# input resolution
#
	if NLONG_OUT == NLONG_IN and NLAT_OUT == NLAT_IN:
#
		DATA_OUT = DATA_IN
#
	else:
#
		DATA_OUT        = np.zeros((NTIMES,NLAT_OUT,NLONG_OUT))
		DATA_OUT[:,:,:] = MISS_DATA
#
# Derive weighted average over relevant grid squares.  Weights from
# number of measurements in grid square
#
		NPTS_LONG = int(NLONG_IN/NLONG_OUT)
		NPTS_LAT  = int(NLAT_IN/NLAT_OUT)
#
		print(NTIMES,NLONG_OUT,NLONG_IN,NPTS_LONG,NLAT_OUT,NLAT_IN,NPTS_LAT)
#
# Loop over time
#
		for iTIME in range(NTIMES):
#
# Loop over output longitude
#
			for iLONG_OUT in range(NLONG_OUT):
#
# Loop over output longitude
#
				for iLAT_OUT in range(NLAT_OUT):
#
					iLONG_ST  = iLONG_OUT*NPTS_LONG
					iLONG_END = (iLONG_OUT+1)*NPTS_LONG-1
					iLAT_ST   = iLAT_OUT*NPTS_LAT
					iLAT_END  = (iLAT_OUT+1)*NPTS_LAT-1
#
					SUM    = 0.0
					COUNT  = 0.0
#
					for iLONG_IN in range(iLONG_ST,iLONG_END+1,1):
						for iLAT_IN in range(iLAT_ST,iLAT_END+1,1):
#
							CH4_COL  = DATA_IN[iTIME,iLAT_IN,iLONG_IN]
							CH4_NUM  = NUM_IN[iTIME,iLAT_IN,iLONG_IN]
#
							if DEBUG == 'Y':
								print(iLONG_OUT+1,iLAT_OUT+1,iLONG_IN+1,iLAT_IN+1, \
									CH4_COL,CH4_NUM)
#
# Ignore missing data points (-999900000).  Initially assumed that these would have
# NUM_IN(iTIMEiLAT_IN,iLONG_IN) = 0.  However, there are instances
# where this is not the case (e.g., December 2005).  Only use positive values of CH4
# column.
#
							if CH4_COL > 0.0:
#
# Sum CH4 column measurement, weighted by number of points per grid cell
#
								SUM      = SUM   + CH4_COL*CH4_NUM
#
# Sum number of points per grid cell 
#
								COUNT    = COUNT + CH4_NUM
#
							 	if DEBUG == 'Y':
									print(iLONG_OUT+1,iLAT_OUT+1,iLONG_IN+1,iLAT_IN+1, \
									CH4_NUM,COUNT,CH4_COL,SUM)	
#
#  Calculate weighted average column CH4 concentration for output box
#
								if COUNT > 0.0:
			   						DATA_OUT[iTIME,iLAT_OUT,iLONG_OUT] = float('%.2f' % (SUM/COUNT))
#
					if DEBUG == 'Y':
						print(iLONG_OUT+1,iLAT_OUT+1,DATA_OUT[iTIME,iLAT_OUT,iLONG_OUT])
#
# Return to calling routine
#
	return DATA_OUT
#
def data_regrid_TRANSCOM(DEBUG,DATA_IN,NUM_IN,NTIMES,MISS_DATA, \
	NLONG_IN,NLAT_IN,NLONG_HGH,NLAT_HGH,NLONG_OUT,NLAT_OUT):
#
#
	DATA_INT    = np.zeros((NTIMES,NLAT_HGH,NLONG_HGH))
	NUM_INT     = np.zeros((NTIMES,NLAT_HGH,NLONG_HGH))
#
	DATA_OUT    = np.zeros((NTIMES,NLAT_OUT,NLONG_OUT))
	DATA_OUT[:,:,:] = MISS_DATA
#
# Regrid to regular output array
# - regrid to intermediate high resolution array with common factor to both input and output array sizes
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_IN)
	NPTS_LAT    = int(NLAT_HGH/NLAT_IN)
#
	print(NTIMES,NLONG_IN,NLONG_HGH,NPTS_LONG,NLAT_IN,NLAT_HGH,NPTS_LAT)
#
	for iLAT in range(NLAT_IN):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_IN):
#
			iLONG_START = iLONG*NPTS_LONG
			if DEBUG == 'Y': print(iLAT+1,iLONG+1,iLAT_START,iLONG_START)
#
			for jLAT in range(NPTS_LAT):
				for jLONG in range(NPTS_LONG):
					DATA_INT[:,jLAT+iLAT_START,jLONG+iLONG_START] = DATA_IN[:,iLAT,iLONG]
					NUM_INT[:,jLAT+iLAT_START,jLONG+iLONG_START]  = NUM_IN[:,iLAT,iLONG]
#
# - regrid to output array
#
	NPTS_LONG   = int(NLONG_HGH/NLONG_OUT)
	NPTS_LAT    = int(NLAT_HGH/NLAT_OUT)
#
	if DEBUG == 'Y':
		print(NTIMES,NLONG_OUT,NLONG_HGH,NPTS_LONG,NLAT_OUT,NLAT_HGH,NPTS_LAT)
		print(NUM_IN.min(),NUM_IN.max(),DATA_IN.min(),DATA_IN.max())
		print(NUM_INT.min(),NUM_INT.max(),DATA_INT.min(),DATA_INT.max())
		print(DATA_IN[np.where(DATA_IN >= 0)].min(),DATA_IN[np.where(DATA_IN >= 0)].max())
		print(DATA_INT[np.where(DATA_INT >= 0)].min(),DATA_INT[np.where(DATA_INT >= 0)].max())
#
	MAX_VALUE   = DATA_IN.max()
#
	for iLAT in range(NLAT_OUT):
#
		iLAT_START  = iLAT*NPTS_LAT
#
		for iLONG in range(NLONG_OUT):
#
			iLONG_START = iLONG*NPTS_LONG
#
			for iTIME in range(NTIMES):
#
				DATA_LOCAL = []
				NUM_LOCAL  = []
#
				for jLAT in range(NPTS_LAT):
					for jLONG in range(NPTS_LONG):
						DATA_LOCAL.append(int(DATA_INT[iTIME,jLAT+iLAT_START,jLONG+iLONG_START]))
						NUM_LOCAL.append(int(NUM_INT[iTIME,jLAT+iLAT_START,jLONG+iLONG_START]))
#
				MODE = statistics.mode(DATA_LOCAL,NUM_LOCAL)
#
				DATA_OUT[iTIME,iLAT,iLONG] = float(MODE[0][0])
				if len(MODE) != 1:
					print(iLAT,iLONG,DATA_LOCAL,NUM_LOCAL,MODE,DATA_OUT[iTIME,iLAT,iLONG])
#
	return DATA_OUT
#
