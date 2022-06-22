#
# Python module to calculate regression parameters
#
# Garry Hayman
# Centre for Ecology and Hydrology
# April 2012
#
import math
import numpy as np
#
# Contains
#
def calculate_Regress_Coefficients(X,Y):
#
# Subroutine to calculate linear regression coefficients
#
	NPOINTS       = len(X)
#
	SUMX          = X.sum()
	SUMY          = Y.sum()
	SUMXY         = (X*Y).sum()
	SUMX2         = (X*X).sum()
#
	XMEAN         = SUMX/NPOINTS
	YMEAN         = SUMY/NPOINTS
#
	MSE           = ((X-Y)*(X-Y)).sum()
	RMSE          = math.sqrt(MSE/NPOINTS)
	CVRMSE        = RMSE/XMEAN
#
	DIFFX         = X-XMEAN
	DIFFY         = Y-YMEAN
	DIFFXY        = DIFFX*DIFFY
	DIFFX2        = DIFFX*DIFFX
	DIFFY2        = DIFFY*DIFFY
#
	SUMDXY        = DIFFXY.sum()
	SUMDX2        = DIFFX2.sum()
	SUMDY2        = DIFFY2.sum()
#
	SLOPE	      =  SUMDXY/SUMDX2
	INTERCEPT     = (SUMY*SUMX2-SUMX*SUMXY)/(NPOINTS*SUMX2-SUMX*SUMX)
	CORREL	      =  SUMDXY/math.sqrt(SUMDX2*SUMDY2)
	COEFF_REGRESS =  CORREL*CORREL 
	X_STDDEV      =  math.sqrt(SUMDX2/NPOINTS)
	Y_STDDEV      =  math.sqrt(SUMDY2/NPOINTS)
#
	DIFFEY        =  Y-(SLOPE*X+INTERCEPT)
	DIFFEY2       =  DIFFEY*DIFFEY
#
 	SUMEY         =  DIFFEY.sum()
	SUMEY2        =  DIFFEY2.sum()
#
	SUM_TOTAL     =  SUMDY2
	SUM_RESID     =  SUMEY2
	SUM_REGRESS   =  SUM_TOTAL-SUM_RESID
	ERROR	      =  math.sqrt(SUM_RESID/(NPOINTS-2))
	ERROR_SLOPE   =  ERROR/math.sqrt(SUMDX2)
	ERROR_INTER   =  ERROR*math.sqrt(SUMX2/(NPOINTS*SUMDX2))
#
	TAYLOR_R      =  Y_STDDEV/X_STDDEV
	TAYLOR_A      =  math.acos(CORREL)
#
	if SLOPE < 0.0:
		CORREL = -CORREL
#
	return  SLOPE,INTERCEPT,CORREL,COEFF_REGRESS,X_STDDEV,Y_STDDEV, \
		ERROR,ERROR_SLOPE,ERROR_INTER,TAYLOR_R,TAYLOR_A,RMSE,CVRMSE
#
def calculate_Model_Performance(X,Y,REL_ERROR):
#
# Subroutine to calculate the model performance parameters
# Mean bias, factor of 2 and index of agreement
#
# Based on Appendix E in WMO-GAW Report 181
#
	NPOINTS       = len(X)
#
	XMEAN         = np.mean(X)
	YMEAN         = np.mean(Y)
#
	BIAS          = (YMEAN-XMEAN)
	MNB           = ((Y-X)/X).sum()/NPOINTS
#
	FAC2          = float(len(X[((Y/X) <= 2.0) & ((Y/X) >= 0.5)]))/float(NPOINTS)
	HIT_RATE      = float(len(X[np.abs((Y-X)/X) <= REL_ERROR]))/float(NPOINTS)
#	print(REL_ERROR,HIT_RATE,float(len(X[np.abs((Y-X)/X) <= REL_ERROR])),NPOINTS,X,Y,np.abs((Y-X)/X))
#
	MSE           = ((Y-X)*(Y-X)).sum()
#
# Garry Hayman (December 2013)
# Error in WMO report for derivation of IOA (range is 0-1 but it was giving negative values on occasion)
# See http://www.rforge.net/doc/packages/hydroGOF/d.html for correct definition
#
#	DEN           = np.abs(X-XMEAN)+np.abs(Y-YMEAN)
	DEN           = np.abs(X-XMEAN)+np.abs(Y-XMEAN)
	IOA           = 1.0-MSE/(DEN*DEN).sum()
#
	return  NPOINTS,XMEAN,YMEAN,BIAS,MNB,FAC2,IOA,HIT_RATE
#
def calculate_Model_Performance_IOA2(X,Y):
#
# Subroutine to calculate revised index of agreement
#
# Revised IOA: Willmott et al., International Journal of Climatology, 2011
#
	XMEAN         = np.mean(X)
	YMEAN         = np.mean(Y)
#
	ABS1          = (np.abs(X-XMEAN)).sum()*2.0
	ABS2          = (np.abs(X-Y)).sum()
#
	if ABS2 <= ABS1:
		IOA           = 1.0-ABS2/ABS1
	else:
		IOA           = ABS1/ABS2-1.0

	return  IOA
#
def mode(DATA,NUM):
#
	FREQUENCIES = {}
	ICOUNT      = 0
#
	for X in DATA:
#
		if (X in FREQUENCIES):
			FREQUENCIES[X] += NUM[ICOUNT]
		else:
			FREQUENCIES[X]  = NUM[ICOUNT]
#
		ICOUNT     += 1
#
	MODE        = max(FREQUENCIES.values())
	MODE_ALL    = [(X, MODE) for X in FREQUENCIES if (MODE == FREQUENCIES[X])]
#
	return MODE_ALL
#
def calculate_all(NCRITICAL,X_ALL,Y_ALL,REL_ERROR,STAT_SELECT,MISS_DATA,DEBUG):
#
	DEBUG       = 'Y'
	X_ALL       = np.array(X_ALL)
	Y_ALL       = np.array(Y_ALL)
#
	NRUNS       = X_ALL.shape[0]
#
# Define arrays
#
	SLOPE         = np.zeros(NRUNS)
	INTERCEPT     = np.zeros(NRUNS)
	CORREL        = np.zeros(NRUNS)
	COEFF_REGRESS = np.zeros(NRUNS)
	X_STDDEV      = np.zeros(NRUNS)
	Y_STDDEV      = np.zeros(NRUNS)
	ERROR         = np.zeros(NRUNS)
	ERROR_SLOPE   = np.zeros(NRUNS)
	ERROR_INTER   = np.zeros(NRUNS)
	TAYLOR_R      = np.zeros(NRUNS)
	TAYLOR_A      = np.zeros(NRUNS)
	RMSE          = np.zeros(NRUNS)
	CVRMSE        = np.zeros(NRUNS)
	NPOINTS       = np.zeros(NRUNS)
	X_MEAN        = np.zeros(NRUNS)
	Y_MEAN        = np.zeros(NRUNS)
	BIAS          = np.zeros(NRUNS)
	MNB           = np.zeros(NRUNS)
	FAC2          = np.zeros(NRUNS)
	IOA           = np.zeros(NRUNS)
	HIT_RATE      = np.zeros(NRUNS)
#
	for iRUN in range(NRUNS):
#
		X           = X_ALL[iRUN,:]
		Y           = Y_ALL[iRUN,:]
#
# Remove any missing data
#
		XREGRESS    = X[(~np.isnan(X)) & (~np.isnan(Y))]
		YREGRESS    = Y[(~np.isnan(X)) & (~np.isnan(Y))]
#
		if DEBUG == 'Y':
			print iRUN,NRUNS,len(XREGRESS)
			print X
			print Y
			print XREGRESS
			print YREGRESS
#
		if len(XREGRESS) >= NCRITICAL:
#
			SLOPE[iRUN],INTERCEPT[iRUN],CORREL[iRUN],COEFF_REGRESS[iRUN], \
				X_STDDEV[iRUN],Y_STDDEV[iRUN],ERROR[iRUN],ERROR_SLOPE[iRUN], \
				ERROR_INTER[iRUN],TAYLOR_R[iRUN],TAYLOR_A[iRUN],RMSE[iRUN],CVRMSE[iRUN] = \
				calculate_Regress_Coefficients(XREGRESS,YREGRESS)
#
			NPOINTS[iRUN],X_MEAN[iRUN],Y_MEAN[iRUN],BIAS[iRUN],MNB[iRUN],FAC2[iRUN],IOA[iRUN],HIT_RATE[iRUN]= \
				calculate_Model_Performance(XREGRESS,YREGRESS,REL_ERROR)
#
			if STAT_SELECT == 'IOA2':
				IOA[iRUN] = calculate_Model_Performance_IOA2(X,Y) 
#
		else:
#
			SLOPE[iRUN]         = MISS_DATA
			INTERCEPT[iRUN]     = MISS_DATA
			CORREL[iRUN]        = MISS_DATA
			COEFF_REGRESS[iRUN] = MISS_DATA
			X_STDDEV[iRUN]      = MISS_DATA
			Y_STDDEV[iRUN]      = MISS_DATA
			ERROR[iRUN]         = MISS_DATA
			ERROR_SLOPE[iRUN]   = MISS_DATA
			ERROR_INTER[iRUN]   = MISS_DATA
			TAYLOR_R[iRUN]      = MISS_DATA
			TAYLOR_A[iRUN]      = MISS_DATA
			RMSE[iRUN]          = MISS_DATA
			CVRMSE[iRUN]        = MISS_DATA
			NPOINTS[iRUN]       = MISS_DATA
			X_MEAN[iRUN]        = MISS_DATA
			Y_MEAN[iRUN]        = MISS_DATA
			BIAS[iRUN]          = MISS_DATA
			MNB[iRUN]           = MISS_DATA
			FAC2[iRUN]          = MISS_DATA
			IOA[iRUN]           = MISS_DATA
			HIT_RATE[iRUN]      = MISS_DATA
#	
	return [ NPOINTS,SLOPE,INTERCEPT,CORREL,COEFF_REGRESS,X_STDDEV,Y_STDDEV,ERROR,ERROR_SLOPE, \
		ERROR_INTER,TAYLOR_R,TAYLOR_A,RMSE,CVRMSE,X_MEAN,Y_MEAN,BIAS,MNB,FAC2,IOA,HIT_RATE ]
