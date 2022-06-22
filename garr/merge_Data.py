#
# Python module to generate a single dataset from multiple datasets
#
# Garry Hayman
# Centre for Ecology and Hydrology
# April 2012
#
import numpy as np
from numpy import arange
#
# Contains
#
def merge_Data(X1,Y1,X2,Y2,MISS_DATA,DEBUG):
#
#
	if DEBUG == 'Y':
		print('X1 = ',X1)
		print('Y1 = ',Y1)
		print('X2 = ',X2)
		print('Y2 = ',Y2)
#
	INTERVAL   = 1.0/24.0
#
	NDATA1     = len(X1)
	NDATA2     = len(X2)
#
	if NDATA1 > NDATA2:
		NDATA      = NDATA2
		X1NEW      = X2
		Y1NEW      = Y2
		X2NEW      = X1
		Y2NEW      = Y1
	else:
		NDATA      = NDATA1
		X1NEW      = X1
		Y1NEW      = Y1
		X2NEW      = X2
		Y2NEW      = Y2
#
	X1MERGE    = []
	Y1MERGE    = []
	X2MERGE    = []
	Y2MERGE    = []
#
	for iDATA in range(NDATA):
#
		INDEX,     = np.where(np.abs(X1NEW[iDATA]-X2NEW)<INTERVAL)
#
		if len(INDEX) > 0 and int(Y1NEW[iDATA]) != int(MISS_DATA) and int(Y2NEW[INDEX]) != int(MISS_DATA):
			X1MERGE.append(X1NEW[iDATA])
			Y1MERGE.append(Y1NEW[iDATA])
			X2MERGE.append(np.asscalar(X2NEW[INDEX]))
			Y2MERGE.append(np.asscalar(Y2NEW[INDEX]))
#
	if DEBUG == 'Y':
		print('X1MERGE = ',X1MERGE)
		print('Y1MERGE = ',Y1MERGE)
		print('X2MERGE = ',X2MERGE)
		print('Y2MERGE = ',Y2MERGE)
#
	NDATA     = len(X1MERGE)
#
	if NDATA1 > NDATA2:
		X1FINAL    = np.array(X2MERGE)
		Y1FINAL    = np.array(Y2MERGE)
		X2FINAL    = np.array(X1MERGE)
		Y2FINAL    = np.array(Y1MERGE)
	else:
		X1FINAL    = np.array(X1MERGE)
		Y1FINAL    = np.array(Y1MERGE)
		X2FINAL    = np.array(X2MERGE)
		Y2FINAL    = np.array(Y2MERGE)
#
	return NDATA,X1FINAL,Y1FINAL,X2FINAL,Y2FINAL
#
def insert_Data(X,Y,DEBUG):
#
	NPOINTS  = len(X)
#
	XNEW     = []
	YNEW     = []
#
	INDICES  = [(i,j) for i in range(NPOINTS) for j in range(i+1,NPOINTS)]
#
	MIN_SEP  = 999999999
#
	for i,j in INDICES: MIN_SEP   = min(MIN_SEP,X[j]-X[i])
#
	for i in range(NPOINTS-1):
#
		XNEW.append(X[i])
		YNEW.append(Y[i])
#
		if X[i+1]-X[i]>2*MIN_SEP:
			XNEW.append(float('nan'))
			YNEW.append(float('nan'))
#
	if len(X) > 0:
		XNEW.append(X[-1])
		YNEW.append(Y[-1])
#
	XNEW     = np.array(XNEW)
	YNEW     = np.array(YNEW)
#
	if DEBUG == 'Y':
		print X
		print Y
		print XNEW
		print YNEW
#
	NPOINTS  = len(XNEW)
#
	return NPOINTS,XNEW,YNEW
#
def insert_Data_multi(NDATA,X,Y,DEBUG):
#
	XALL     = []
	YALL     = []
	NPTS_ALL = np.zeros(NDATA)
#
	for iDATA in range(NDATA):
#
		NPOINTS  = len(X[iDATA,:])
#
		XNEW     = []
		YNEW     = []
#
		INDICES  = [(i,j) for i in range(NPOINTS) for j in range(i+1,NPOINTS)]
#
		MIN_SEP  = 999999999
#
		for i,j in INDICES: MIN_SEP   = min(MIN_SEP,X[iDATA,j]-X[iDATA,i])
#
		for i in range(NPOINTS-1):
#
			XNEW.append(X[iDATA,i])
			YNEW.append(Y[iDATA,i])
#
			if X[iDATA,i+1]-X[iDATA,i]>2*MIN_SEP:
				XNEW.append(float('nan'))
				YNEW.append(float('nan'))
#
		XNEW.append(X[iDATA,-1])
		YNEW.append(Y[iDATA,-1])
#
		XALL.append(XNEW)
		YALL.append(YNEW)
#
		NPTS_ALL[iDATA] = len(XNEW)
#
	if DEBUG == 'Y': print X.shape,NPTS_ALL,NPTS_ALL.max()
	XOUT     = np.zeros((NDATA,NPTS_ALL.max()))
	YOUT     = np.zeros((NDATA,NPTS_ALL.max()))
#
	XOUT[:,:] = float('nan')
	YOUT[:,:] = float('nan')
#
	for iDATA in range(NDATA):
		XOUT[iDATA,0:NPTS_ALL[iDATA]] = np.array(XALL[iDATA][:])
		YOUT[iDATA,0:NPTS_ALL[iDATA]] = np.array(YALL[iDATA][:])
#
	if DEBUG == 'Y':
		print X
		print Y
		print XOUT
		print YOUT
#
	return XOUT,YOUT
#
def climatology_single(X,Y,MISS_DATA,DEBUG):
#
# Assumed to be presented as annual dataset
#
#
	DEBUG      = 'Y'
	NMONTHS    = 12
#
	XCLIM      = np.zeros(NMONTHS)
	YCLIM      = np.zeros(NMONTHS)
#
	for iMONTH in range(NMONTHS):
#
		TEMP            = Y[iMONTH::12]
		YCLIM[iMONTH]   = np.average(TEMP[~np.isnan(TEMP)])
#
		if DEBUG == 'Y':
			print iMONTH,YCLIM[iMONTH],Y[iMONTH::12]
#
	XCLIM[:] = 0.5+arange(NMONTHS)
#
	return XCLIM,YCLIM
#
def climatology_multi(X,Y,MISS_DATA,DEBUG):
#
	X1         = X[0,:]
	Y1         = Y[0,:]
	X2         = X[1,:]
	Y2         = Y[1,:]
#
	INTERVAL   = 1.0/24.0
#
	NDATA      = len(X2[~np.isnan(X2)])
	NMONTHS    = 12
#
	Y1CLIM     = []
	Y2CLIM     = []
	XCLIM      = np.zeros((2,NMONTHS))
	YCLIM      = np.zeros((2,NMONTHS))
#
	for iDATA in range(NDATA):
#
		INDEX,     = np.where(np.abs(X2[iDATA]-X1)<INTERVAL)
		if DEBUG == 'Y': print iDATA,INDEX,X2[iDATA],X1[INDEX]
#
		Y2CLIM.append(Y2[iDATA])
		if len(INDEX) > 0:
			Y1CLIM.append(np.asscalar(Y1[INDEX]))
		else:
			Y1CLIM.append(float('nan'))
#
	Y1CLIM     = np.array(Y1CLIM)
	Y2CLIM     = np.array(Y2CLIM)
#
	if DEBUG == 'Y':
		print 'Y1CLIM = ',len(Y1CLIM),Y1CLIM
		print 'Y2CLIM = ',len(Y2CLIM),Y2CLIM
#
	for iMONTH in range(NMONTHS):
#
		TEMP            = Y1CLIM[iMONTH::12]
		YCLIM[0,iMONTH] = np.average(TEMP[~np.isnan(TEMP)])
		TEMP            = Y2CLIM[iMONTH::12]
		YCLIM[1,iMONTH] = np.average(TEMP[~np.isnan(TEMP)])
#
		if DEBUG == 'Y':
			print iMONTH,YCLIM[0,iMONTH],Y1CLIM[iMONTH::12]
			print iMONTH,YCLIM[1,iMONTH],Y2CLIM[iMONTH::12]
#
	XCLIM[0,:] = 0.5+arange(NMONTHS)
	XCLIM[1,:] = 0.5+arange(NMONTHS)
#
	return XCLIM,YCLIM
#
def get_increment(MIN,MAX):
#
# Derive increment for figure axes given min and max
#
	RANGE      = MAX-MIN
#
	if RANGE >   400:
		INC        = 100
	elif RANGE > 200:
		INC        =  50
	elif RANGE > 100:
		INC        =  20
	elif RANGE >  40:
		INC        =  20
	elif RANGE >  20:
		INC        =   5
	elif RANGE >   6:
		INC        =   2
	else:
		INC        =   1
#
	NTICKS      = int((MAX-MIN)/INC)+1
	TICKS       = MIN+INC*arange(NTICKS)
#	print MIN,MAX,RANGE,INC,NTICKS,TICKS
#
	return INC,NTICKS,TICKS
#
def climatology_single_regress(X,Y,MISS_DATA,DEBUG):
#
# Assumed to be presented as annual dataset
#
	NMONTHS    = 12
#
	XCLIM      = np.zeros(NMONTHS)
	YCLIM      = np.zeros(NMONTHS)
#
	for iMONTH in range(NMONTHS):
#
		TEMP            = X[iMONTH::12]
		XCLIM[iMONTH]   = np.average(TEMP[~np.isnan(TEMP)])
#
		TEMP            = Y[iMONTH::12]
		YCLIM[iMONTH]   = np.average(TEMP[~np.isnan(TEMP)])
#
		if DEBUG == 'Y':
			print iMONTH,XCLIM[iMONTH],X[iMONTH::12]
			print iMONTH,YCLIM[iMONTH],Y[iMONTH::12]
#
	return XCLIM,YCLIM
#
