# Python script to analyse output from a JULES run
#
# python JULES_map_netCDF_v3.py N Y 0.5 225 GRADS
#
import os
import numpy as np
from numpy import arange,dtype
#
import data_netCDF
#
# ##################################################################################################
# (1) JULES offline - get variables in dataset
# ##################################################################################################
#
def get_JULES_variables(FILE_CDF_IN,JULES_VARINFO):
#
	VAR_NAMES      = data_netCDF.data_netCDF_getVARNAMES(FILE_CDF_IN)
#
	for iVAR in range(len(VAR_NAMES)):
		VAR_NAME     = VAR_NAMES[iVAR]
		VAR_NAME_MOD   = VAR_NAME.lower()
#		print iVAR,VAR_NAMES[iVAR],JULES_VARINFO[VAR_NAME_MOD][1]
		TEXT         = ('%4d %-20s spin-up diagnostic: %s' \
			% (iVAR,VAR_NAMES[iVAR],JULES_VARINFO[VAR_NAME_MOD][1]))
		print(TEXT)
#
# Select variable
	print 'Input index for variable (select those with Y for valid spin-up diagnostic): '
	iVAR           = int(input())
	VAR_NAME       = VAR_NAMES[iVAR]
#
# Assign attributes
#
	NLEVELS        = 10
	iFLAG          =  1
	iSOIL          = -1
	VAR_NAME_MOD   = VAR_NAME.lower()
	PLOT_LABEL     = JULES_VARINFO[VAR_NAME_MOD][5]+JULES_VARINFO[VAR_NAME_MOD][6]
	PLOT_MIN       = float(JULES_VARINFO[VAR_NAME_MOD][3])
	PLOT_MAX       = float(JULES_VARINFO[VAR_NAME_MOD][4])
	PLOT_INC       = (PLOT_MAX-PLOT_MIN)/NLEVELS
	CLEVELS        = arange(PLOT_MIN,PLOT_MAX+PLOT_INC,PLOT_INC)
#
	print JULES_VARINFO[VAR_NAME_MOD][2]
#
	if 'z' in JULES_VARINFO[VAR_NAME_MOD][2]:
		print 'Input soil level (1-4): '
		iSOIL          = int(input())-1
#
	return VAR_NAME,VAR_NAME_MOD,PLOT_LABEL,PLOT_MIN,PLOT_MAX,PLOT_INC,CLEVELS,iSOIL,iFLAG
#
# ##################################################################################################
# (2.1) JULES offline - single global netCDF dataset
# ##################################################################################################
#
# ##################################################################################################
# (2.2) JULES offline - multi-annual global netCDF datasets
# ##################################################################################################
#
def get_JULES_data_global_multi(NYEARS,START_YEAR,FILE_PREFIX,FILE_SUFFIX,JULES_VARINFO, \
	MISS_DATA,iFLAG,DEBUG):
#
	COUNT       = 0
#
	for iYEAR in range(NYEARS):
#
		SYEAR       = ('%4s' % (START_YEAR+iYEAR))
		FILE_CDF_IN   = FILE_PREFIX+SYEAR+FILE_SUFFIX
		print ('%3d: ' % iYEAR)+FILE_CDF_IN
#
		if iFLAG == 0:
			VAR_NAME,VAR_NAME_MOD,PLOT_LABEL,PLOT_MIN, \
				PLOT_MAX,PLOT_INC,CLEVELS,iSOIL,iFLAG= \
				get_JULES_variables(FILE_CDF_IN,JULES_VARINFO)
#		
# Skip if file missing
#
		if os.path.exists(FILE_CDF_IN):
			if os.path.getsize(FILE_CDF_IN) > 0:
#
				DIMS,VAR_DATA = \
					data_netCDF.data_netCDF_array_var(FILE_CDF_IN,VAR_NAME)
				if DEBUG == 'Y':
					print VAR_DATA.shape,VAR_DATA.min(),VAR_DATA.max()
#
				if iYEAR == 0:
					NTIMES         = DIMS[0]*NYEARS
					DIMS_MOD       = [NTIMES]
					for i in range(1,len(DIMS)):
						DIMS_MOD.append(DIMS[i])
					VAR_DATA_ALL   = np.zeros(DIMS_MOD)
#
				NCOUNT        = DIMS[0]
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:]     = VAR_DATA[:,:]
				elif len(DIMS) == 3:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:]   = VAR_DATA[:,:,:]
				elif len(DIMS) == 4:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:,:] = VAR_DATA[:,:,:,:]
#
			else:
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:]     = MISS_DATA
				elif len(DIMS) == 3:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:]   = MISS_DATA
				elif len(DIMS) == 4:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:,:] = MISS_DATA
#
# Increment as NCOUNT timepoints
#
			COUNT += NCOUNT
#
	VAR_DATA_ALL = np.squeeze(VAR_DATA_ALL)
#
# Return to calling routine
#
	return VAR_DATA_ALL,VAR_NAME,VAR_NAME_MOD,PLOT_LABEL,PLOT_MIN,PLOT_MAX,PLOT_INC,CLEVELS,iSOIL
#
# ##################################################################################################
#
def get_JULES_data_global_multi2(NYEARS,START_YEAR,FILE_PREFIX,FILE_SUFFIX,VAR_NAME, \
	MISS_DATA,DEBUG):
#
	COUNT       = 0
#
	for iYEAR in range(NYEARS):
#
		SYEAR       = ('%4s' % (START_YEAR+iYEAR))
		FILE_CDF_IN   = FILE_PREFIX+SYEAR+FILE_SUFFIX
		print ('%3d: ' % iYEAR)+FILE_CDF_IN
#
# Skip if file missing
#
		if os.path.exists(FILE_CDF_IN):
			if os.path.getsize(FILE_CDF_IN) > 0:
#
				DIMS,VAR_DATA = \
					data_netCDF.data_netCDF_array_var(FILE_CDF_IN,VAR_NAME)
				if DEBUG == 'Y':
					print VAR_DATA.shape,VAR_DATA.min(),VAR_DATA.max()
#
				if iYEAR == 0:
					NTIMES         = DIMS[0]*NYEARS
					DIMS_MOD       = [NTIMES]
					for i in range(1,len(DIMS)):
						DIMS_MOD.append(DIMS[i])
					VAR_DATA_ALL   = np.zeros(DIMS_MOD)
#
				NCOUNT        = DIMS[0]
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:]     = VAR_DATA[:,:]
				elif len(DIMS) == 3:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:]   = VAR_DATA[:,:,:]
				elif len(DIMS) == 4:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:,:] = VAR_DATA[:,:,:,:]
#
			else:
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:]     = MISS_DATA
				elif len(DIMS) == 3:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:]   = MISS_DATA
				elif len(DIMS) == 4:
					VAR_DATA_ALL[COUNT:COUNT+NCOUNT,:,:,:] = MISS_DATA
#
# Increment as NCOUNT timepoints
#
			COUNT += NCOUNT
#
	VAR_DATA_ALL = np.squeeze(VAR_DATA_ALL)
#
# Return to calling routine
#
	return VAR_DATA_ALL
#
# ##################################################################################################
# (2.3) JULES offline - single multi-annual regional netCDF datasets
# ##################################################################################################
#
def get_JULES_data_regional_single(NREGIONS,NPOINTS,FILE_PREFIX,FILE_SUFFIX,JULES_VARINFO, \
	MISS_DATA,iFLAG,DEFAULT,DEBUG):
#
	COUNT       = 0
#
	for iREGION in range(NREGIONS):
#
		if 'm46' in FILE_PREFIX:
			FILE_CDF_IN   = FILE_PREFIX+('%04d' % (iREGION+1))+FILE_SUFFIX
		elif 'g03' in FILE_PREFIX:
			FILE_CDF_IN   = FILE_PREFIX+('%s' % iREGION)+FILE_SUFFIX
		print ('%3d %6d: ' % (iREGION,COUNT))+FILE_CDF_IN
#
		if iFLAG == 0:
			VAR_NAME,VAR_NAME_MOD,PLOT_LABEL,PLOT_MIN, \
				PLOT_MAX,PLOT_INC,CLEVELS,iSOIL,iFLAG= \
				get_JULES_variables(FILE_CDF_IN,JULES_VARINFO)
#		
# Skip if file missing
#
		if os.path.exists(FILE_CDF_IN):
			if os.path.getsize(FILE_CDF_IN) > 0:
#
				DIMS,VAR_DATA = \
					data_netCDF.data_netCDF_array_var(FILE_CDF_IN,VAR_NAME)
				if DEBUG == 'Y' and len(VAR_DATA) == 1:
					print VAR_DATA.shape,VAR_DATA.min(),VAR_DATA.max()
#
				if iREGION == 0:
					DIMS_MOD       = []
					for i in range(0,len(DIMS)-1):
						DIMS_MOD.append(DIMS[i])
					DIMS_MOD.append(NPOINTS)
					VAR_DATA_ALL   = np.zeros(DIMS_MOD)
#
				if len(DIMS) == 1:
					NCOUNT       = DEFAULT
				else:
					NCOUNT       = DIMS[-1]
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[:,COUNT:COUNT+NCOUNT]     = VAR_DATA[:,:]
				elif len(DIMS) == 3:
					VAR_DATA_ALL[:,:,COUNT:COUNT+NCOUNT]   = VAR_DATA[:,:,:]
				elif len(DIMS) == 4:
					VAR_DATA_ALL[:,:,:,COUNT:COUNT+NCOUNT] = VAR_DATA[:,:,:,:]
#
		if not os.path.exists(FILE_CDF_IN) or os.path.getsize(FILE_CDF_IN) <= 0:
#
			if len(DIMS) == 2:
				VAR_DATA_ALL[:,COUNT:COUNT+NCOUNT]     = MISS_DATA
			elif len(DIMS) == 3:
				VAR_DATA_ALL[:,:,COUNT:COUNT+NCOUNT]   = MISS_DATA
			elif len(DIMS) == 4:
				VAR_DATA_ALL[:,:,:,COUNT:COUNT+NCOUNT] = MISS_DATA
#
		print iREGION,COUNT,NCOUNT,NPOINTS
#
# Increment regions
#
		COUNT += NCOUNT
#
	VAR_DATA_ALL = np.squeeze(VAR_DATA_ALL)
	print VAR_DATA_ALL.shape
#
# Return to calling routine
#
	return VAR_DATA_ALL,VAR_NAME,VAR_NAME_MOD,PLOT_LABEL,PLOT_MIN,PLOT_MAX,PLOT_INC,CLEVELS,iSOIL
#
# ##################################################################################################
#
def get_JULES_data_regional_single2(NREGIONS,NPOINTS,FILE_PREFIX,FILE_SUFFIX,VAR_NAME, \
	MISS_DATA,DEFAULT,DEBUG):
#
	COUNT       = 0
#
	for iREGION in range(NREGIONS):
#
		if 'm46' in FILE_PREFIX or 'MAMM' in FILE_PREFIX:
			FILE_CDF_IN   = FILE_PREFIX+('%04d' % (iREGION+1))+FILE_SUFFIX
		elif 'g03' in FILE_PREFIX:
			FILE_CDF_IN   = FILE_PREFIX+('%s' % iREGION)+FILE_SUFFIX
		print ('%3d %6d: ' % (iREGION,COUNT))+FILE_CDF_IN
#
# Skip if file missing
#
		if os.path.exists(FILE_CDF_IN):
			if os.path.getsize(FILE_CDF_IN) > 0:
#
				DIMS,VAR_DATA = \
					data_netCDF.data_netCDF_array_var(FILE_CDF_IN,VAR_NAME)
				if DEBUG == 'Y' and len(VAR_DATA) == 1:
					print VAR_DATA.shape,VAR_DATA.min(),VAR_DATA.max()
#
				if iREGION == 0:
					DIMS_MOD       = []
					for i in range(0,len(DIMS)-1):
						DIMS_MOD.append(DIMS[i])
					DIMS_MOD.append(NPOINTS)
					VAR_DATA_ALL   = np.zeros(DIMS_MOD)
#
				if len(DIMS) == 1:
					NCOUNT       = DEFAULT
				else:
					NCOUNT       = DIMS[-1]
#
				if len(DIMS) == 2:
					VAR_DATA_ALL[:,COUNT:COUNT+NCOUNT]     = VAR_DATA[:,:]
				elif len(DIMS) == 3:
					VAR_DATA_ALL[:,:,COUNT:COUNT+NCOUNT]   = VAR_DATA[:,:,:]
				elif len(DIMS) == 4:
					VAR_DATA_ALL[:,:,:,COUNT:COUNT+NCOUNT] = VAR_DATA[:,:,:,:]
#
		if not os.path.exists(FILE_CDF_IN) or os.path.getsize(FILE_CDF_IN) <= 0:
#
			if len(DIMS) == 2:
				VAR_DATA_ALL[:,COUNT:COUNT+NCOUNT]     = MISS_DATA
			elif len(DIMS) == 3:
				VAR_DATA_ALL[:,:,COUNT:COUNT+NCOUNT]   = MISS_DATA
			elif len(DIMS) == 4:
				VAR_DATA_ALL[:,:,:,COUNT:COUNT+NCOUNT] = MISS_DATA
#
		print iREGION,COUNT,NCOUNT,NPOINTS
#
# Increment regions
#
		COUNT += NCOUNT
#
	VAR_DATA_ALL = np.squeeze(VAR_DATA_ALL)
	print VAR_DATA_ALL.shape
#
# Return to calling routine
#
	return VAR_DATA_ALL
#
# ##################################################################################################
# (3) JULES offline - regird to regular lat-lon grid
# ##################################################################################################
#
def JULES_data_regrid(VAR_DATA,NLON,NLAT,LON_DATA,LAT_DATA,LONG_START,LONG_END, \
	LAT_START,LAT_END,RESOL,iSOIL,MISS_DATA,DEBUG):
#
# Assign land points to regular lat-lon grid
#
	DIMS         = VAR_DATA.shape
	NPOINTS      = DIMS[-1]
	print VAR_DATA.shape,DIMS,NPOINTS
#
	if iSOIL == -1:
		VAR_DATA_REGRID = np.zeros((DIMS[0],NLAT,NLON))
		VAR_DATA_REGRID[:,:,:]   = MISS_DATA
	else:
		VAR_DATA_REGRID = np.zeros((DIMS[0],DIMS[1],NLAT,NLON))
		VAR_DATA_REGRID[:,:,:,:] = MISS_DATA
#
	for INDEX in range(NPOINTS):
#
		if LON_DATA[INDEX] > LONG_END:
			LON_DATA[INDEX] = LON_DATA[INDEX]-360.0
#
		iLON           = int((LON_DATA[INDEX]-LONG_START-RESOL/2.0)/RESOL)
		iLAT           = int((LAT_DATA[INDEX]-LAT_START -RESOL/2.0)/RESOL)
#
# Skip if out of range
#
		if iLON >= 0 and iLON < NLON and iLAT >= 0 and iLAT < NLAT:
#
			if iSOIL == -1:
				VAR_DATA_REGRID[:,iLAT,iLON]   = VAR_DATA[:,INDEX]
			else:
				VAR_DATA_REGRID[:,:,iLAT,iLON] = VAR_DATA[:,:,INDEX]
#
	return VAR_DATA_REGRID
#
# Return to calling routine
#
# ##################################################################################################
#
def JULES_data_regrid_mask(VAR_DATA,NLON,NLAT,LON_DATA,LAT_DATA,LONG_START,LONG_END, \
	LAT_START,LAT_END,RESOL,iSOIL,MISS_DATA,DEBUG):
#
# Assign land points to regular lat-lon grid
#
	DIMS           = VAR_DATA.shape
	NPOINTS        = DIMS[-1]
	print VAR_DATA.shape,DIMS,NPOINTS
#
	LAND_MASK      = np.empty((NLAT,NLON))
	LAND_MASK[:,:] = 1
#
	if iSOIL == -1:
		VAR_DATA_REGRID = np.zeros((DIMS[0],NLAT,NLON))
		VAR_DATA_REGRID[:,:,:]   = MISS_DATA
	else:
		VAR_DATA_REGRID = np.zeros((DIMS[0],DIMS[1],NLAT,NLON))
		VAR_DATA_REGRID[:,:,:,:] = MISS_DATA
#
	for INDEX in range(NPOINTS):
#
		if LON_DATA[INDEX] > LONG_END:
			LON_DATA[INDEX] = LON_DATA[INDEX]-360.0
#
		iLON           = int((LON_DATA[INDEX]-LONG_START-RESOL/2.0)/RESOL)
		iLAT           = int((LAT_DATA[INDEX]-LAT_START -RESOL/2.0)/RESOL)
#
		LAND_MASK[iLAT,iLON] = 0
#
# Skip if out of range
#
		if iLON >= 0 and iLON < NLON and iLAT >= 0 and iLAT < NLAT:
#
			if iSOIL == -1:
				VAR_DATA_REGRID[:,iLAT,iLON]   = VAR_DATA[:,INDEX]
			else:
				VAR_DATA_REGRID[:,:,iLAT,iLON] = VAR_DATA[:,:,INDEX]
#
	return VAR_DATA_REGRID,LAND_MASK
#
# Return to calling routine
