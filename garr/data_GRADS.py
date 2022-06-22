# Python module to selected datasets
#
# Garry Hayman
# Centre for Ecology and Hydrology
# December 2011
#
# Contains
#
def setup_GRADS(iSYSTEM):
#
	import os
	import sys
#
	GRADS_DIR      = '/data/grp/eow/garr/Projects/Methane/CODE/GRADS'
	GADDIR         = '/users/global/dbcl/packages/grads/grads'
	GASCRP         = '/users/global/dbcl/packages/grads/lib'
#
	os.chdir(GRADS_DIR)
#
	if iSYSTEM=='UNIX':
		os.putenv('GADDIR',GADDIR)
		os.putenv('GASCRP',GASCRP)
		PATH           = os.getenv('PATH')
		PATH           = PATH+':'+GADDIR+':'+GASCRP
		os.putenv('PATH',PATH)
#
# Return to calling routine
#
	return
#
def data_GRADS(FILE_CTL,DEBUG):
#
	import grads
	import numpy
	from grads.ganum import GaNum
#
	ga            = grads.GrADS(Bin='grads',Window=False)
#	ga            = GaNum(Bin='gradsnc')
#
	try:
#
		FILE_HAND     = ga.open(FILE_CTL)
#
		print FILE_HAND.title
		print 'File Id              : ', FILE_HAND.fid
		print 'Dataset type         : ', FILE_HAND.type
		print 'No. of Time steps    : ', FILE_HAND.nt
		print 'No. of Longitudes    : ', FILE_HAND.nx
		print 'No. of Latitudes     : ', FILE_HAND.ny
		print 'No. of Levels        : ', FILE_HAND.nz
		print 'Variable  names      : ', FILE_HAND.vars
		print 'Variable levels      : ', FILE_HAND.var_levs
		print 'Variable titles      : ', FILE_HAND.var_titles
		print 
		print ">>> OK <<< open CTL file"
	except:
		print ">>> NOT OK <<< cannot open CTL file"
#
	VAR_NAMES     = []
	print('Input variable names: ')
	INPUT         = input()
	for i in range(len(INPUT)): VAR_NAMES.append(INPUT[i])
#
	ga('set t 1 last')
#
	NUM_VARS      = len(VAR_NAMES)
	VAR_DATA_ALL  = []
#
	for iVAR in range (NUM_VARS):
#
		VAR_NAME      = VAR_NAMES[iVAR]
		print(" ")
		print('%4d %s' % (iVAR,VAR_NAME))
#
		VAR_DATA      = ga.exp(VAR_NAME)
		VAR_DATA_ALL.append(VAR_DATA)
#
		if iVAR == 0:
			LONG          = VAR_DATA.grid.lon
			LAT           = VAR_DATA.grid.lat
			STIME         = VAR_DATA.grid.time
#
		if DEBUG=='Y':
			print(VAR_DATA.shape)
			print(len(LONG),LONG[0],LONG[len(LONG)-1])
			print(len(LAT),LAT[0],LAT[len(LAT)-1])
			print(len(STIME),STIME[0],STIME[len(STIME)-1])
#
# Return to calling routine
#
	return LONG,LAT,STIME,VAR_NAMES,VAR_DATA_ALL
#
def data_GRADS_var(FILE_CTL,VAR_NAME,DEBUG):
#
	import grads
	import numpy
	from grads.ganum import GaNum
#
	ga            = grads.GrADS(Bin='grads',Window=False)
#	ga            = GaNum(Bin='gradsnc')
#
	try:
#
		FILE_HAND     = ga.open(FILE_CTL)
#
		print FILE_HAND.title
		print 'File Id              : ', FILE_HAND.fid
		print 'Dataset type         : ', FILE_HAND.type
		print 'No. of Time steps    : ', FILE_HAND.nt
		print 'No. of Longitudes    : ', FILE_HAND.nx
		print 'No. of Latitudes     : ', FILE_HAND.ny
		print 'No. of Levels        : ', FILE_HAND.nz
		print 'Variable  names      : ', FILE_HAND.vars
		print 'Variable levels      : ', FILE_HAND.var_levs
		print 'Variable titles      : ', FILE_HAND.var_titles
		print 
		print ">>> OK <<< open CTL file"
	except:
		print ">>> NOT OK <<< cannot open CTL file"
#
	ga('set t 1 last')
#
	iVAR          = 0 
	print(" ")
	print('%4d %s' % (iVAR,VAR_NAME))
#
	VAR_DATA      = ga.exp(VAR_NAME)
#
	LONG          = VAR_DATA.grid.lon
	LAT           = VAR_DATA.grid.lat
	STIME         = VAR_DATA.grid.time
#
	if DEBUG=='Y':
		print(VAR_DATA.shape)
		print(len(LONG),LONG[0],LONG[len(LONG)-1])
		print(len(LAT),LAT[0],LAT[len(LAT)-1])
		print(len(STIME),STIME[0],STIME[len(STIME)-1])
#
# Return to calling routine
#
	return LONG,LAT,STIME,VAR_DATA
