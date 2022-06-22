# Python module to derive scale factors for ACCMIP anthropogenic emissions
#
import numpy as np
#
def read_file(FILE,NTRANS,START_YEAR,END_YEAR,DEBUG):
#
	iYEAR      = -1
        iTRAN      = -1
	iREAD      = 'N'
#
	NYEARS     = END_YEAR-START_YEAR+1
	EMISSIONS  = np.zeros((NTRANS+2,NYEARS))
#
	for LINE in open(FILE):
#
#		if DEBUG == 'Y': print LINE.replace('\n','')
#
		if iREAD == 'Y' and (iTRAN >= 0 and iTRAN < NTRANS+2):
			EMISSIONS[iTRAN,iYEAR] = float(LINE.replace('\n',''))
			if DEBUG == 'Y': print iYEAR,iTRAN,EMISSIONS[iTRAN,iYEAR]
			iTRAN    +=  1
#
		if iTRAN == NTRANS+2:
			iTRAN     =  0
			iREAD     = 'N' 
#
		if 'Domain' in LINE:
			YEAR       = int(LINE[28:32])
			if DEBUG == 'Y': print YEAR
			if YEAR >= START_YEAR and YEAR <= END_YEAR:
				iREAD     = 'Y'
				iTRAN     =  0
				iYEAR    +=  1
#
	return EMISSIONS			
# 
def derive_factors(FILES,NTRANS,START_YEAR,END_YEAR,DEBUG):
#
	NFILES     = len(FILES)
	NYEARS     = END_YEAR-START_YEAR+1
	FACTORS    = np.zeros((NTRANS))
	EMISSIONS  = np.zeros((NFILES,NTRANS+2,NYEARS))
#
	for iFILE in range(NFILES):
		EMISSIONS[iFILE,:,:] = read_file(FILES[iFILE],NTRANS,START_YEAR,END_YEAR,DEBUG)
#
	SCALE_TOT0 = EMISSIONS[0,0,:].sum()/EMISSIONS[0,-1,:].sum()
	SCALE_TOT1 = EMISSIONS[1,0,:].sum()/EMISSIONS[1,-1,:].sum()
#	SCALE_TOT0 = EMISSIONS[0,0,:].sum()
#	SCALE_TOT1 = EMISSIONS[1,0,:].sum()
#
	if DEBUG == 'Y':
		print EMISSIONS
		print SCALE_TOT0,SCALE_TOT1
#
# Assign factors for each TRANSCOM region
#
	for iTRAN in range(NTRANS):
		FACTORS[iTRAN] = EMISSIONS[0,iTRAN+1,:].sum()*SCALE_TOT0/ \
				(EMISSIONS[1,iTRAN+1,:].sum()*SCALE_TOT1) 
#
	return FACTORS
#
def read_scale_factors(FILE,NTRANS,iSPECIES,START_YEAR,END_YEAR,DEBUG):
#
        iTRAN      = -1
	iREAD      = 'N'
#
	NYEARS     = END_YEAR-START_YEAR+1
	FACTORS    = np.zeros((NTRANS))
#
	for LINE in open(FILE):
#
#		if DEBUG == 'Y': print LINE.replace('\n','')
#
		if iREAD == 'Y':
			INPUT_LINE             = LINE.replace('\n','').split()
			FACTORS[iTRAN]         = float(INPUT_LINE[3+iSPECIES])
			if DEBUG == 'Y': print iTRAN+1,FACTORS[iTRAN]
			iREAD     = 'N'
#
		if 'Region' in LINE:
				iREAD     = 'Y'
				iTRAN    +=  1
#
	return FACTORS
# 
