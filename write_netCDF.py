#
# Python module to write to a netCDF file
#
# Garry Hayman
# Centre for Ecology and Hydrology
# February 2011
#
import sys
from netCDF4 import Dataset
from numpy   import dtype
#
# Contains
#
def write_netCDF(FILE_CDF,NDIMS,DIM_INFO,DIM_DATA,VAR_INFO,VAR_DEP,VAR_DATA,MISS_DATA):
#
# Write out netCDF file
#
	print FILE_CDF
	NC_FID     = Dataset(FILE_CDF,'w',format='NETCDF3_CLASSIC')
	setattr(NC_FID,'history','netCDF file from Python')
#
# Define the co-ordinate variables
#
	for iDIM in range(NDIMS):
		NC_FID.createDimension(DIM_INFO[iDIM][1],DIM_INFO[iDIM][0])
		DIM_CDF       = NC_FID.createVariable(DIM_INFO[iDIM][1],dtype(DIM_INFO[iDIM][2]).char,(DIM_INFO[iDIM][1]))
		DIM_CDF.units = DIM_INFO[iDIM][3]
		DIM_CDF[:]    = DIM_DATA[iDIM]
#
# Create and write data to variable.
#
	print VAR_DEP
	DATA_CDF    = NC_FID.createVariable(VAR_INFO[0],dtype(VAR_INFO[1]).char,(VAR_DEP))
#
	print DIM_INFO[0][0]
#	for iDIM in range(DIM_INFO[0][0]):
#		DATA_CDF[iDIM,:,:] = VAR_DATA[iDIM,:,:]
	DATA_CDF[:] = VAR_DATA
#
	setattr(DATA_CDF,'missing_value',MISS_DATA)
#
# Close the file.
#
	NC_FID.close()
#
	TEXT        = '*** SUCCESS writing file '+FILE_CDF
	print TEXT
#
	return
#
def write_netCDF_multi(FILE_CDF,NDIMS,DIM_INFO,DIM_DATA,NVAR,VAR_INFO_ALL,VAR_DEP,VAR_DATA_ALL,MISS_DATA):
#
# Write out netCDF file
#
	NC_FID     = Dataset(FILE_CDF,'w',format='NETCDF3_CLASSIC')
	setattr(NC_FID,'history','netCDF file from Python')
#
# Define the co-ordinate variables
#
	for iDIM in range(NDIMS):
		print DIM_INFO[iDIM]
		NC_FID.createDimension(DIM_INFO[iDIM][1],DIM_INFO[iDIM][0])
		DIM_CDF       = NC_FID.createVariable(DIM_INFO[iDIM][1],dtype(DIM_INFO[iDIM][2]).char,(DIM_INFO[iDIM][1]))
		DIM_CDF.units = DIM_INFO[iDIM][3]
		DIM_CDF[:]    = DIM_DATA[iDIM]
#
# Create and write data to variable.
#
	for iVAR in range(NVAR):
#
#		VAR_DEP       = VAR_DEP_ALL[iVAR]
		print VAR_DEP
#
		VAR_INFO      = VAR_INFO_ALL[iVAR]
		print VAR_INFO
		DATA_CDF      = NC_FID.createVariable(VAR_INFO[0],dtype(VAR_INFO[1]).char,(VAR_DEP))
#
		print DIM_INFO[0][0]
		DATA_CDF[:] = VAR_DATA_ALL[iVAR]
#
		setattr(DATA_CDF,'missing_value',MISS_DATA)
		if len(VAR_INFO[2]) != 0:
			setattr(DATA_CDF,'units',VAR_INFO[2])
#
# Close the file.
#
	NC_FID.close()
#
	TEXT        = '*** SUCCESS writing file '+FILE_CDF
	print TEXT
#
	return
#
def write_netCDF_multi_new(FILE_CDF,NDIMS,DIM_INFO,DIM_DATA,NVAR,VAR_INFO_ALL,VAR_DEP_ALL,VAR_DATA_ALL,MISS_DATA_ALL):
#
# Write out netCDF file
#
	NC_FID     = Dataset(FILE_CDF,'w',format='NETCDF3_CLASSIC')
	setattr(NC_FID,'history','netCDF file from Python')
#
# Define the co-ordinate variables
#
	for iDIM in range(NDIMS):
		print DIM_INFO[iDIM]
		NC_FID.createDimension(DIM_INFO[iDIM][1],DIM_INFO[iDIM][0])
		DIM_CDF       = NC_FID.createVariable(DIM_INFO[iDIM][1],dtype(DIM_INFO[iDIM][2]).char,(DIM_INFO[iDIM][1]))
		DIM_CDF.units = DIM_INFO[iDIM][3]
		DIM_CDF[:]    = DIM_DATA[iDIM]
#
# Create and write data to variable.
#
	for iVAR in range(NVAR):
#
		print iVAR
		VAR_DEP       = VAR_DEP_ALL[iVAR]
		print VAR_DEP
#
		VAR_INFO      = VAR_INFO_ALL[iVAR]
		print VAR_INFO
		DATA_CDF      = NC_FID.createVariable(VAR_INFO[0],dtype(VAR_INFO[1]).char,(VAR_DEP))
#
		print DIM_INFO[0][0]
		DATA_CDF[:] = VAR_DATA_ALL[iVAR]
#
		setattr(DATA_CDF,'missing_value',MISS_DATA_ALL[iVAR])
		if len(VAR_INFO[2]) != 0:
			setattr(DATA_CDF,'units',VAR_INFO[2])
#
# Close the file.
#
	NC_FID.close()
#
	TEXT        = '*** SUCCESS writing file '+FILE_CDF
	print TEXT
#
	return
#
def write_netCDF_variable_append(FILE_CDF,VAR_INFO,VAR_DEP,VAR_DATA,MISS_DATA):
#
# Write out netCDF file
#
	NC_FID     = Dataset(FILE_CDF,'a',format='NETCDF3_CLASSIC')
#
# Create and write data to variable.
#
	print VAR_DEP
	DATA_CDF    = NC_FID.createVariable(VAR_INFO[0],dtype(VAR_INFO[1]).char,(VAR_DEP))
	DATA_CDF[:] = VAR_DATA
#
	setattr(DATA_CDF,'missing_value',MISS_DATA)
	setattr(DATA_CDF,'units',VAR_INFO[2])
#
# Close the file.
#
	NC_FID.close()
#
	TEXT        = '*** SUCCESS writing file '+FILE_CDF
	print TEXT
#
	return
#
