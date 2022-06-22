# Python module to convert points to a different projection

from osgeo import ogr
from osgeo import osr
import numpy as np
import shapefile
from matplotlib.patches import Polygon

def define_projections(EPSG_IN=4326,EPSG_OUT=27700):

# Define the source projection, e.g. WGS lat/lon
 
	PROJ_IN     = osr.SpatialReference() # Define a SpatialReference object
	PROJ_IN.ImportFromEPSG( EPSG_IN )

# Define the target projection, e.g., OS national grid

	PROJ_OUT    = osr.SpatialReference() # Define the SpatialReference object
	PROJ_OUT.ImportFromEPSG( EPSG_OUT )

# Define a coordinate transformtion object, *from* EPSG_IN *to* EPSG_OUT

	TX          = osr.CoordinateTransformation(PROJ_IN,PROJ_OUT)

# Return to calling routine

	return PROJ_IN,PROJ_OUT,TX

def convert_points(PROJ_IN,PROJ_OUT,TX,X_IN,Y_IN,DEBUG):

	NX          = X_IN.shape[0]
	NY          = X_IN.shape[1]
	X_OUT       = np.zeros((NX,NY))
	Y_OUT       = np.zeros((NX,NY))

# Loop over X_IN and Y_IN

	for IX in range(NX):
		for IY in range(NY):

# Do the transformation using the TransformPoint method

			XPOINT_IN,YPOINT_IN          = float(X_IN[IX,IY]), float(Y_IN[IX,IY])
    			X_OUT[IX,IY], Y_OUT[IX,IY],Z = TX.TransformPoint ( XPOINT_IN,YPOINT_IN )

# Print out
			if DEBUG == 'Y':
				TEXT    = '%6d %6d %12.6f %12.6f %12.6f %12.6f' % \
					(IX+1,IY+1,X_IN[IX,IY],Y_IN[IX,IY],X_OUT[IX,IY],Y_OUT[IX,IY])
    				print TEXT

	return X_OUT,Y_OUT
#
def get_shapefile(FILE_SHAPE,SCALE,DEBUG):
#
# Routine to read and parse shapefile
#
	SHAPE_FID   = shapefile.Reader(FILE_SHAPE)
#
	SHAPES      = SHAPE_FID.shapes()
	NSHAPES     = len(SHAPES)
#
	PATCHES_ALL = []
	for iSHAPE in xrange(NSHAPES):
		PATCHES     = []
		POINTS      = np.array(SHAPES[iSHAPE].points)
		PART        = SHAPES[iSHAPE].parts
		PAR         = list(PART)+[POINTS.shape[0]]
#
		for INDEX in xrange(len(PART)):
			PATCHES.append(Polygon(POINTS[PAR[INDEX]:PAR[INDEX+1]]/SCALE))
			if DEBUG == 'Y':
				print NSHAPES,INDEX,PART,PAR
				print POINTS[PAR[INDEX]:PAR[INDEX+1]]
#
		PATCHES_ALL.append(PATCHES)
#
	return PATCHES_ALL
