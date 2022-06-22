#!/usr/bin/python

import gdal
import numpy as np
import netCDF4 as nc
from scipy import stats

def ReadBilFile(bil):
    gdal.GetDriverByName('EHdr').Register()
    img = gdal.Open(bil)
    band = img.GetRasterBand(1)
    data = band.ReadAsArray()
    #geotransform = img.GetGeoTransform()
    #xstart = geotransform[0]
    #xstep  = geotransform[1]
    #xsize  = img.RasterXSize
    #ystart = geotransform[3]
    #ystep  = geotransform[5]
    #ysize  = img.RasterYSize
    #lons=(np.arange(xsize)*xstep)+xstart
    #lats=(np.arange(ysize)*ystep)+ystart
    #return data, lats, lons
    return data


HWSD_dir   = '/users/eow/edwcom/data/HWSD/'
HWSD_file  = 'hwsd.bil'
HWSD_param_file = 'HWSD_DATA.txt'
HWSD_SOC_name = '"S_OC"'
HWSD_MU_name  = '"MU_GLOBAL"'

CRU_DIR    = '/users/eow/edwcom/CRUNCEP/'
CRU_SC_outfile = 'cruncep-HWSD-soilcarbon.nc'
CRUNCEP_gfile = 'cru_ncep_land.nc'

#Read in HWSD data using function above
HWSD_data=ReadBilFile(HWSD_dir+HWSD_file)
#HWSD_data,HWSD_lats,HWSD_lons=ReadBilFile(HWSD_dir+HWSD_file)

HWSD_param_inf = open(HWSD_dir+HWSD_param_file)
HWSD_param_data = HWSD_param_inf.readlines()
HWSD_param_inf.close()

HWSD_param_hdrs = HWSD_param_data.pop(0)
HWSD_param_hdrs = HWSD_param_hdrs.split(',')
SOC_index      = HWSD_param_hdrs.index(HWSD_SOC_name)
MU_index       = HWSD_param_hdrs.index(HWSD_MU_name)

HWSD_SOC_param = []
HWSD_MU_index  = []
for line in HWSD_param_data:
    split=line.split(',')
    HWSD_MU_index.append(float(split[MU_index]))
    HWSD_SOC_param.append(float(split[SOC_index]))


# Create 0.5 degree HWSD grid and SC map based on modal HWSD class
# if HWSD class has invalid SC then choose second modal value
HWSD_0p5deg        = np.zeros([360,720])  # Array of the most common HWSD calss at 0.5 degrees
pseudo_HWSD_0p5deg = np.zeros([360,720])-999.0  # Array of the HWSD class used when allocating SC
                         # This will not be the above when the most common class does not have a SC
HWSD_SC_0p5deg     = np.zeros([360,720])-999.0 # Array to store the SC
HWSD_SCmean_0p5deg     = np.zeros([360,720])-999.0 # Array to store the SC

# loop around each point
for i in range(720):
    for j in range(360):
        # cut out relevent section
        section=HWSD_data[j*60.:(j*60.)+60.,i*60.:(i*60.)+60.].flat
        
        # find the modal value and store in array
        HWSD_0p5deg[j,i]=stats.mode(section)[0]
        
        
        # if not water then get soil carbon value:
        if (HWSD_0p5deg[j,i]!=0):
            index = HWSD_MU_index.index(HWSD_0p5deg[j,i])
            HWSD_SC_0p5deg[j,i]=HWSD_SOC_param[index]
            
            # if soil carbon is invalid we take second most common soil class 
            # (if there is one with more than 10% cover)
            if ((HWSD_SC_0p5deg[j,i]==-9999.0) & \
                (len(section[np.where(section!=HWSD_0p5deg[j,i])])>360)):

                pseudo_HWSD_0p5deg[j,i]=stats.mode(section[np.where(section!=HWSD_0p5deg[j,i])])[0]
                # check it's not water
                if (pseudo_HWSD_0p5deg[j,i]!=0):
                    index = HWSD_MU_index.index(pseudo_HWSD_0p5deg[j,i])
                    HWSD_SC_0p5deg[j,i]=HWSD_SOC_param[index]
                else:
                    #if it is water re-set to -999.0
                    pseudo_HWSD_0p5deg[j,i]=-999.0
                    HWSD_SC_0p5deg[j,i]=-999.0
                
                # if soil carbon is still invalid we take third soil class
                if (HWSD_SC_0p5deg[j,i]==-9999.0):
                    sub_section_index=np.where((section!=HWSD_0p5deg[j,i]) & \
                                               (section!=pseudo_HWSD_0p5deg[j,i]))
                    # only do if greater than 10% of box, i.e. 60*60*0.1
                    if (len(sub_section_index)>360):
                        pseudo_HWSD_0p5deg[j,i]=stats.mode(section[sub_section_index])
                        index = HWSD_MU_index.index(pseudo_HWSD_0p5deg[j,i])
                        HWSD_SC_0p5deg[j,i]=HWSD_SOC_param[index]
                    # else we set to 0.0
                    else:
                        pseudo_HWSD_0p5deg[j,i]=HWSD_0p5deg[j,i]
                        HWSD_SC_0p5deg[j,i]=0.0
            
            # if the second soil class is less than 10%, we set soil carbon to zero
            elif (HWSD_SC_0p5deg[j,i]==-9999.0):
                HWSD_SC_0p5deg[j,i]=0.0
                pseudo_HWSD_0p5deg[j,i]=HWSD_0p5deg[j,i]
            else:
                pseudo_HWSD_0p5deg[j,i]=HWSD_0p5deg[j,i]
            
        
HWSD_0p5deg[np.where(HWSD_0p5deg==0)]=-999.0
HWSD_0p5deg=np.ma.masked_equal(HWSD_0p5deg,-999.0)
pseudo_HWSD_0p5deg=np.ma.masked_equal(pseudo_HWSD_0p5deg,-999.0)
HWSD_SC_0p5deg=np.ma.masked_equal(HWSD_SC_0p5deg,-999.0)





