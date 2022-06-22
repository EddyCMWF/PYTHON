#!/bin/env python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os

# User input:
COMMON_tag='_ECOSSElayers'
EnW_tag=COMMON_tag+''
Scot_tag=COMMON_tag+''
out_tag=COMMON_tag+''

fill_value=-9999.

##################################################################################
# FUNCTIONS
####################################
# Define Merge Function
def MERGE_DATA(data1,data2,fill_value=fill_value):
    print('Merging Data')
    MERGED_DATA=np.zeros_like(data1.data)+fill_value
    # Do each layer independently
    for iSD in range(data1.shape[0]):
        data1_lay=data1[iSD,:]
        data2_lay=data2[iSD,:]
        data1_Mask=(data1_lay.mask==False)&(data2_lay.mask==True)
        data2_Mask=(data1_lay.mask==True)&(data2_lay.mask==False)
        BOTH_Mask=(data1_lay.mask==False)&(data2_lay.mask==False)
        
        # ratio of overlapping grid cells
        ratio=np.mean(data1_lay[BOTH_Mask])/np.mean(data2_lay[BOTH_Mask])
                
        MERGED_DATA[iSD,data1_Mask]=data1_lay[data1_Mask].data
        MERGED_DATA[iSD,data2_Mask]=data2_lay[data2_Mask].data*ratio

        MERGED_DATA[iSD,BOTH_Mask]= ( (data1_lay[BOTH_Mask].data ) + \
                                       data2_lay[BOTH_Mask].data*ratio )        \
                                     / 2.
    
    MERGED_DATA=np.ma.masked_equal(MERGED_DATA,fill_value)

    #plt.subplot(1,3,1)
    #plt.imshow(data1[0,:],origin='bottom')
    #plt.colorbar()
    #plt.subplot(1,3,2)
    #plt.imshow(data2[0,:],origin='bottom')
    #plt.colorbar()
    #plt.subplot(1,3,3)
    #plt.imshow(MERGED_DATA[0,:],origin='bottom')
    #plt.colorbar()
    #plt.show()

    return MERGED_DATA

# Define Fill missing values laterally function:
def FILL_GAPS_LATERALLY(data_dict,LAND_MASK):
    print('Filling Gaps Laterally')
    for var in data_dict:
        print(var)
        os.system('date')
        # loop over soil depths
        for iSD in range(data_dict[var].shape[0]):
            # Index of bad land points
            badex=np.where((data_dict[var].mask[iSD,:]==True)&(LAND_MASK==False))
            # Index of good land points
            goodex=np.where((data_dict[var].mask[iSD,:]==False)&(LAND_MASK==False))
            # convert good index into floats
            goodex = [ good.astype(float) for good in goodex ]
            # calculate the closest good point to each bad point
            #closest=[]
            #for iBAD in zip(badex[0],badex[1]):
            #    tempclosest=np.argmin(  np.abs(goodex[0]-iBAD[0]) \
            #                          + np.abs(goodex[1]-iBAD[1]) ) 
            #    closest.append(int(tempclosest))
            closest = [ int(np.argmin( ((goodex[0]-iBAD[0])**2) + ((goodex[1]-iBAD[1])**2) ) ) \
                            for iBAD in zip(badex[0],badex[1])  ]

            #convert good index back to integers
            goodex = [ good.astype(int) for good in goodex ]
            # replace bad values with closest good value and change mask
            data_dict[var].data[iSD,badex[0],badex[1]] =   \
                    data_dict[var].data[iSD,goodex[0][closest],goodex[1][closest]]
            data_dict[var].mask[iSD,badex[0],badex[1]] = False

        #plt.imshow(data_dict[var][0,:],origin='bottom')
        #plt.colorbar()
        #plt.show()
    return data_dict
#####################################################################################
# Directories and Options

SOIL_DIR='/prj/GREENHOUSE/SOIL_PROPERTIES/datasets/'

# Data to be merged 
Merged_COMPnames=['sand','silt','clay','org_carb','ph','Bulk_Density']

EnW_file=SOIL_DIR+'England_Wales_Soil_Data/LDE16_12_SRUC_Tarsitano/' + \
                    'EnW_Soil_WeightedComposition_CHESSgrid'+EnW_tag+'.nc'
Scot_file=SOIL_DIR+'Scotland_Soil_Data/' +   \
                    'Scot_Soil_WeightedComposition_CHESSgrid'+Scot_tag+'.nc'

COMP_outfile=SOIL_DIR+'Merged_Soil_WeightedComposition_CHESSgrid'+out_tag+'.nc'

CHESS_landcover_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'

####################################################################################
#Read in latlon/xy data from chess_landcover
LLinf=nc.Dataset(CHESS_landcover_file,'r')
landcover=LLinf.variables['frac'][:]
LLinf.close()
LAND_MASK=landcover.mask[0,:]

#Read in all EnW data
print('Reading: '+EnW_file)
EnW_datadict={}
EnWinf=nc.Dataset(EnW_file,'r')
for var in Merged_COMPnames:
    EnW_datadict[var]=EnWinf.variables[var][:]
EnWinf.close()

#Read in all Scot data
print('Reading: '+Scot_file)
Scot_datadict={}
Scotinf=nc.Dataset(Scot_file,'r')
for var in Merged_COMPnames:
    try:
        Scot_datadict[var]=Scotinf.variables[var][:]
    except:
        pass
Scotinf.close()
Scot_datadict['Bulk_Density']=1.772 - (0.4127*np.log(Scot_datadict['org_carb']))

#print(type(Scot_datadict['Bulk_Density']))
#for iSD in range(4):
#    plt.subplot(2,2,iSD+1)
#    plt.imshow(Scot_datadict['Bulk_Density'][iSD,:],origin='bottom')
#    plt.colorbar()
#plt.show()

####################################################################################
# Merge datasets, using function defined at top of script
COMP_datadict={}
for param in Merged_COMPnames:    
    COMP_datadict[param]=MERGE_DATA( EnW_datadict[param], \
                                     Scot_datadict[param] ).copy()

# Fill gaps laterally using function at top of script
COMP_datadict=FILL_GAPS_LATERALLY(COMP_datadict,LAND_MASK)

# Renormalise  sand+silt+clay to 100% 
SSC_total=COMP_datadict['sand']+COMP_datadict['silt']+COMP_datadict['clay']
for var in ['sand','silt','clay']:
    COMP_datadict[var] = COMP_datadict[var]/(SSC_total*0.01)

# Open template file to copy over attribute and dimensions
templatef=nc.Dataset(EnW_file,'r')

# Output Composition data
print("Writing composition data to: "+COMP_outfile)
outf=nc.Dataset(COMP_outfile,'w')
# create dimensions
for dim in templatef.dimensions:
    outf.createDimension(str(dim),len(templatef.dimensions[dim]))
    outvar=outf.createVariable(str(dim),'float32',str(dim))
    for att in templatef.variables[dim].ncattrs():
        outvar.setncattr(str(att),templatef.variables[dim].getncattr(str(att)))
    outvar[:]=templatef.variables[dim][:]
    
    
for var in COMP_datadict:
    invar=templatef.variables[var]
    outvar=outf.createVariable(var,'float32',invar.dimensions)
    for att in invar.ncattrs():
        outvar.setncattr(str(att),invar.getncattr(str(att)))
    outvar[:]=COMP_datadict[var]
    
outf.title='Merged Soil Composition map for the UK'
outf.data1='Cranfield England and Wales Soil Map'
outf.data2='Scotland Soils'
outf.method='Soil Composition aggregated on to the CHESS grid'
outf.MergedBy='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.organisation='Centre of Ecology and Hydrology'

outf.close()


templatef.close()




