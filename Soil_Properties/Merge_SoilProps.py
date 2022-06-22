#!/bin/env python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fill_value=-9999.

##################################################################################
# FUNCTIONS
####################################
# Define Merge Function
def MERGE_DATA(data1,data2,fill_value=fill_value):
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

    return np.ma.masked_equal(MERGED_DATA,fill_value)

# Define Fill missing values laterally function:
def FILL_GAPS_LATERALLY(data_dict,LAND_MASK):
    for var in data_dict:
        # loop over soil depths
        for iSD in range(data_dict[var].shape[0]):
            # Index of bad land points
            badex=np.where((data_dict[var].mask[iSD,:]==True)&(LAND_MASK==False))
            # Index of good land points
            goodex=np.where((data_dict[var].mask[iSD,:]==False)&(LAND_MASK==False))
            # convert good index into floats
            goodex = [ good.astype(float) for good in goodex ]
            # calculate the closest good point to each bad point
            closest=[]
            for iBAD in zip(badex[0],badex[1]):
                tempclosest=np.argmin(  np.abs(goodex[0]-iBAD[0]) \
                                      + np.abs(goodex[1]-iBAD[1]) ) 
                closest.append(int(tempclosest))
            #convert good index back to integers
            goodex = [ good.astype(int) for good in goodex ]
            # replace bad values with closest good value and change mask
            data_dict[var].data[iSD,badex[0],badex[1]] =   \
                    data_dict[var].data[iSD,goodex[0][closest],goodex[1][closest]]
            data_dict[var].mask[iSD,badex[0],badex[1]] = False
    return data_dict
#####################################################################################
# Directories and Options
EnW_tag=''
Scot_tag=''
out_tag=''

SOIL_DIR='/prj/GREENHOUSE/SOIL_PROPERTIES/datasets/'

# Data to be merged 
Merged_COMPnames=['sand','silt','clay','org_carb','ph','Bulk_Density']
Merged_BCnames=['BC_sm_sat','BC_hcon','BC_hcap','BC_sm_wilt', \
                'BC_bexp','BC_sathh','BC_satcon','BC_sm_crit']
Merged_VGnames=['VG_oneovernminusone','VG_ksat','VG_hcap','VG_sm_wilt', \
                'VG_oneoveralpha','VG_sm_crit','VG_hcon','VG_sm_sat'   ]

EnW_file=SOIL_DIR+'England_Wales_Soil_Data/LDE16_12_SRUC_Tarsitano/' + \
                    'EnW_Soil_WeightedCompositionProperties_CHESSgrid'+EnW_tag+'.nc'
Scot_file=SOIL_DIR+'Scotland_Soil_Data/' +   \
                    'Scot_Soil_WeightedCompositionProperties_CHESSgrid'+Scot_tag+'.nc'

COMP_outfile=SOIL_DIR+'Merged_Soil_WeightedComposition_CHESSgrid'+out_tag+'.nc'
BC_outfile=SOIL_DIR+'Merged_Soil_WeightedBCproperties_CHESSgrid'+out_tag+'.nc'
VG_outfile=SOIL_DIR+'Merged_Soil_WeightedVGproperties_CHESSgrid'+out_tag+'.nc'

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
for var in EnWinf.variables:
    EnW_datadict[var]=EnWinf.variables[var][:]
EnWinf.close()

#Read in all Scot data
print('Reading: '+Scot_file)
Scot_datadict={}
Scotinf=nc.Dataset(Scot_file,'r')
for var in Scotinf.variables:
    Scot_datadict[var]=Scotinf.variables[var][:]
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

BC_datadict={}
for param in Merged_BCnames:    
    BC_datadict[param]=MERGE_DATA( EnW_datadict[param], \
                                   Scot_datadict[param] ).copy()

VG_datadict={}
for param in Merged_VGnames:
    VG_datadict[param]=MERGE_DATA( EnW_datadict[param], \
                                   Scot_datadict[param] ).copy()
    

# Fill gaps laterally using function at top of script
COMP_datadict=FILL_GAPS_LATERALLY(COMP_datadict,LAND_MASK)

VG_datadict=FILL_GAPS_LATERALLY(VG_datadict,LAND_MASK)

BC_datadict=FILL_GAPS_LATERALLY(BC_datadict,LAND_MASK)

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


# Output BC Properties data 
print("Writing Brooks and Corey data to: "+BC_outfile)
outf=nc.Dataset(BC_outfile,'w')
# create dimensions
for dim in templatef.dimensions:
    outf.createDimension(str(dim),len(templatef.dimensions[dim]))
    outvar=outf.createVariable(str(dim),'float32',str(dim))
    for att in templatef.variables[dim].ncattrs():
        outvar.setncattr(str(att),templatef.variables[dim].getncattr(str(att)))
    outvar[:]=templatef.variables[dim][:]
    
    
for var in BC_datadict:
    invar=templatef.variables[var]
    outvar=outf.createVariable(var,'float32',invar.dimensions)
    outvar[:]=BC_datadict[var]
    
outf.title='Merged Brooks and Corey Soil Properties map for the UK'
outf.data1='Cranfield England and Wales Soil Map'
outf.data2='Scotland Soils'
outf.method='Soil Properties calcualted at the soil type level, then aggregated on to the CHESS grid'
outf.MergedBy='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.organisation='Centre of Ecology and Hydrology'

outf.close()



# Output VG Properties data 
print("Writing Van Genuchten data to: "+VG_outfile)
outf=nc.Dataset(VG_outfile,'w')
# create dimensions
for dim in templatef.dimensions:
    outf.createDimension(str(dim),len(templatef.dimensions[dim]))
    outvar=outf.createVariable(str(dim),'float32',str(dim))
    for att in templatef.variables[dim].ncattrs():
        outvar.setncattr(str(att),templatef.variables[dim].getncattr(str(att)))
    outvar[:]=templatef.variables[dim][:]
    
    
for var in VG_datadict:
    invar=templatef.variables[var]
    outvar=outf.createVariable(var,'float32',invar.dimensions)
    outvar[:]=VG_datadict[var]
    
outf.title='Merged Van Genuchten Soil Properties map for the UK'
outf.data1='Cranfield England and Wales Soil Map'
outf.data2='Scotland Soils'
outf.method='Soil Properties calcualted at the soil type level, then aggregated on to the CHESS grid'
outf.MergedBy='Edward Comyn-Platt, edwcom@ceh.ac.uk'
outf.organisation='Centre of Ecology and Hydrology'

outf.close()


templatef.close()




