#!/bin/env python
import netCDF4 as nc
import numpy as np
from SoilTools import BrooksCorey as BC
from SoilTools import VanGenuchten as VG
from SoilTools import Thermal as TH

#######################################################################################
# Directories and Filenames
tag='_CultAR'
Input_Dir='/users/eow/edwcom/GREENHOUSE/SOIL_PROPERTIES/datasets/'
Comp_File=Input_Dir+'Merged_Soil_WeightedComposition_CHESSgrid'+tag+'.nc'

Output_Dir=Input_Dir
OUT_BC_file=Output_Dir+'Merged_Soil_WeightedComposition_BCproperties_CHESSgrid'+tag+'.nc'
OUT_VG_file=Output_Dir+'Merged_Soil_WeightedComposition_VGproperties_CHESSgrid'+tag+'.nc'

CHESS_landcover_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'

#######################################################################################
# Options
fill_value=-9999.
Soil_Layer_Thick=np.array([0.1,0.25,0.65,2.0])   #(metres)
Soil_Layer_Depth=np.array([0.1,0.23,1.0,3.0]) #(metres)
nSD=len(Soil_Layer_Thick)

#######################################################################################
# Read in latlon/xy data from chess_landcover
LLinf=nc.Dataset(CHESS_landcover_file,'r')
landcover=LLinf.variables['frac'][:]
LLinf.close()
LAND_MASK=landcover.mask[0,:]


#######################################################################################
# Read Soil Composition File
INPUT_vars=['sand','silt','clay','org_carb']
COMP_data={}
COMPinf=nc.Dataset(Comp_File,'r')
for var in INPUT_vars:
    COMP_data[var]=COMPinf.variables[var][:]


#######################################################################################
# Compute the Brooks and Corey Soil Properties
BC_Properties=BC.get_BC_soil_properties(COMP_data['clay'], \
                                        COMP_data['sand'], \
                                        COMP_data['silt']  )

BC_Properties['hcap']=TH.hcap(COMP_data['clay'], \
                              COMP_data['sand'], \
                              COMP_data['silt'], \
                              sm_sat=BC_Properties['sm_sat'] )

BC_Properties['hcon']=TH.hcon_Farouki(COMP_data['clay'], \
                                      COMP_data['sand'], \
                                      COMP_data['silt'], \
                                      sm_sat=BC_Properties['sm_sat'] )

print('Outputting BC data to: '+OUT_BC_file)
# Output BC data
outf=nc.Dataset(OUT_BC_file,'w')

# create dimensions
for dim in COMPinf.dimensions:
    outf.createDimension(str(dim),len(COMPinf.dimensions[dim]))
    outvar=outf.createVariable(str(dim),'float32',str(dim))
    for att in COMPinf.variables[dim].ncattrs():
        outvar.setncattr(str(att),COMPinf.variables[dim].getncattr(str(att)))
    outvar[:]=COMPinf.variables[dim][:]

for var in BC_Properties:
    outvar=outf.createVariable('BC_'+str(var),'float32',('z','y','x'))
    outvar[:]=BC_Properties[var]

outf.Title='Brooks and Corey Soil Properties on the British National Grid'
outf.datasource1='Cranfield Soil data for England and Wales'
outf.datasource2='Scotland Soils data for Scotland'
outf.method='Created using pedotransfer functions'
outf.creator='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
outf.organisation='Centre of Ecology and Hydrology, Wallingford'

outf.close()
#######################################################################################

#######################################################################################
# Create boolean array of top/surface mask
SUB_SOIL=np.ones_like(COMP_data['clay'].data,dtype=bool)
SUB_SOIL[0,:]=False
# Compute the Van Genuchten Soil Properties
VG_Properties=VG.get_VG_soil_props_from_comp(COMP_data['clay'],           \
                                             COMP_data['sand'],           \
                                             COMP_data['silt'],           \
                                             OC_PC=COMP_data['org_carb'], \
                                             SUB_SOIL=SUB_SOIL,           \
                                             RETURN_SOIL_TEXTURE=True,    \
                                             )


VG_Properties['hcap']=TH.hcap(COMP_data['clay'],              \
                              COMP_data['sand'],              \
                              COMP_data['silt'],              \
                              sm_sat=VG_Properties['sm_sat'], \
                              l_vg=True                       )

VG_Properties['hcon']=TH.hcon_Farouki(COMP_data['clay'],      \
                                      COMP_data['sand'],      \
                                      COMP_data['silt'],      \
                                      sm_sat=VG_Properties['sm_sat'],\
                                      l_vg=True             )



VG_Properties['Soil_Texture']=VG_Properties['Soil_Texture'].astype(float)

# Output VG data
print('Outputting VG data to: '+OUT_VG_file)
outf=nc.Dataset(OUT_VG_file,'w')

# create dimensions
for dim in COMPinf.dimensions:
    outf.createDimension(str(dim),len(COMPinf.dimensions[dim]))
    outvar=outf.createVariable(str(dim),'float32',str(dim))
    for att in COMPinf.variables[dim].ncattrs():
        outvar.setncattr(str(att),COMPinf.variables[dim].getncattr(str(att)))
    outvar[:]=COMPinf.variables[dim][:]

for var in VG_Properties:
    outvar=outf.createVariable('VG_'+str(var),'float32',('z','y','x'))
    outvar[:]=VG_Properties[var]

outf.Title='Van Genuchten Soil Properties on the British National Grid'
outf.datasource1='Cranfield Soil data for England and Wales'
outf.datasource2='Scotland Soils data for Scotland'
outf.method='Created using pedotransfer functions'
outf.creator='Edward Comyn-Platt (edwcom@ceh.ac.uk)'
outf.organisation='Centre of Ecology and Hydrology, Wallingford'

outf.close()

#######################################################################################

COMPinf.close()


