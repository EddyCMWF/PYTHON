# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 11:58:35 2016

@author: Eward Comyn-Platt
"""


import netCDF4 as nc
from SoilTools import BrooksCorey as BC
from SoilTools import VanGenucthen as VG


Input_Dir='/users/eow/edwcom/GREENHOUSE/SOIL_PROPERTIES/datasets/'
Comp_File=Input_Dir+'Merged_Soil_Composition_CHESSgrid.nc'

Output_Dir=Input_Dir
OUT_BC_file=Output_Dir+'SoilProperties_BC_CHESSgrid.nc'
OUT_VG_file=Output_Dir+'SoilProperties_VG_CHESSgrid.nc'

# Read Soil Composition File
INPUT_vars=['sand','silt','clay','carbon']
COMP_data={}
COMPinf=nc.Dataset(Comp_File,'r')
for var in INPUT_vars:
    COMP_data[var]=COMPinf.variables[var][:]
COMPinf.close()


# Compute the Brooks and Corey Soil Properties
BC_Properties=BC.get_BC_soil_properties(COMP_data['clay'],\
                                        COMP_data['sand'],\
                                        COMP_data['silt'] )




# Compute the Brooks and Corey Soil Properties
VG_Properties=VG.get_VG_soil_properties(COMP_data['clay'],        \
                                        COMP_data['sand'],        \
                                        COMP_data['silt'],        \
                                        OC_PC=COMP_data['carbon'],\
                                        RETURN_SOIL_TEXTURE=True  )



print(BC_Properties)

import matplotlib.pyplot as plt


plotdata=COMP_data['sand']

plt.figure(figsize=(15,15))
for iSD in range(4):
    plt.subplot(2,2,iSD+1)
    plt.imshow(plotdata[iSD,:],origin='bottom')
    plt.colorbar()

plt.show()



