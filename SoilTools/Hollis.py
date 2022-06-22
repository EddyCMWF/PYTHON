##############################################################
#
# Python module for Hollis PedoTransfer Functions
#    For Bulk Density
#   See J.M. Hollis 2012
#  
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np


##############################################################
# Function: Cult_TopSoil_BD
# Purpose: Calculate the Bulk density for Cultivate Top Soil
#            from the soil compistion. (Hollis, 2012)
# 
# Input: ORG_CARB_PC = Percentage organic carbon
#        SAND_PC = Percentage sand
#        CLAY_PC = Percentage clay
# Output: BULK_DENSITY = Bulk Density 
##############################################################
def Cult_TopSoil_BD(ORG_CARB_PC,SAND_PC,CLAY_PC):

    return    0.80806             \
           + (0.823844*np.exp(-0.27993*ORG_CARB_PC)) \
           + (0.0014065*SAND_PC)  \
           - (0.0010299*CLAY_PC)   


##############################################################
# Function: Comp_SubSoil_BD
# Purpose: Calculate the Bulk density for Compact Sub-Soil
#            from the soil compistion. (Hollis, 2012)
# 
# Input: ORG_CARB_PC = Percentage organic carbon
#        DEPTH = Depth of Soil Layer in cm
#        SAND_PC = Percentage sand
# Output: BULK_DENSITY = Bulk Density 
##############################################################
def Comp_SubSoil_BD(ORG_CARB_PC,DEPTH,SAND_PC):

    return    1.1257              \
           - (0.1140245*np.log(ORG_CARB_PC)) \
           + (0.0555*np.log(DEPTH))  \
           - (0.002248*SAND_PC)   


