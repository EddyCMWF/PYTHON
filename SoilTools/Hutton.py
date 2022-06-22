##############################################################
#
# Python module for James Hutton Institute PedoTransfer Functions
#   See A. Lilly 2013
#  
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
#import numpy as np


##############################################################
# Function: bulk_density
# Purpose: Calculate the bulk_density from the soil compistion. 
#           (A. Lilly, 2013)
# 
# Input: ORG_CARB_PC = Percentage organic carbon
#        CLAY_PC = Percentage clay
#        SILT_PC = Percentage silt
# Optional In: ORG_CARB_MAX = 30, Maximum value for organic carbon
#                                 to avoid negative bulk densities
# Output: BULK_DENSITY = Bulk Density 
##############################################################
def bulk_density(ORG_CARB_PC,CLAY_PC,SILT_PC,\
                 ORG_CARB_MAX=15.):
    ORG_CARB_PC_wk = ORG_CARB_PC.copy()
    ORG_CARB_PC_wk[ORG_CARB_PC_wk>ORG_CARB_MAX]=ORG_CARB_MAX

    return    1.5031              \
           - (0.0418*ORG_CARB_PC) \
           + (0.01489*CLAY_PC)    \
           - (0.00523*SILT_PC)   


