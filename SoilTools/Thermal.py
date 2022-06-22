##############################################################
#
# Python module for Soil Thermal Property Functions
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np
from SoilTools.Soil_Global_Functions import *
from SoilTools import BrooksCorey as BC
from SoilTools import VanGenuchten as VG

##############################################################
# Function: hcap
# Purpose: Calculate the Heat Capacity 
#            from the soil compistion or a precalculated sm_sat
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
#        SILT_PC = Percentage silt
# Optional Input: sm_sat=None - Volumetric Soil Moisture at saturation
#                                If not provide will calculate using BC/VG equations
#                 l_vg=False - logical to use Van Genuchten equations to calculate
#                               sm_sat, if not provided.
#                 CAP_[CLAY/SAND/SILT]=cap_[clay/sand/silt] 
#                      - Heat Capacities of clay/sand/silt, default is the value 
#                          provided in SOIL_GLOBAL_FUNCTIONS
# Output: hcap = heat capacity 
##############################################################
def hcap(CLAY_PC,SAND_PC,SILT_PC,\
            sm_sat=None,l_vg=False, \
            CAP_CLAY=cap_clay,CAP_SAND=cap_sand,CAP_SILT=cap_silt):
    if sm_sat==None:
        if l_vg:
            sm_sat=VG.sm_sat_rel(CLAY_PC,SAND_PC)
        else:
            sm_sat=BC.sm_sat(CLAY_PC,SAND_PC)

    hcap= (1.-sm_sat) * ( (CLAY_PC*CAP_CLAY)+  \
                          (SAND_PC*CAP_SAND)+  \
                          (SILT_PC*CAP_SILT) ) \
                *0.01
    return hcap


##############################################################
# Function: hcap_Johansen
# Purpose: Calculate the Heat Capacity for a mean dry heat capacity 
#               Given by Johansen 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: sm_sat=None - Volumetric Soil Moisture at saturation
#                                If not provide will calculate using BC/VG equations
#                 CLAY_PC=None - Percentage Clay used to calculate sm_sat if not provided
#                 SAND_PC=None - Percentage Sand - Used to calculate sm_sat if not provided
#                 l_vg=False - logical to use Van Genuchten equations to calculate
#                               sm_sat, if not provided.
#                 CAP_[CLAY/SAND/SILT]=cap_[clay/sand/silt] 
#                      - Heat Capacities of clay/sand/silt, default is the value 
#                          provided in SOIL_GLOBAL_FUNCTIONS
# Output: hcap = heat capacity 
##############################################################
def hcap_Johansen(  sm_sat=None,l_vg=False, \
                    CLAY_PC=None, SAND_PC=None, \
                    CAP_SOIL=Johansen_cap_soil ):
    if sm_sat==None:
        if (CLAY_PC==None)|(SAND_PC==None):
            print('Error- Must Provide sm_sat or (CLAY_PC and SAND_PC) to hcap_Johansen')
            return np.nan
        if l_vg:
            sm_sat=VG.sm_sat_rel(CLAY_PC,SAND_PC)
        else:
            sm_sat=BC.sm_sat(CLAY_PC,SAND_PC)

    hcap= (1.-sm_sat) * CAP_SOIL
    
    return hcap


##############################################################
# Function: hcon_Farouki
# Purpose: Calculate the Farouki Thermal Conductivity 
#            from the soil compistion 
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
#        SILT_PC = Percentage silt
# Optional Input: sm_sat=None - A precalculated Volumetric Soil Moisture at saturation
#                                If not provide will calculate using BC/VG equations
#                 l_vg=False - logical to use Van Genuchten equations to calculate
#                               sm_sat, if not provided.
#                 CAP_[AIR/CLAY/SAND/SILT]=cap_[air/clay/sand/silt] 
#                      - Heat Capacities of air/clay/sand/silt, default is the value 
#                          provided in SOIL_GLOBAL_FUNCTIONS
# Output: hcon = thermal conductivity
##############################################################
def hcon_Farouki(CLAY_PC,SAND_PC,SILT_PC,\
                 sm_sat=None,l_vg=False, \
                 CON_AIR=con_air,CON_CLAY=con_clay,CON_SAND=con_sand,CON_SILT=con_silt):
    if sm_sat==None:
        if l_vg:
            sm_sat=VG.sm_sat_rel(CLAY_PC,SAND_PC)
        else:
            sm_sat=BC.sm_sat(CLAY_PC,SAND_PC)

    hcon=   (CON_AIR**sm_sat)                \
          * (CON_CLAY**(CLAY_PC*0.01*sm_sat)) \
          * (CON_SAND**(SAND_PC*0.01*sm_sat)) \
          * (CON_SILT**(SILT_PC*0.01*sm_sat))

    return hcon

##############################################################
# Function: hcon_LU
# Purpose: Calculate the LU Thermal Conductivity 
#            from the soil compistion or sm_sat
# Optional Input: sm_sat=None - A precalculated Volumetric Soil Moisture at saturation
#                                If not provide will calculate using BC/VG equations
#                 CLAY_PC=None - Percentage Clay used to calculate sm_sat if not provided
#                 SAND_PC=None - Percentage Sand - Used to calculate sm_sat if not provided
#                 l_vg=False - logical to use Van Genuchten equations to calculate
#                               sm_sat, if not provided.
# Output: hcon = thermal conductivity
##############################################################
def hcon_LU(sm_sat=None,CLAY_PC=None,SAND_PC=None,l_vg=False):
    if sm_sat==None:
        if (CLAY_PC==None)|(SAND_PC==None):
            print('Error- Must Provide sm_sat or (CLAY_PC and SAND_PC) to hcap_Johansen')
        if l_vg:
            sm_sat=VG.sm_sat_rel(CLAY_PC,SAND_PC)
        else:
            sm_sat=BC.sm_sat(CLAY_PC,SAND_PC)
    hcon= ( -0.56 * sm_sat ) +0.51 
    return hcon


