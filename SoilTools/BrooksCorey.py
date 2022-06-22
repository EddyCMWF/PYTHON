##############################################################
#
# Python module for Brooks and Corey Functions
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np
from SoilTools.Soil_Global_Functions import *

##############################################################
# Function: bexp
# Purpose: Calculate the Brooks and Corey b exponent 
#            from the soil compistion. (Marthews, 2014)
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Output: b = b exponent
##############################################################
def bexp(CLAY_PC,SAND_PC):
    b=3.1 + (0.157*CLAY_PC) - (0.003*SAND_PC)
    return b
    
##############################################################
# Function: lambd
# Purpose: Calculate the Brooks and Corey lambda exponent 
#            from the soil compistion. (Marthews, 2014)
#             (lambda=1/b)
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Output: lamb = lambda exponent
##############################################################
def lambd(CLAY_PC,SAND_PC):
    lamb = 1./bexp(CLAY_PC,SAND_PC) 
    return lamb
    
##############################################################
# Function: sathh
# Purpose: Calculate the Brooks and Corey Soil Matric Suction
#            at Saturation in meters, (h_e in Marthews, 2014)
#             
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: exponential=False - Flag to use base e 
#                                       instead of base 10
# Output: sathh = soil matric suctionat saturation in meters
##############################################################
def sathh(CLAY_PC,SAND_PC,exponential=False):
    if exponential:
        sathh = 0.01 * ( np.exp(2.17-(0.0063*CLAY_PC)-(0.0158*SAND_PC)) )
    else:
        sathh = 0.01 * ( 10.**(2.17-(0.0063*CLAY_PC)-(0.0158*SAND_PC)) )
    return sathh
        
##############################################################
# Function: psisat
# Purpose: Calculate the Brooks and Corey Soil Matric Suction
#            at Saturation in Pa, (psi_e in Marthews, 2014)
#              (psisat=sathh/(-rho*g))
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Inpiut: exponential=False - Flag to use base e 
#                                       instead of base 10
# Output: psisat = soil matric suctionat saturation in Pascals
##############################################################
def psisat(CLAY_PC,SAND_PC,exponential=False):
    psisat = InvPascalsLaw(sathh(CLAY_PC,SAND_PC,exponential=exponential))
    return psisat
        
##############################################################
# Function: satcon
# Purpose: Calculate the Brooks and Corey Hydraulic Conductivity
#            at Saturation in kg m^-2 s^-1, (k_sat in Marthews, 2014)
#             
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: exponential=False - Flag to use base e 
#                                       instead of base 10
# Output: satcon - hydraulic conductivity at saturation
##############################################################
def satcon(CLAY_PC,SAND_PC,exponential=False):
    if exponential:
        satcon =  (25.4/3600.)* ( np.exp(-0.6-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
    else:
        satcon =  (25.4/3600.)* ( 10**(-0.6-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
    return satcon
        
##############################################################
# Function: satconMO
# Purpose: Calculate the Brooks and Corey Hydraulic Conductivity
#            at Saturation using the MO formulation in kg m^-2 s^-1,
#             
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: exponential=False - Flag to use base e 
#                                       instead of base 10
# Output: satcon - hydraulic conductivity at saturation
##############################################################
def satconMO(CLAY_PC,SAND_PC,exponential=False):
    if exponential:
        satcon =  ( np.exp(-5.55-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
    else:
        satcon =  ( 10**(-2.75-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
    return satcon
    

##############################################################
# Function: sm_sat
# Purpose: Calculate the Brooks and Corey Soil Moisture 
#            at Saturation m^3 m^-3,
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: 
# Output: sm_sat - soil moisture at saturation
##############################################################
def sm_sat(CLAY_PC,SAND_PC):
    sm_sat= 0.01 * ( 50.5 - (0.037*CLAY_PC) - (0.142*SAND_PC) )
    return sm_sat


##############################################################
# Function: sm_wilt
# Purpose: Calculate the Brooks and Corey Soil Moisture 
#            at Wilting point in m^3 m^-3,
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
#                 SM_Residual = 0. (m^3/m^-3) - Soil Moisture Residual 
#                 exponential = False  - compute psisat with base e
# Output: sm_wilt - soil moisture at wilting point
##############################################################
def sm_wilt(CLAY_PC,SAND_PC,\
            Perm_Wilt_Point=-1500.,SM_Residual=0.,exponential=False):
    sm_wilt =  (sm_sat(CLAY_PC,SAND_PC)-SM_Residual) \
             * ( (psisat(CLAY_PC,SAND_PC,exponential=exponential)/(Perm_Wilt_Point*1e3) )\
                 ** lambd(CLAY_PC,SAND_PC) ) \
             + SM_Residual 
    return sm_wilt
    
##############################################################
# Function: sm_crit
# Purpose: Calculate the Brooks and Corey Soil Moisture 
#            at Critical point in m^3 m^-3,
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: Critical_Point=-33 (kPa) - Critical Point Pressure
#                 SM_Residual = 0. (m^3/m^-3) - Soil Moisture Residual 
#                 exponential = False  - compute psisat with base e
# Output: sm_crit - soil moisture at critical point
##############################################################
def sm_crit(CLAY_PC,SAND_PC,\
            Critical_Point=-33.,SM_Residual=0.,exponential=False):
    sm_crit =  (sm_sat(CLAY_PC,SAND_PC)-SM_Residual) \
             * ( (psisat(CLAY_PC,SAND_PC,exponential=exponential)/(Critical_Point*1e3) )\
                 ** lambd(CLAY_PC,SAND_PC) ) \
             + SM_Residual 
    return sm_crit

##############################################################
# Function: effsat_frompsi
# Purpose: Calculate the Brooks and Corey Effective Saturation
#           as a fraction  (Marthews, 2014)
# Input: psi = Soil Matric Potential
#        CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: exponential=False (Compute matric potential using base e)
# Output: effsat - effective saturation calculated using Matric Potential (psi)
##############################################################
def effsat_frompsi(psi,CLAY_PC,SAND_PC,exponential=False):
    effsat  = (psi/psisat(CLAY_PS,SAND_PC,exponential=exponential)) ** \
               lambd(CLAY_PC,SAND_PC) 
    return effsat
    
##############################################################
# Function: effsat_fromSM 
# Purpose: Calculate the Brooks and Corey Effective Saturation
#           as a fraction  (Marthews, 2014)
# Input: SM = Volumetric Soil Moisture (m^3/m^-3)
#        CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
# Optional Input: SM_Residual=0 - Soil Moisture Residual
# Output: effsat - effective saturation calculated using the SM
##############################################################
def effsat_fromSM(SM,CLAY_PC,SAND_PC,SM_Residual=0.):
    effsat  = (SM-SM_Residual) / \
              (sm_sat(CLAY_PC,SAND_PC)-SM_Residual) 
    return effsat
    
    
#########################################################################
# Function: get_BC_soil_properties
# Purpose: Fetch all Hydraulic Brooks and Corey Soil Properties from the
#           sand, silt, clay [and OC] percentage maps
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
#        SILT_PC = Percentage silt
# 
# Optional Input:
#        OC_PC   = Percentage organic carbon
#                   If not provide, organic content not considered
#        fill_value=-1
#        Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
#        Critical_Point=-33 (kPa) - Critical Point Pressure
#
# Output: Properties - Dictionary containg the Soil Hydraulic 
#                       Properties required by JULES
#######################################################################
def get_BC_soil_properties(CLAY_PC,SAND_PC,SILT_PC, \
                           exponential=False, SM_Residual=0.,\
                           Critical_Point=-33.,Perm_Wilt_Point=-1500.,\
                           ):

    Properties = { 'bexp':bexp(CLAY_PC,SAND_PC),                       \
                   'sathh':sathh(CLAY_PC,SAND_PC,       \
                                 exponential=exponential),             \
                   'satcon':satcon(CLAY_PC,SAND_PC,       \
                                   exponential=exponential),           \
                   'sm_sat':sm_sat(CLAY_PC,SAND_PC),                   \
                   'sm_crit':sm_crit(CLAY_PC,SAND_PC,              \
                                     Critical_Point=Critical_Point,\
                                     SM_Residual=SM_Residual,      \
                                     exponential=exponential),         \
                   'sm_wilt':sm_wilt(CLAY_PC,SAND_PC,                \
                                     Perm_Wilt_Point=Perm_Wilt_Point,\
                                     SM_Residual=SM_Residual,        \
                                     exponential=exponential),         \
                   }
    return Properties


    
    
