##############################################################
#
# Python module for Van Genuchten 
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np
from SoilTools.Soil_Global_Functions import *

# Van Genuchten Soil Limits:
clay_H=60
clay_M=35
clay_L=18
silt_H=50
sand_H=65
oc_H=10
oc_M=3

VG_Params_top = { 'ksat': [6.944e-02,1.396e-02,2.630e-03,2.870e-02,1.736e-02,9.259e-03], \
                  'sm_sat': [0.403,0.439,0.430,0.520,0.614,0.766], \
                  'sm_res': [0.025,0.010,0.010,0.010,0.010,0.010], \
                  'alpha': [3.830,3.140,0.830,3.670,2.650,1.300], \
                  'n': [1.3774,1.1804,1.2539,1.1012,1.1033,1.2039], \
                  }

VG_Params_sub = { 'ksat': [8.102e-02,1.245e-02,4.630e-03,9.838e-03,9.531e-03,9.259e-03], \
                  'sm_sat': [0.366,0.392,0.412,0.481,0.538,0.766], \
                  'sm_res': [0.025,0.010,0.010,0.010,0.010,0.010], \
                  'alpha': [4.300,2.490,0.820,1.980,1.680,1.300], \
                  'n': [1.5206,1.1689,1.2179,1.0861,1.0730,1.2039], \
                  }

##############################################################
# Function: Soil_Texture 
# Purpose: Classify the soil based on the sand/silt/clay and 
#            organic carbon percentages
#          This has been done in the same order as ELR function
#            who extracted from the MO CAP
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
#        SILT_PC = Percentage silt
# Optional Input:
#        OC_PC   = Percentage organic carbon
#                   If not provide, organic content not considered
# Output: Soil_Text - 6 bit array
#               1 = S_CR - Coarse Soil
#               2 = S_MD - Medium Soil
#               3 = S_MF - Medium Fine Soil
#               4 = S_FI - Fine Soil
#               5 = S_VF - Very Fine Soil
#               6 = S_OR - Organic Soil
#              
##############################################################
def Soil_Texture(CLAY_PC,SAND_PC,SILT_PC,OC_PC=None):
    # Convert any input to np.arrays
    scalar=False
    if not hasattr(CLAY_PC, "__len__"):
        CLAY_PC=np.array([CLAY_PC])
        scalar=True  # Flag to convert back to scalar
    if not hasattr(SAND_PC, "__len__"):
        SAND_PC=np.array([SAND_PC])
    if not hasattr(SILT_PC, "__len__"):
        SILT_PC=np.array([SILT_PC])

    Soil_Text=np.zeros_like(CLAY_PC,dtype='byte')
    # Very Fine Soil:
    VF_index=   (CLAY_PC>=clay_H)
    Soil_Text[VF_index] = 5

    # Medium Fine Soil:
    FI_index=   (CLAY_PC<clay_H)  \
              & (CLAY_PC>=clay_M)
    Soil_Text[FI_index] = 4

    # Medium Soil
    MD_index =   (CLAY_PC<clay_M)  \
               & (CLAY_PC>=0.)     \
               & ((SAND_PC<sand_H)|(CLAY_PC>clay_L))
    Soil_Text[MD_index] = 2
    
    # Medium Fine Soil
    MF_index =   (SILT_PC>=silt_H)  \
               & (CLAY_PC<clay_M)
    Soil_Text[MF_index] = 3
    
    # Coarse Soil
    CR_index =   (SAND_PC>=sand_H)  \
               & (CLAY_PC<=clay_L)  
    Soil_Text[CR_index] = 1

    if OC_PC!=None:
        if not hasattr(OC_PC, "__len__"):
            OC_PC=np.array([OC_PC])

        # Organic Soil
        OR_index = (OC_PC>=oc_H)
        Soil_Text[OR_index]=6

        # Additional Medium Soil from Coarse soils with Organic Content > 3%
        MD_index_2 =   (OC_PC>=oc_M)  \
                     & (Soil_Text==1)
        Soil_Text[MD_index_2]=2

    if (Soil_Text==0).any():
        print('Van Genuchten Error: Unclassified soil in dataset')
   
    if scalar:
        Soil_Text=Soil_Text[0]

    return Soil_Text
    
##############################################################
# Function: VG_param_alloc 
# Purpose: Allocate the Van Genucther parameter values to a 
#           soil texture array or value
# 
# Input: param - VG parameter value [ksat,theta_r,theta_s,alpha,n]
#        SOIL_TEXTURE = Soil Texture, corresponding to:
#                           1 = S_CR - Coarse Soil
#                           2 = S_MD - Medium Soil
#                           3 = S_MF - Medium Fine Soil
#                           4 = S_FI - Fine Soil
#                           5 = S_VF - Very Fine Soil
#                           6 = S_OR - Organic Soil
# Optional Input: SUB_SOIL=False - Flag for the Sub Surface Soil Parameters
#                 fill_value=-1
# Output: data - Allocated Van Genuchten data
##############################################################
def VG_param_alloc(param,SOIL_TEXTURE,SUB_SOIL=False,fill_value=-1):
    # Convert any input to np.arrays
    scalar=False
    if not hasattr(SOIL_TEXTURE, "__len__"):
        SOIL_TEXTURE=np.array([SOIL_TEXTURE])
        scalar=True  # Flag to convert back to scalar
    data = np.zeros_like(SOIL_TEXTURE,dtype='float32')+fill_value
    mask = (SOIL_TEXTURE>=1)&(SOIL_TEXTURE<=6)
    index= SOIL_TEXTURE[mask].astype('int')-1
    # if SUB_SOIL flag is a boolean apply sub/top to entire array
    if type(SUB_SOIL)==bool:
        if SUB_SOIL:
            data[mask]=np.array(VG_Params_sub[param])[index]
        else:
            data[mask]=np.array(VG_Params_top[param])[index]
    else:
        # Assume SUb_SOIL it is a Mask, True=sub_soil; False=top_soil
        # Combine SUB_SOIL and SOIL_TEXTURE masks
        SUB_mask  = (mask)&(SUB_SOIL)
        SUB_index = SOIL_TEXTURE[SUB_mask].astype('int')-1
        TOP_mask  = (mask)&(SUB_SOIL==False)
        TOP_index = SOIL_TEXTURE[TOP_mask].astype('int')-1
        # allocate sub and top soil components independantly
        data[SUB_mask]=np.array(VG_Params_sub[param])[SUB_index]
        data[TOP_mask]=np.array(VG_Params_top[param])[TOP_index]

    if scalar:
        data=data[0]
    return data
   

##############################################################
# Function: sm_sat_rel
# Purpose: Calculate the relative soil moisture at saturation:
#            sm_sat_rel = sm_sat-sm_res
# 
# Input: sm_sat - VG_parameter sm_sat
#        sm_res - VG parameter sm_res 
# 
# Output: sm_sat_rel - relative soil moisture at satration point
##############################################################
def sm_sat_rel(sm_sat,sm_res):
    return sm_sat-sm_res

##############################################################
# Function: sm_crit
# Purpose: Calculate the soil moisture at the critical point 
# 
# Input: alpha - VG parameter alpha 
#        n - VG parameter n 
#        sm_sat - VG_parameter sm_sat
#        sm_res - VG parameter sm_res 
# 
# Optional Input:
#        Critical_Point=-33 (kPa) - Critical Point Pressure
#
# Output: sm_wilt - Soil Moisture at Wilting Point
##############################################################
def sm_crit(alpha,n,sm_sat,sm_res, \
               Critical_Point=-33.):
    # Convert Matric Potential (Pa) to a pressure head (m)
    H_c=PascalsLaw(Critical_Point)*1e3
    # Apply equation from CAP
    sm_crit = ( ( 1 + (alpha*H_c)**n) ** ((1/n)-1) ) \
             * (sm_sat-sm_res)

    return sm_crit

##############################################################
# Function: sm_wilt
# Purpose: Calculate the soil moisture at the wilting point 
# 
# Input: alpha - VG parameter alpha 
#        n - VG parameter n 
#        sm_sat - VG_parameter sm_sat
#        sm_res - VG parameter sm_res 
# 
# Optional Input: 
#        Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
#
# Output: sm_wilt - Soil Moisture at Critical Point
##############################################################
def sm_wilt(alpha,n,sm_sat,sm_res, \
             Perm_Wilt_Point=-1500.):
    # Convert Matric Potential (Pa) to a pressure head (m)
    H_w=PascalsLaw(Perm_Wilt_Point)*1e3
    # Apply equation from CAP
    sm_wilt = ( ( 1 + (alpha*H_w)**n) ** ((1/n)-1) ) \
             * (sm_sat-sm_res)

    return sm_wilt


##############################################################
# Function: get_VG_soil_props_from_comp
# Purpose: Fetch all Hydraulic Van Genuchten Soil Properties from the
#           sand, silt, clay [and OC] percentage maps
# 
# Input: CLAY_PC = Percentage clay
#        SAND_PC = Percentage sand
#        SILT_PC = Percentage silt
# 
# Optional Input:
#        OC_PC   = Percentage organic carbon
#                   If not provide, organic content not considered
#        SUB_SOIL=False - Flag for the Sub Surface Soil Parameters
#        fill_value=-1
#        Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
#        Critical_Point=-33 (kPa) - Critical Point Pressure
#
# Output: Properties - Dictionary containg the Soil Hydraulic 
#                       Properties required by JULES
##############################################################
def get_VG_soil_props_from_comp(CLAY_PC,SAND_PC,SILT_PC, \
                                OC_PC=None, SUB_SOIL=False,fill_value=-1,\
                                Critical_Point=-33.,Perm_Wilt_Point=-1500.,\
                                RETURN_SOIL_TEXTURE=False,\
                                ):

    # Get Soil Textures From Composition
    SOIL_TEXTURE = Soil_Texture(CLAY_PC,SAND_PC,SILT_PC,OC_PC=OC_PC)

    Raw_Properties={}
    for param in VG_Params_top:
        Raw_Properties[param]= VG_param_alloc(param,SOIL_TEXTURE,\
                              SUB_SOIL=SUB_SOIL,fill_value=fill_value)

    # Now calcluate the sm_sat_rel, sm_crit and sm_wilt
    Properties = { 'oneovernminusone':1./(Raw_Properties['n']-1.),     \
                   'oneoveralpha':1./Raw_Properties['alpha'],          \
                   'ksat':Raw_Properties['ksat'],                      \
                   'sm_sat':sm_sat_rel(Raw_Properties['sm_sat'],\
                                       Raw_Properties['sm_res']),     \
                   'sm_crit':sm_crit(Raw_Properties['alpha'],     \
                                     Raw_Properties['n'],         \
                                     Raw_Properties['sm_sat'],    \
                                     Raw_Properties['sm_res'],    \
                                     Critical_Point=Critical_Point),   \
                   'sm_wilt':sm_wilt(Raw_Properties['alpha'],       \
                                     Raw_Properties['n'],           \
                                     Raw_Properties['sm_sat'],      \
                                     Raw_Properties['sm_res'],      \
                                     Perm_Wilt_Point=Perm_Wilt_Point), \
                   }
    
    if RETURN_SOIL_TEXTURE:
        Properties['Soil_Texture']=SOIL_TEXTURE
    
    return Properties



##############################################################
# Function: get_VG_soil_props_from_text
# Purpose: Fetch all Hydraulic Van Genuchten Soil Properties from the
#           SOIL_TEXTURE array/value
# 
# Input: SOIL_TEXTURE = Soil Texture array/value
# 
# Optional Input:
#        SUB_SOIL=False - Flag for the Sub Surface Soil Parameters
#        fill_value=-1
#        Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
#        Critical_Point=-33 (kPa) - Critical Point Pressure
#
# Output: Properties - Dictionary containg the Soil Hydraulic 
#                       Properties required by JULES
##############################################################
def get_VG_soil_props_from_text(SOIL_TEXTURE, \
                                SUB_SOIL=False,fill_value=-1,\
                                Critical_Point=-33.,Perm_Wilt_Point=-1500.,\
                                ):

    Raw_Properties={}
    for param in VG_Params_top:
        Raw_Properties[param]= VG_param_alloc(param,SOIL_TEXTURE,\
                              SUB_SOIL=SUB_SOIL,fill_value=fill_value)

    # Now calcluate the sm_sat_rel, sm_crit and sm_wilt
    Properties = { 'oneovernminusone':1./(Raw_Properties['n']-1.),     \
                   'oneoveralpha':1./Raw_Properties['alpha'],          \
                   'ksat':Raw_Properties['ksat'],                      \
                   'sm_sat':sm_sat_rel(Raw_Properties['sm_sat'],\
                                       Raw_Properties['sm_res']),     \
                   'sm_crit':sm_crit(Raw_Properties['alpha'],     \
                                     Raw_Properties['n'],         \
                                     Raw_Properties['sm_sat'],    \
                                     Raw_Properties['sm_res'],    \
                                     Critical_Point=Critical_Point),   \
                   'sm_wilt':sm_wilt(Raw_Properties['alpha'],       \
                                     Raw_Properties['n'],           \
                                     Raw_Properties['sm_sat'],      \
                                     Raw_Properties['sm_res'],      \
                                     Perm_Wilt_Point=Perm_Wilt_Point), \
                   }
    
    return Properties


