##############################################################
#
# Python module for Soil Property Functions
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np

rho=1000.  # Volumetric density of water (1000 kg m^-3)
g=9.81     # Gravitational Acceleration at Earth's surface

def PascalsLaw(psi):
    h=-psi/(rho*g)
    return h
def InvPascalsLaw(h):
    psi=-h*rho*g
    return psi

##############################################################
# Class: BrooksCorey
# 
# Purpose: Class of Brooks and Corey Functions
# 
##############################################################
class BrooksCorey:

    ##############################################################
    # Function: bexp
    # Purpose: Calculate the Brooks and Corey b exponent 
    #            from the soil compistion. (Marthews, 2014)
    # 
    # Input: CLAY_PC = Percentage clay
    #        SAND_PC = Percentage sand
    # Output: b = b exponent
    ##############################################################
    def bexp(self,CLAY_PC,SAND_PC):
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
    def lambd(self,CLAY_PC,SAND_PC):
        lamb = 1./self.bexp(CLAY_PC,SAND_PC) 
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
    def sathh(self,CLAY_PC,SAND_PC,exponential=False):
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
    def psisat(self,CLAY_PC,SAND_PC,exponential=False):
        psisat = PascalsLaw(self.sathh(CLAY_PC,SAND_PC,exponential=exponential))
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
    def satcon(self,CLAY_PC,SAND_PC,exponential=False):
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
    def satconMO(self,CLAY_PC,SAND_PC,exponential=False):
        if exponential:
            satcon =  ( np.exp(-5.55-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
        else:
            satcon =  ( 10**(-2.75-(0.0064*CLAY_PC)+(0.0126*SAND_PC)) )
        return satcon
    
    ##############################################################
    # Function: sm_wilt
    # Purpose: Calculate the Brooks and Corey Soil Moisture 
    #            at Wilting point in m^3 m^-3,
    # Input: CLAY_PC = Percentage clay
    #        SAND_PC = Percentage sand
    # Optional Input: Perm_Wilt_Point=-1500 (kPa) - Permament wilting point
    # Output: sm_wilt - soil moisture at wilting point
    ##############################################################
    def sm_wilt(self,CLAY_PC,SAND_PC,Perm_Wilt_Point=-1500.):
        sm_wilt =  self.satcon(CLAY_PC,SAND_PC) * \
                   (self.sathh(CLAY_PC,SAND_PC)/PascalsLaw(Perm_Wilt_Point*-1e3) ) ** \
                       self.bexp(CLAY_PC,SAND_PC) 
        return sm_wilt
    
    
    
    #V_WILT(I)=V_SAT(I)*(SATHH(I)/HEAD_W)**(1.0/B(I))
