#############################################################################
#
#   JULES photosynthesis functions.
#   
#   Owner: Edward Comyn-Platt
#   
#   Pupose: For sensitivity studies of the JULES photosynthesis behaviour 
#
#   Description: Based on Clark et al JULES model description. All equation numbers
#                 Refer to this paper.
#
#############################################################################
import numpy as np

#########################################################################
# Function: Qfunc
# Author: Edward Comyn-Platt
# Purpose: Simple Q Function
# Inputs: Q - Q Factor
#         T - Temperature (from reference point) (degrees C)
# Optional Inputs - Trate=10 - for change in rate of the Q function (degrees C)
# Outputs: fT - Q10 function for given temperature
#########################################################################
def Qfunc(Q,T,Trate=10):
    return Q**(Trate*T)

#########################################################################
# Function: Vcmax25
# Author: Edward Comyn-Platt
# Purpose: Calculation of Vcmax25 based on leaf nitrogen concentration 
#           (equation 6)
# Inputs: ne - constant relating leaf nitrogen to rubisco carboxylation capacity
#                (mol CO2 m^-2 s^-1 kg C kg^-1 N)
#         nl - Leaf nitrogen concentration (kg N kg^-1 C)
# Outputs: Vcmax25 - The maximun rate of Rubisco carboxylation at 25 degrees
#                    [sort of] (mol CO2 m^-2 s^-1)
#########################################################################
def Vcmax25(ne,nl): 
    return ne*nl

#########################################################################
# Function: TAUfunc
# Author: Edward Comyn-Platt
# Purpose: Calculation of tau based on canopy (leaf) temperature 
#           (equation 8)
# Input: Tc - Canopy(/leaf) temperature (degrees C)
# Optional Input: Q10=0.57 - Q10 factor for tau (-)
# Outputs: tau - Rubisco specificity (-)
#########################################################################
def TAUfunc(Tc,Q10=0.57): 
    return 2600.*Qfunc(Q10,Tc)

#########################################################################
# Function: KCfunc
# Author: Edward Comyn-Platt
# Purpose: Calculation of Kc (Michaelis-Menten parameter for CO2) 
#                   based on canopy (leaf) temperature 
#           (equation 9)
# Input: Tc - Canopy(/leaf) temperature (degrees C)
# Optional Input: Q10=2.1 - Q10 factor for Kc (-)
# Outputs: Kc - Michaelis-Menten parameter for CO2 (Pa)
#########################################################################
def KCfunc(Tc,Q10=2.1): 
    return 30.*Qfunc(Q10,Tc)

#########################################################################
# Function: KOfunc
# Author: Edward Comyn-Platt
# Purpose: Calculation of Ko (Michaelis-Menten parameter for O2) 
#                   based on canopy (leaf) temperature 
#           (equation 9)
# Input: Tc - Canopy(/leaf) temperature (degrees C)
# Optional Input: Q10=1.2 - Q10 factor for Ko (-)
# Outputs: Ko - Michaelis-Menten parameter for O2 (Pa)
#########################################################################
def KOfunc(Tc,Q10=1.2): 
    return 3e4*Qfunc(Q10,Tc)

#########################################################################
# Function: OAfunc
# Author: Edward Comyn-Platt
# Purpose: Calculation of Oa (Atmospheric pressure of O2) 
#           (From JULES source code, sf_stom_jls.F90 in vn4.5)
# Input: Psurf - Surface pressure (Pa)
# Optional Input: o2=0.23 - Atmopheric concentration of O2 (kg O2 kg^-1 Air)
#                 epo2=1.106 -Ratio of molecular weights of O2 and dry air (-)
# Outputs: Oa - Atmorpheric pressure of O2 (Pa)
#########################################################################
def OAfunc(Psurf,o2=0.23,epo2=1.106):
    return (Psurf**-1) * (o2/epo2)

#########################################################################
# Function: CO2_compensation_point
# Author: Edward Comyn-Platt
# Purpose: Calculation of CO2_compensate based on Rubisco specificity and 
#                  the partial pressure of Oxyger
#           (equation 7)
# Inputs: Oa  - the partial pressure of atmospheric oxygen (Pa)
#         tau -  Rubisco specificity (-)
#         l_C3 - logical stating whether or not plant is a C3 grass (-)
# Outputs: CO2_compensate - CO2 compensation point in the absence 
#                            of mitochondrial respiration (Pa)
#########################################################################
def CO2_compensation_point(Oa,tau,l_C3):
    # assume it is an array
    try:
        CO2_compensate=np.zeros_like(Oa)
        CO2_compensate[l_C3==1] = Oa/(2*tau)
    except:
        if l_C3==1:
            CO2_compensate=Oa/(2*tau)
        else:
            CO2_compensate=0
    return CO2_compensate
    
##########################################################################
# Function: Vcmax_func
# Author: Edward Comyn-Platt
# Purpose:  Calculate Vcmax as a function of temperature (equation 4)
# Input: Vcmax25 - the maximum rate of carboxylation of Rubisco at 
#                   25 degrees C [sort of]  (mol CO2 m^-2 s^-1)
#        Tc - Canopy(/leaf) temperature (degrees C)
#        Tupp - upper temperature parameter (degrees C)
#        Tlow - lower temperature parameter (degrees C)
# Optional Input: Q10=2.0 - Q10 factor for photosynthesis (-)
# Output: Vcmax - the maximum rate of carboxylation of Rubisco (mol CO2 m^-2 s^-1) 
##########################################################################
def Vcmax_func(Vcmax25,Tc,Tupp,Tlow,Q10=2.):
    fT = Qfunc(Q10,Tc-25)
    numerator = Vcmax25*fT
    denominator = (1+np.exp(0.3*(Tc-Tupp))) * \
                  (1+np.exp(0.3*(Tlow-Tc)))   
    return numerator/denominator

##########################################################################
# Function: rubisco_limited_rate
# Author: Edward Comyn-Platt
# Purpose:  Calculate the rubisco limited rate of photosynthesis - equation 1
# Input: Vcmax - the maximum rate of carboxylation of Rubisco (mol CO2 m^-2 s^-1)
#        CO2_leaf - the leaf internal CO2 partial pressure (Pa)
#        CO2_compensate - the CO2 compensation point in the absence of 
#                            mitochondrial respiration (Pa)
#        Oa   - the partial pressure of atmospheric oxygen (Pa)
#        Kc   - Michaelis-Menten parameters for CO2 (Pa)
#        Ko   - Michaelis-Menten parameters for O2 (Pa) 
#        l_C3 - logical stating whether or not plant is a C3 grass (-)
# Output: W_c - the Rubisco limited rate photosynthesis (mol CO2 m^-2 s^-1)
##########################################################################
def rubisco_limited_rate(Vcmax,CO2_leaf,CO2_compensate,Oa,Kc,Ko,l_C3):
    # Assume input data is an array:
    try:
        # create array of ones same shape as Vcmax
        factor = np.ones_like(Vcmax)
        # For elements where l_C3==1 calculate the scale factor according to equation 1
        factor[l_C3==1] = ( CO2_leaf-CO2_compensate )/  \
                          ( CO2_leaf+ (Kc*(1+ (Oa/Ko))) )
    except:
        # if a scalar value then:
        if (l_C3==1):
            # if l_C3==1 calculate the scale factor according to equation 1
            factor=( CO2_leaf-CO2_compensate )/  \
                   ( CO2_leaf+ (Kc*(1+ (Oa/Ko))) )
        else:
            factor=1
    # Equation 1:
    W_c = Vcmax * factor
    return W_c

##########################################################################
# Function: light_limited_rate
# Author: Edward Comyn-Platt
# Purpose:  Calculate the light limited rate of photosynthesis - equation 2
# Input: Ipar - the maximum rate of carboxylation of Rubisco (mol PAR m^-2 s^-1)
#        CO2_leaf - the leaf internal CO2 partial pressure (Pa)
#        CO2_compensate - the CO2 compensation point in the absence of 
#                          mitochondrial respiration (Pa)
#        alpha - Quantum efficiency of photosynthesis (mol CO2 mol^-1 PAR)
#        omega - Leaf scattering coefficient for PAR (-)
#        l_C3 - logical stating whether or not plant is a C3 grass (-)
# Output: W_l - the Light limited rate photosynthesis (mol CO2 m^-2 s^-1)
##########################################################################
def light_limited_rate(Ipar,CO2_leaf,CO2_compensate,alpha,omega,l_C3):
    # Assume input data is an array:
    try:
        # create array of ones same shape as Vcmax
        factor = np.ones_like(Ipar)
        # For elements where l_C3==1 calculate the scale factor according to equation 2
        factor[l_C3==1] = ( CO2_leaf-CO2_compensate )/  \
                          ( CO2_leaf+(2*CO2_compensate) )
    except:
        # if a scalar value then:
        if (l_C3==1):
            #if l_C3==1 calculate the scale factor according to equation 2
            factor=( CO2_leaf-CO2_compensate )/  \
                   ( CO2_leaf+(2*CO2_compensate) )
        else:
            factor=1
    # Equation 2:
    W_l = Ipar * alpha * (1-omega) * factor
    return W_l

##########################################################################
# Function: trans_limited_rate
# Author: Edward Comyn-Platt
# Purpose:  Calculate the transport/PEPCarboxylate limited rate of photosynthesis - equation 3
# Input: Vcmax - the maximum rate of carboxylation of Rubisco (mol CO2 m^-2 s^-1)
#        CO2_leaf - the leaf internal CO2 partial pressure (Pa)
#        Psurf - The surface air pressure (Pa)
#        l_C3 - logical stating whether or not plant is a C3 grass (-)
# Output: W_e - the transport/PEPCarboxylate limited rate photosynthesis (mol CO2 m^-2 s^-1)
##########################################################################
def light_limited_rate(Vcmax,CO2_leaf,Psurf,l_C3):
    # Assume input data is an array:
    try:
        # create array of zeros same shape as Vcmax
        W_e = np.ones_like(Vcmax)
        # Apply equation 3 depending on whether or not it is a C3 or C4 plant
        W_e[l_C3==1] = W_e*0.5
        W_e[l_C3==0] = W_e*2e4 * (CO2_leaf/Psurf)
    except:
        # if a scalar value then use an if statement:
        if (l_C3==1):
            W_e = W_e*0.5
        else:
            W_e = W_e*2e4 * (CO2_leaf/Psurf)
    return W_e

##########################################################################
# Function: trans_limited_rate
# Author: Edward Comyn-Platt
# Purpose:  Calculate the rate of gross photosynthesis (W) as a smoothed
#            minimum of the three limiting rates (equations 11 and 12)
# Input: W_c - rubisco limited rate of photosynthesis (mol CO2 m^-2 s^-1)
#        W_l - light limited rate of photosynthesis (mol CO2 m^-2 s^-1)
#        W_e - transport/PEPCarboxylate limited rate of photosynthesis (mol CO2 m^-2 s^-1)
# Optional Input:  beta1=0.83 - colimitation coefficient 1 (-)
#                  beta2=0.93 - colimitation coefficient 2 (-)
# Output: W - the rate of gross photosynthesis (mol CO2 m^-2 s^-1)
##########################################################################
def gross_photosynthesis_rate(W_c,W_l,W_e,\
                               beta1=0.83,beta2=0.93):


    # first find the smoothed minimum of the light and rubisco limiting rates
    # equation 11:
    # beta1*Wp**2 + Wp*(Wc+Wl) + (Wc*Wl)
    # Quadratic coefficients:
    a = np.zeros_like(W_c)+beta1
    b = -(W_c+W_l)
    c = W_c*W_l
    Wp= ( -b - np.sqrt( (b**2)-(4*a*c) ) )  \
            / (2*a)
    
    
    # first find the smoothed minimum of Wp and the trasnport limiting rate
    # equation 11:
    # beta2*W**2 + W*(Wp+We) + (Wp*We)
    # Quadratic coefficients:
    a = np.zeros_like(W_c)+beta2
    b = -(W_p+W_e)
    c = W_p*W_e
    W = ( -b - np.sqrt( (b**2)-(4*a*c) ) )  \
            / (2*a)
    
    return W

#########################################################################
# Function: leaf_dark_resp
# Author: Edward Comyn-Platt
# Purpose: Calculation of Leaf Dark Respiration Rd (equation 13) 
# Input: Vcmax - the maximum rate of carboxylation of Rubisco (mol CO2 m^-2 s^-1)
#        f_dr  - dark respiration coefficient (-)
# Outputs: Rd - Leaf Dark Respiration (mol CO2 m^-2 s^-1)
#########################################################################
def leaf_dark_resp(VCmax,f_dr):
    return Vcmax*f_dr

#########################################################################
# Function: leaf_photosynthesis
# Author: Edward Comyn-Platt
# Purpose: Calculation of Leaf Level Photosynthesis (equation 14) 
# Input: W - the rate of gross photosynthesis (mol CO2 m^-2 s^-1)
#        Rd - Leaf Dark Respiration (mol CO2 m^-2 s^-1)
# Outputs: Ap - Leaf Photosynthesis (mol CO2 m^-2 s^-1)
#########################################################################
def leaf_photosynthesis(W,Rd):
    return W - Rd

#########################################################################
# Function: leaf_level_photosynthesis
# Author: Edward Comyn-Platt
# Purpose: Calculation of Leaf Level Photosynthesis (equation 15)
# Input: Ap - Leaf Photosynthesis (mol CO2 m^-2 s^-1)
#        beta - soil moisture stress factor (fsmc in source code)
# Outputs: Al - Leaf Level Photosynthesis (mol CO2 m^-2 s^-1)
#########################################################################
def leaf_level_photosynthesis(Ap,beta):
    return Ap*beta 

#########################################################################
# Function: beta_func
# Author: Edward Comyn-Platt
# Purpose: Calculate beta - the soil moisture stress factor (FSMC)
# Input: theta - mean soil moisture concentration in the root zone
#        theta_w - wilting point soil moisture concentration
#        theta_c - critical point soil moisture concentration
# Output: beta - the soil moiture stress factor
#########################################################################
def beta_func(theta,theta_w,theta_c):
    beta = (theta-theta_w)/(theta_c-theta_w)
    # First try assumption that data is in an array and limit between 1 and 0
    try:
       beta[beta>1.]=1. 
       beta[beta<0.]=0.
    except:
        # if a scalar use an if statement
        if beta>1.
            beta=1.
        elif beta<0.:
            beta=0.

    return beta










