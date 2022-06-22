##############################################################
#
# Python module for Soil Property Functions
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
#import numpy as np

#Constants
rho=1000.  # Volumetric density of water (1000 kg m^-3)
g=9.81     # Gravitational Acceleration at Earth's surface

# heat capacities of sand silt and clay from table 4.1 of Frozen Earth
cap_clay=2.373e6
cap_sand=2.133e6
cap_silt=2.133e6

# thermal conductivities of sand silt clay and air from table 4.1 of Frozen Earth
con_sand=1.57
con_silt=1.57
con_clay=1.16
con_air =0.025

# Johansen mean soil heat capacity
Johansen_cap_soil=1.942e6

#Functions
def PascalsLaw(psi):
    h=-psi/(rho*g)
    return h
def InvPascalsLaw(h):
    psi=-h*rho*g
    return psi

