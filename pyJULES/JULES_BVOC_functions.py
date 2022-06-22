#
#
# Python module of JULES BVOC functions
#
# Edward Comyn-Platt, edwcom@ceh.ac.uk
# Centre for Ecology and Hydrology
# 2015
#
import numpy as np

###############################################################################
# Function: F_CO2
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the internal CO2 factor of BVOC emission
#               Pacifico 2010
#
# Required Inputs:
#  Ci - Internal leaf CO2 partial pressure
#         if multiple pft, then pft must be the first dimension
#
# Optional Inputs:
#  Ci_st - Internal leaf CO2 partial pressure under standard conditions
#                 Default = [33.46,33.46,34.26,29.98,34.26],JULES default)
#
# Output:
#  f_co2 - Internal CO2 factor of BVOC emission
# 
###############################################################################
def F_CO2( Ci,Ci_st=[33.46,33.46,34.26,29.98,34.26] ):
    
    ###########################################################################
    # 1. Convert all data to np.array
    ###################################
    Ci=np.array(Ci)
    Ci_st=np.array(Ci_st)
        
    ###########################################################################
    # 2. check dimensions
    ###########################
    # 2.1 nPFTs is set to number of gamma_v values
    nPFTs = len(Ci_st)
    #
    # 2.2 check PFTf has correct number of PFTs
    if (Ci.shape[0]!=nPFTs):
        print( 'ERROR: Ci and Ci_st have inconsistent nPFTs')
        print( 'Ci.shape = ', Ci.shape )
        print( 'Ci_st.shape = ',Ci_st.shape)
        print( 'Returning NaN')
        return np.nan
    #
    ##############################################
    # Compute the factor
    #    Ci/Ci_st
    #  reshpae is used to maintain the dimensions
    f_co2=(Ci.reshape(5,-1)/Ci_st[:,None]).reshape(Ci.shape) 
    ########################
    return f_co2


###############################################################################
# Function: F_T_isop
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the Temperature factor of isoprene emission
#               Pacifico 2010
#
# Required Inputs:
#  Tstar - Leaf/Surface temperature 
#
# Optional Inputs:
#  f_tmax - Empirical estimate that makes isoprene emission level off
#            at temperature close to 40 C (Arneth personal communication)
#             default = 2.3
#  atau - Scaling parameter (isoprene) (eq A4b Arneth et al., 2007)
#             default=0.1
#  T_ref - Reference temperature (Arneth 2007)
#             default=303.15
#
# Output:
#  f_T_isop - Temperature factor of Isoprene emissions
# 
###############################################################################
def F_T_isop( Tstar, f_tmax=2.3, atau=0.1, T_ref=303.15 ):
    
    ###########################################################################
    # 1. Convert all data to np.array
    ###################################
    t_star=np.array(t_star)
        
    ##############################################
    # Compute the factor
    #    Minimum ( f_tmax, exp[ atau* (tstar-t_ref) 
    f_t_isop = np.exp( atau * (Tstar-T_ref) )
    f_t_isop[f_t_isop<f_tmax]=f_tmax
    ########################
    return f_t_isop


###############################################################################
# Function: F_GPP_isop
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the GPP factor of isoprene emission
#               Pacifico 2010
#
# Required Inputs:
#  GPP - GPP 
#
# Optional Inputs:
#  GPP_st - GPP under standard conditions
#         Default = [1.29e-7,2.58e-8,2.07e-7,3.42e-7,1.68e-7], JULES default 
#
# Output:
#  f_gpp_isop - GPP factor of isoprene emission
# 
###############################################################################
def F_GPP_isop( GPP,GPP_st=[1.29e-7,2.58e-8,2.07e-7,3.42e-7,1.68e-7]):
    
    ###########################################################################
    # 1. Convert all data to np.array
    ###################################
    GPP=np.array(GPP)
    GPP_st=np.array(GPP_st)
        
    ###########################################################################
    # 2. check dimensions
    ###########################
    # 2.1 nPFTs is set to number of gamma_v values
    nPFTs = len(GPP_st)
    #
    # 2.2 check PFTf has correct number of PFTs
    if (GPP.shape[0]!=nPFTs):
        print( 'ERROR: GPP and GPP_st have inconsistent nPFTs')
        print( 'GPP.shape = ', GPP.shape )
        print( 'GPP_st.shape = ',GPP_st.shape)
        print( 'Returning NaN')
        return np.nan
    #
    ##############################################
    # Compute the factor
    #    GPP/GPP_st
    #  reshpae is used to maintain the dimensions
    f_gpp_isop=(GPP.reshape(5,-1)/GPP_st[:,None]).reshape(GPP.shape) 
    ########################
    return f_gpp_isop


