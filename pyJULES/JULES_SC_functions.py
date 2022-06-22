#
#
# Python module of JULES soil carbon functions
#
# Edward Comyn-Platt, edwcom@ceh.ac.uk
# Centre for Ecology and Hydrology
# 2015
#
import numpy as np

###############################################################################
# Function: JULES_TOTAL_LITTERFALL
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the total carbon from litterfall
#           based on the JULES formulation in D. Clark et al (2011)
#
#          Calculates the litter fall for each PFT then sums.
#
# Required Inputs:
#  PFTf - fraction of each PFT for each point.
#           list or np.array, dimension = [ nPFTs, nPoints ]
#  Local_Litterfall - The Litterfall for each PFT for each point
#                     list or np.array, dimension = [ nPFTs, nPoints ]
#                     (calculate using JULES_LOCAL_LITTERFALL)
#                      units - kgC m^-2 year^-1
#  Carbon_vege - Vegetation Carbon Density for each PFT for each point
#                list or np.array, dimension = [ nPFTs, nPoints ]
#                (calculate using JULES_LOCAL_LITTERFALL with get_Cv set to True)
#                 units - kgC m^-2
# 
# Optional Inputs:
#  Competition - Competition component of the Litterfall equation (63)
#                 If included must be on the same [ nPFTs,nPoints ] grid as
#                 Local_Litterfall.
#                 Default = None, i.e. no competition
#  gamma_v    - the Disturbance turnover rate for each PFT.
#                 Default = [0.005,0.007,0.20,0.20,0.05] (JULES default)
#  Per_PFT    - logical switch to output total_litterfall for each PFT
#                 Default = False
#
# Output:
#  TOTAL_Litterfall - Total grid box litterfall.
#                        np.array( nPoints ) 
#                     or np.array( nPFTs, nPoints ) if Per_PFT=True
#                      units - kgC m^-2 year^-1
#
###############################################################################
def JULES_TOTAL_LITTERFALL( PFTf, Local_Litterfall, Carbon_vege,  \
                            Competition=None,                     \
                            gamma_v=[0.010,0.004,0.10,0.10,0.05], \
                            Per_PFT=False):
    
    ###########################################################################
    # 1. Convert all data to np.array
    ###################################
    PFTf=np.array(PFTf)
    Local_Litterfall=np.array(Local_Litterfall)
    Carbon_vege=np.array(Carbon_vege)
    gamma_v=np.array(gamma_v)
        
    ###########################################################################
    # 2. check dimensions
    ###########################
    # 2.1 nPFTs is set to number of gamma_v values
    nPFTs = len(gamma_v)
    #
    # 2.2 check PFTf has correct number of PFTs
    if (PFTf.shape[0]!=nPFTs):
        print( 'ERROR: PFTf and gamma_v have inconsistent nPFTs')
        print( 'PFTf.shape = ', PFTf.shape )
        print( 'gamma_v.shape = ',gamma_v.shape)
        print( 'Returning NaN')
        return np.nan
    #
    # 2.3 check PFTf and Local_Litterfall are same shape
    if (Local_Litterfall.shape!=PFTf.shape):
        print( 'ERROR: Local_Litterfall and PFTf have inconsistent dimensions')
        print( 'Local_Litterfall.shape = ', Local_Litterfall.shape )
        print( 'PFTf.shape = ', PFTf.shape)
        print( 'Returning NaN')
        return np.nan
    #
    # 2.4 check Carbon_vege and PFTf are same shape
    if (PFTf.shape!=Carbon_vege.shape):
        print( 'ERROR: PFTf and Carbon_vege have inconsistent dimensions')
        print( 'PFTf.shape = ', PFTf.shape )
        print( 'Carbon_vege.shape = ', Carbon_vege.shape)
        print( 'Returning NaN')
        return np.nan
    #
    # 2.5 if competition exists check dimensions, else create dummy array of zeros
    if (Competition!=None):
        if (PFTf.shape!=Competition.shape):
            print( 'ERROR: PFTf and Competition have inconsistent dimensions')
            print( 'PFTf.shape = ', PFTf.shape )
            print( 'Competition.shape = ', PFTf.shape)
            print( 'Setting Competition Component To Zero And Continuing')
            Competition=np.zeros_like(Local_Litterfall)
        else:
            print( 'Competition term included in calculation')
        #
    else:
        Competition=np.zeros_like(PFTf)
    #  
    # 2.5 If PFTf has only one dimension (i.e. one point),
    #       add additional dimensions for array indexing
    if (len(PFTf.shape)==1):
        PFTf=PFTf.reshape(PFTf.shape[0],1)
        Local_Litterfall=Local_Litterfall.reshape(Local_Litterfall.shape[0],1)
        Carbon_vege=Carbon_vege.reshape(Carbon_vege.shape[0],1)
        Competition=Competition.reshape(Competition.shape[0],1)
    
    ###########################################################################
    # 3. Calculate total Litterfall by summing over all PFTs
    ########################################################
    # 3.1 if Per_PFT not set then sum over all PFTs
    if not (Per_PFT):
        # 3.1.1 Create TOTAL_Litterfall array of dimensions nPoints
        TOTAL_Litterfall = np.zeros( Local_Litterfall.shape[1:] )
        #
        # 3.1.2 Loop over each PFT, summing to TOTAL as we go, (Eqn 63)
        for PFT in range(nPFTs):
            TOTAL_Litterfall +=    PFTf[PFT,:] * \
                                   ( Local_Litterfall[PFT,:]           + \
                                     (gamma_v[PFT]*Carbon_vege[PFT,:]) + \
                                     Competition[PFT,:]                  )
            #
        #
    # 3.2 if Per_PFT set keep on PFT grid
    else:
        # 3.2.1 Create TOTAL_Litterfall array of dimensions [nPFTs,nPoints]
        TOTAL_Litterfall = np.zeros( Local_Litterfall.shape )
        #
        # 3.2 Loop over each PFT, (Eqn 63 but not summed)
        for PFT in range(nPFTs):
            TOTAL_Litterfall[PFT,:] =    PFTf[PFT,:] * \
                                         ( Local_Litterfall[PFT,:]           + \
                                           (gamma_v[PFT]*Carbon_vege[PFT,:]) + \
                                           Competition[PFT,:]                  )
            #
        #
        
    ############################################################################
    # 3. Return output
    ########################
    return TOTAL_Litterfall

###############################################################################
# Function: JULES_LOCAL_LITTERFALL
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the local carbon litterfall 
#           based on the JULES formulation in D. Clark et al (2011)
#           equations 47-59
#
# Required Inputs:
#  LAI - LAI for each PFT for each point.
#          list or np.array, dimensions = [ nPFTs, nPoints ]
#           m^2 / m^2
#  Tleaf - Leaf temperature for each PFT for each point or just for each point.
#          if given for each point only Tleaf is replicated onto each PFT
#          list or np.array, dimensions = [ nPoints ] or [ nPFTs, nPoints ]
#           units - K or C
# 
# Optional Inputs:
#  gamma_0 - Minimum leaf turnover rate.
#             Default = [0.25,0.25,0.25,0.25,0.25] (JULES defaults)
#             units = (360 days)^-1
#  dSM     - Rate of change of turnover with Soil Moisture
#             Default = [0.,0.,0.,0.,0.] (JULES defaults, i.e. currently ignored)
#             units = (360 days)^-1       
#  dT      - Rate of change of turnover with temperature
#             Default = [9.0,9.0,9.0,9.0,9.0] (JULES defaults)
#             units = (360 days K)^-1
#  Toff    - Threshold temperature for leaf mortality change
#             Default = [278.15,233.15,278.15,278.15,233.15] (JULES defaults)
#             units = K
#  p_stat   - phenological status, if provided on the same grid as LAI,
#             it is used to describe the phenological status of the vegetation, 
#             thus improved leaf mortality estimate. eqns 48-50
#             Default = None, 
#             units = none
#  LAI_b   - Balanced (or seasonal maximum) LAI, if provided on the same grid as LAI,
#             it is used to calculate p_stat (see above), 
#             Default = None, 
#             units = m^2 / m^2
#  CanHt   - Canopy height, if included on the same grid as LAI, it is used
#             to calculate LAI_b which in turn is used to improve 
#             leaf mortality estimate. eqns 58 and 61
#             Default = None
#             units = m
#  gamma_p - Rate of leaf growth, required for p_stat calculation.
#             Default = [15.0,20.0,20.0,20.0,20.0] (JULES defaults)
#             units = (360 days)^-1
#  a_wl    - Allometric coefficient,
#             Default = [0.65,0.65,0.005,0.005,0.10] (JULES defaults)
#             units = kgC m^-2
#  a_ws    - Ratio of total to resp stem carb, only required for CanHt to LAI_b calc
#             Default = [10.0,10.0,1.0,1.0,10.0] (JULES defaults)
#             units = none
#  b_wl    - Allometric exponent, 
#             Default = [1.667,1.667,1.667,1.667,1.667] (JULES defaults)
#             units = none
#  nu_sl   - Live stemwood coefficient, only required for CanHt to LAI_b calc
#             Default = [0.01,0.01,0.01,0.01,0.01] (JULES defaults)
#             units = kgC m^-2 per unit LAI
#  sigma_l - Specific Leaf Density
#             Default = [0.0375,0.1,0.025,0.05,0.05] (JULES defaults)
#             units = kgC m^-2 per unit LAI
#  gamma_r - Root turnover rate
#             Default = [0.25,0.15,0.25,0.25,0.25] (JULES_defaults)
#             units = (360 days)^-1
#  gamma_w - Wood/Stem turnover rate
#             Default = [0.005,0.005,0.20,0.20,0.05] (JULES_defaults)
#             units = (360 days)^-1
#  get_Cv  - logical operator for optional output.
#             set to True to return the vegetation carbon density for each PFT
#             if set to true output returned as:
#              Local_Litterfall, Carbon_vege = JULES_LOCAL_LITTERFALL([input])
#             Default = False
# 
# Outputs:
#  Local_Litterfall - Local litter fall rate for each PFT at each point
#                     np.array, dimensions = [ nPFTs, nPoints ]
#                     units - kgC m^-2 year^-1
# 
# Optional Output:
#  Carbon_vege - Vegetation Carbon Density for each PFT ast each point
#                  np.array, dimensions = [ nPFTs, nPoints ]
#                     units - kgC m^-2
#
# Notes: 
#  1. If neither LAI_b or CanHt are provided then LAI is assumed to equal LAI_b
#  2. 
#
 
def JULES_LOCAL_LITTERFALL( LAI, Tleaf,                                \
                            gamma_0=[0.25,0.25,0.25,0.25,0.25],        \
                            dSM=[0.,0.,0.,0.,0.],                      \
                            dT=[9.0,9.0,0.0,0.0,9.0],                  \
                            Toff=[278.15,243.15,258.15,258.15,243.15], \
                            LAI_b=None, CanHt=None, p_stat=None,       \
                            gamma_p=[20.0,20.0,20.0,20.0,20.0],        \
                            a_wl=[0.65,0.65,0.005,0.005,0.10],         \
                            a_ws=[10.0,10.0,1.0,1.0,10.0],             \
                            b_wl=[1.667,1.667,1.667,1.667,1.667],      \
                            nu_sl=[0.01,0.01,0.01,0.01,0.01],          \
                            sigma_l=[0.0375,0.1,0.025,0.05,0.05],      \
                            gamma_r=[0.25,0.25,0.25,0.25,0.25],        \
                            gamma_w=[0.01,0.01,0.20,0.20,0.05],        \
                            get_Cv=False,get_components=False,         \
                            with_phenol=True, get_gamma_lm=False       ):
    
    ####################################################################
    # 1. convert input data to np.arrays
    ######################################
    LAI=np.array(LAI)
    Tleaf=np.array(Tleaf)
    gamma_0=np.array(gamma_0)
    dSM=np.array(dSM)
    dT=np.array(dT)
    Toff=np.array(Toff)
    sigma_l=np.array(sigma_l)
    gamma_r=np.array(gamma_r)
    gamma_w=np.array(gamma_w)
    a_wl=np.array(a_wl)
    b_wl=np.array(b_wl)
    
    ####################################################################
    # 2. Check PFT parameter array dimensions:
    ##########################################
    # 2.1 Essential PFT params (gamm_0, dT and Toff)
    if ((gamma_0.shape!=dT.shape)      or \
        (gamma_0.shape!=Toff.shape)    or \
        (gamma_0.shape!=sigma_l.shape) or \
        (gamma_0.shape!=a_wl.shape)    or \
        (gamma_0.shape!=b_wl.shape)    or \
        (gamma_0.shape!=gamma_r.shape) or \
        (gamma_0.shape!=gamma_w.shape)    ):
        print( 'ERROR in JULES_LOCAL_LITTERFALL')
        print( 'Essential PFT parameters not consistent in shape')
        print( 'gamma_0.shape = ', gamma_0.shape)
        print( 'dT.shape = ',dT.shape)
        print( 'Toff.shape = ', Toff.shape)
        print( 'sigma_l.shape = ', sigma_l.shape)
        print( 'a_wl.shape = ', a_wl.shape)
        print( 'b_wl.shape = ', b_wl.shape)
        print( 'gamma_r.shape = ', gamma_r.shape)
        print( 'gamma_w.shape = ', gamma_w.shape)
        print( 'Returning NaN, sort it out mate!')
        return np.nan
    # Set nPFTs to length of gamma_0
    nPFTs = len(gamma_0)
    #
    # 2.2 Check first LAI dimension is equal to nPFTs
    if (LAI.shape[0]!=nPFTs):
        print( 'ERROR in JULES_LOCAL_LITTERFALL')
        print( 'LAI has different number of PFTs than PFT parameters!')
        print( 'Returning NaN, sort it out mate!')
        return np.nan
    #
    # 2.3 Check Tleaf dims
    # if equal to LAI do nothing
    if (LAI.shape[:]!=Tleaf.shape):
        if (LAI.shape[1:]==Tleaf.shape):
            print( 'Replicating Tleaf onto LAI grid')
            Tleaf = np.array([ Tleaf for PFT in range(nPFTs)])
        else:
            print( 'ERROR in JULES_LOCAL_LITTERFALL')
            print( 'LAI and and Tleaf have different spatial dimensions!')
            print( 'Returning NaN, sort it out mate!')
            return np.nan
    else:
        print( 'LAI and Tleaf on the same grid.')
    #
    # 2.4 if p_stat provided check dimensions
    if (p_stat!=None):
        if (p_stat.shape!=LAI_shape):
            print( 'WARNING in JULES_LOCAL_LITTERFALL')
            print( 'p_stat does not match LAI')
            print( 'setting p_stat to None')
            p_stat=None
    #
    # 2.5 if LAI_b provided check dimensions
    if (LAI_b!=None):
        if (LAI_b.shape!=LAI.shape):
            print( 'WARNING in JULES_LOCAL_LITTERFALL')
            print( 'LAI_b does not match LAI')
            print( 'setting LAI_b to None')
            LAI_b=None
    #
    # 2.6 if CanHt provided check dimension
    if (CanHt!=None):
        if (CanHt_b.shape!=CanHt.shape):
            print( 'WARNING in JULES_LOCAL_LITTERFALL')
            print( 'CanHt does not match LAI')
            print( 'setting CanHt to None')
            CanHt=None
    #
    # 2.7 if LAI_b or CanHt provided set additional params to np.arrays and check dims
    # 2.7.1 gamma_p required for p_stat, LAI_b and CanHt
    if (p_stat!=None) or (LAI_b!=None) or (CanHt!=None):
        gamma_p=np.array(gamma_p)
        if (gamma_p.shape!=gamma_0.shape):
            print( 'ERROR in JULES_LOCAL_LITTERFALL')
            print( 'gamma_p not in consistent in shape')
            print( 'gamma_0.shape = ', gamma_0.shape)
            print( 'gamma_p.shape = ', gamma_p.shape)
            print( 'setting LAI_b and CanHt to None')
            LAI_b=None
            CanHt=None
    #
    # 2.7.2 a_wl, a_ws, b_wl and nu_sl required for CanHt only
    if (CanHt!=None):
        a_ws=np.array(a_ws)
        nu_sl=np.array(nu_sl)
        if ((gamma_0.shape!=a_ws.shape)   or \
            (gamma_0.shape!=nu_sl.shape)    ):
            print( 'ERROR in JULES_LOCAL_LITTERFALL')
            print( 'CanHt PFT parameters not consistent in shape')
            print( 'gamma_0.shape = ', gamma_0.shape)
            print( 'a_ws.shape = ', a_ws.shape)
            print( 'nu_sl.shape = ', nu_sl.shape)
            print( 'setting CanHt to None')
            CanHt=None
        #
    #
    # 2.8 After all check's have been done check for surplus data
    # 2.8.1 If both LAI_b and CanHt provided choose LAI_b 
    if (LAI_b!=None) and (CanHt!=None):
        print( 'WARNING in JULES_LOCAL_LITTERFALL')
        print( 'Both LAI_b and CanHt have been provided')
        print( 'Setting CanHt to None to reduce unnecessary computation')
        CanHt=None
    # 2.8.1 If both p_stat and LAI_b provided choose p_stat
    if p_stat and LAI_b:
        print( 'WARNING in JULES_LOCAL_LITTERFALL')
        print( 'Both p_stat and LAI_b have been provided')
        print( 'Setting LAI_b to None to reduce unnecessary computation')
        LAI_b=None
    
    #####################################################################
    # 3. Calculate leaf mortality rates (eqns 47-50, 58 and 60)
    #################################
    #
    # 3.1 Calculate leaf mortality rate without phenology
    gamma_lm = np.zeros_like(LAI)
    for PFT in range(nPFTs):
        #print( 'PFT - ', PFT)
        #if PFT==0:
        #    for TL in Tleaf[PFT,:].flat:
        #        print( TL)
        # Set all values to min leaf mortality
        gamma_lm[PFT,:]=gamma_0[PFT]
        #print( gamma_lm)
        # where Tleaf is less than Toff apply eqn
        #print( Tleaf[PFT,:].shape)
        #print( gamma_lm[PFT,:].shape)
        #print( np.min(gamma_lm[PFT,:]))
        tempindex=np.where(Tleaf[PFT,:]<=Toff[PFT])
        if (len(tempindex[0])>0):
            gamma_lm[PFT,:][tempindex]       \
                =  gamma_0[PFT] * \
                   ( 1 + (dT[PFT] * ( Toff[PFT] - \
                                      Tleaf[PFT,:][tempindex])) )
            #print( gamma_0[PFT])
            #print( dT[PFT])
            #print( np.max( ( Toff[PFT] - \
            #                         Tleaf[PFT,:][tempindex])) )
            #print( np.min( ( Toff[PFT] - \
            #                          Tleaf[PFT,:][tempindex])) )
        #print( np.min(gamma_lm[PFT,:]) )
    #
    
    # 3.2 If CanHt provided calculate LAI_b using inverse of eqns 58 and 60
    if (CanHt!=None):
        LAI_b=np.zeros_like(CanHt)
        for PFT in range(nPFTs):
            LAI_b[PFT,:]= ( a_ws[PFT] * eta_sl[PFT] * CanHt[PFT,:] / a_wl[PFT] ) \
                          **   ( 1. / (b_wl[PFT]-1. )  )                   
    #
    # 3.3 If LAI_b provided or calculated from CanHt, calculate phenological status, p_stat
    if (LAI_b!=None):
        p_stat = LAI/LAI_b
    # elif p_stat provided then calculate LAI_b from p_stat and LAI
    elif (p_stat!=None):
        LAI_b  = LAI/p_stat
    # else, i.e. nothing has been provided, set LAI_b to LAI and p_stat to 1s
    else:
        LAI_b  = LAI
        p_stat = np.ones_like(LAI)
    #
    # 3.4 Calculate effective leaf turnover rate (eqn 50 (and 49) )
    if with_phenol==True:
        gamma_l = np.zeros_like(gamma_lm)
        for PFT in range(nPFTs):
            # set all vals to max val
            gamma_l[PFT,:]=gamma_p[PFT]
            # set any values where gamma_lm is less than 2*gamma_0 to p*gamma_lm
            index = np.where( gamma_lm[PFT,:]<=(2*gamma_0[PFT]) )
            if (len(index[0])>0):
                gamma_l[PFT,:][index]  \
                    =  p_stat[PFT,:][index] * gamma_lm[PFT,:][index]
    else:
        gamma_l = gamma_lm
    
      
    ###########################################################################
    # 4. Calculate the changes in Carbon Density (Eqns 56-58) 
    ################################
    # 4.1 Leaf Carbon Change (eqn 56):
    Leaf_C = np.array([ sigma_l[PFT]*LAI_b[PFT,:] for PFT in range(nPFTs) ])
    # 4.2 Root Carbon Change (=Leaf CC, eqn 57):
    Root_C = Leaf_C
    # 4.3 Stem Carbon Change (eqn 58):
    Stem_C = np.array([a_wl[PFT] * (LAI_b[PFT,:] ** b_wl[PFT]) \
                          for PFT in range(nPFTs)] )
    
    ###########################################################################
    #  5. Calculate Local_Litterfall (eqn 59)
    ####################################
    # 5.1, Leaf Component
    Leaf_LF = gamma_l*Leaf_C
    # 5.1 Root Component
    Root_LF = np.array([gamma_r[PFT]*Root_C[PFT,:] for PFT in range(nPFTs)])
    # 5.2 Stem Component
    Stem_LF = np.array([gamma_w[PFT]*Stem_C[PFT,:] for PFT in range(nPFTs)])
    #
    # 5.3 Sum components
    Local_Litterfall = Leaf_LF + Root_LF + Stem_LF
    
    ###########################################################################
    #  6. If get_Cv set to true calculate vegetation carbon (eqn 55)
    ####################################
    if get_Cv:
        Carbon_vege = Leaf_C+Root_C+Stem_C
    
    ###########################################################################
    # 6. Return Output
    if get_Cv:
        return Local_Litterfall, Carbon_vege
    elif get_components:
        return Leaf_LF, Root_LF, Stem_LF
    elif get_gamma_lm:
        return Local_Litterfall, gamma_l
    else:
        return Local_Litterfall
    #
#

###############################################################################
# Function: JULES_LITTER_to_SCpool
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the fraction of Litterfall carbon that
#          That goes into each Cpool (RPM and DPM, nPools=2)
# 
# Require Inputs:
#  Litterfall - Litterfall rate for each PFT
#                Dimensions - [ nPFTs, nPoints ]
#                 units - kgC m^-2 year^-1
#
# Optional Inputs:
#  alpha_dr - factor controling the ratios
#              Dimensions - [ nPFTs ]
#              Default = [0.25,0.25,0.67,0.67,0.33] (JULES Default)
#
# Outputs:
#  Soil_Carbon - Soil Carbon in each pool from input litter fall 
#                 Dimension = [ nPools, nPoints ]
#                  units - kgC m^-2 year^-1 
#
###############################################################################
def JULES_LITTER_to_SCpool(Litterfall, \
                           alpha_dr=[0.25,0.25,0.67,0.67,0.33] ):
    #
    #########################################################################
    # 1. Convert inputs to np.arrays
    ################################
    Litterfall=np.array(Litterfall)
    alpha_dr=np.array(alpha_dr)
    nPools=2
    
    #########################################################################
    # 2. Check Dimensions
    ################################
    if (Litterfall.shape[0]!=alpha_dr.shape[0]):
        print( 'ERROR in JULES_LITTER_to_SCpool' )
        print( 'First dimension of Litterfall should correspond to length of alpha_dr!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    #
    if (len(Litterfall.shape)==(1)):
        print( 'Litterfall has only one dimension, adding pseudo dimension for nPoints' )
        Litterfall=Litterfall.reshape(Litterfall.shape[0],1)
    #
    # Create Soil Carbon Dimensions, [nPools,nPoints]
    Litterfall_dims=np.array(Litterfall.shape)
    nPFTs   = Litterfall_dims[0]
    SC_dims = Litterfall_dims
    SC_dims[0]=nPools
    
    #########################################################################
    # 3. Calculate Weighting Fractions
    ###################################
    f_dpm = alpha_dr / (1.+alpha_dr)
    #
    #print( f_dpm )
    
    #########################################################################
    # 4. Apply to Litterfall components
    ####################################
    # 4.1 Create Soil_Carbon array of dimensions [ nPools, nPoints ]
    Soil_Carbon = np.zeros(SC_dims,dtype='float64')
    #
    # 4.2 Sum all the Litterfall components for each PFT going into each pool
    for PFT in range(nPFTs):
        # DPM fraction:
        Soil_Carbon[0,:] += f_dpm[PFT]*Litterfall[PFT,:]
        Soil_Carbon[1,:] += (1-f_dpm[PFT])*Litterfall[PFT,:]
        #print( np.sum(Soil_Carbon[0,:],dtype='float64'), np.sum(Soil_Carbon[1,:],dtype='float64') )
        #print( np.min(Soil_Carbon[0,:]), np.max(Soil_Carbon[1,:]) )
        #print( np.max(Soil_Carbon[0,:]), np.max(Soil_Carbon[1,:]) )
    #
    
    #########################################################################
    # 4. Return Soil_Carbon
    ####################################
    return Soil_Carbon

###############################################################################
# Function: JULES_SOIL_RESPIRATION
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to  calculate the soil respiration given the 
#         soil carbon, soil temperature, soil moisture and vegetation fraction.
#         based on the JULES formulation in D. Clark et al (2011)
#
# Required Inputs: 
#  Cs - Soil Carbon - list, 'list of lists' or np.array 
#                     dimensions [ nPools , nPoints ]
#                      units - kgC m^-2
#  Ts - Soil Temperature - list or np.array 
#                     dimensions [ nPoints ]
#                      units - K or C
#  SM - Soil Moisture  - list or np.array 
#                     dimensions [ nPoints ]
#                      units - fraction of saturation
#  SMw - Soil Moisture at wilting point - list or np.array 
#                     dimensions [ nPoints ]
#                      units - fraction of saturation
#  Vf - Vegetation Fractional Cover - list or np.array
#                     dimensions [ nPoints ]
#
# Optional Inputs:
# kappa - Soil Specific Respiration Rate
#          kappa is used to define nPool, i.e. nPools=len(kappa)
#          list or np.array - Dimensions [ nPools ]
#          Default       = [ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-8 ]
#          Correponding to [ DPM,     RPM,     BIO,     HUM ]
#          See JULES documentation for more info.
# Tfunc - Temperature Function to use.
#          String, Default='Q10'
#          Allowed values: 'Q10', 'RothC_T'
# Q10s  - Q_10_soil factor, only relevent if Tfunc='Q10'.
#          number, Default=2.0
# SMfunc - Soil Moisture function to use
#           String, Default='JULES_std'
#           Allowed values: 'JULES_std'
# 
# OUTPUT_opt - output option:
#                'pools' - (Default) returns respiration for each pool (Rs_pools)
#                'total' - returns total soil respiration (Rs)
#                'both'  - returns both of above in order Rs, Rs_pools
#           (see below)
#
# Outputs:
# Rs_pools - Soil Respiration for each pool, 
#              np.array dimensions [ nPools, nPoints ]
#               units - kgC m^-2 s^-1
#     AND/OR
# Rs       - Total soil respiration
#              np.array, dimensions [ nPoints ]
#               units - kgC m^-2 s^-1
#
# Required functions:
#    JULES_Q10_func
#    JULES_RothC_T_func
#    JULES_std_Vf_fun     
#    JULES_std_Vf_fun
#
# Notes:
#   - nPoints dimension can be multi-dimensional, but must be consistent
#        accross all variables
#
#############################################################################
def JULES_SOIL_RESPIRATION( Cs, Ts, SM, SMw, Vf, \
                            kappa=[3.22e-7,9.65e-9,2.12e-8,6.43e-10], \
                            Tfunc='Q10', Q10=2.0, \
                            SMfunc='JULES_std', \
                            Vffunc='JULES_std', \
                            OUTPUT_opt='pools'  ):
    
    #####################################################################
    # 1. Check parameters
    ######################
    # 1.1 kappa
    # Convert to np.array
    kappa=np.array(kappa)
    # Check is 1D
    if len(kappa.shape)!=1:
        print( 'kappa should have only one dimension!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    # nPools equal to number of kappa values
    nPools=len(kappa)
    #
    #######################
    # 1.2 Cs
    # convert to np.array if neccessary
    Cs=np.array(Cs)
    # check first dim is equal to nPools
    if Cs.shape[0]!=nPools:
        print( 'First dimension of Cs should correspond to nPools!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    #
    #######################
    # 1.3 Ts
    # convert to np.array if neccessary
    Ts=np.array(Ts)
    # check dimensions are equal to Cs (exclusing nPools dim)
    if Ts.shape!=Cs.shape[1:]:
        print( 'Ts dimensions do not match Cs!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    # Check for Kelvin
    if np.max(Ts)<=100.:
        print( 'Max Ts is less than 100, I doubt this is Kelvin' )
        print( 'Converting Ts to Kelvin seen as you were too lazy!' )
        Ts=Ts+273.15
    #
    #######################
    # 1.4 SM
    # convert to np.array if neccessary
    SM=np.array(SM)
    # check dimensions are equal to Cs (exclusing nPools dim)
    if SM.shape!=Cs.shape[1:]:
        print( 'SM dimensions do not match Cs!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    #
    #######################
    # 1.5 Vf
    # convert to np.array if neccessary
    Vf=np.array(Vf)
    # check dimensions are equal to Cs (exclusing nPools dim)
    if Vf.shape!=Cs.shape[1:]:
        print( 'Vf dimensions do not match Cs!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan
    
    #####################################################################
    # 2. Apply relevent T function
    ################################
    if Tfunc=='Q10':
        F_Ts = JULES_Q10_T_func(Ts,Q=Q10)
    elif Tfunc=='RothC':
        F_Ts = JULES_RothC_T_func(Ts)

    
    #####################################################################
    # 3. Apply relevent SM function
    ################################
    if SMfunc=='JULES_std':
        F_SM = JULES_std_SM_func(SM,SMw)
    
    #####################################################################
    # 4. Apply relevent Vf function
    ###############################
    if Vffunc=='JULES_std':
        F_Vf = JULES_std_Vf_func(Vf)
    
    #####################################################################
    # 
    # 5. Calculate the Respiration for each Carbon Pool
    ################################################### 
    # Create Rs array of same dimensions as Cs
    Rs_pools = np.zeros_like( Cs )
    #
    # loop round each soil pool
    for pool in range(len(kappa)):
        Rs_pools[pool,:] =   kappa[pool] \
                           * Cs[pool,:]  \
                           * F_Ts[:]     \
                           * F_SM[:]     \
                           * F_Vf[:]
    
    ####################################################################
    # 6. Calculate the total soil respiration if required
    #     (sum over soil pools)
    #########################################
    if (OUTPUT_opt=='total') or (OUTPUT_opt=='both'):
        Rs = np.sum(Rs_pools,axis=0)
    
    ####################################################################
    # 7. Output requested data
    ##########################################
    if (OUTPUT_opt=='pools'):
        return Rs_pools
    elif (OUTPUT_opt=='total'):
        return Rs
    elif (OUTPUT_opt=='both'):
        return Rs, Rs_pools
    else:
        print( 'Unrecognised OUTPUT_opt, returning: Rs, Rs_pools' )
        print( "P.S. you're a numpty" )
        return Rs, Rs_pools

###############################################################################
# Function: JULES_Q10_T_func
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the JULES q10 function
# Required Inputs: 
#  T - Temperature (K or C, calcualtion done in K) 
#            - single val, list or np.array 
#
# Optional Inputs:
# Q    - Q factor, Default=2.0
# Tref - Reference temperature, Default=298.15 (K)
#
# Outputs:
# F_T - Q_10 function, np.array 
##############################################################################
def JULES_Q10_T_func( T, \
                      Q=2.0, Tref=282.4 ):
    #
    # convert T to np.array
    T=np.array(T)
    #
    # convert to K if max T less than 100
    if (np.max(T)<=100.):
        print( 'Max T is less than 100, I doubt this is Kelvin' )
        print( 'Converting T to Kelvin seen as you were too lazy!' )
        T=T+273.15
    #
    F_T =  Q ** ( (T-Tref)/10. )
    #
    #print( F_T )
    return F_T

###############################################################################
# Function: JULES_RothC_T_func
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the RothC Temperature function
# Required Inputs: 
#  T - Temperature (K or C, calcualtion done in K) 
#            - single val, list or np.array 
#
# Outputs:
# F_T - RothC Temperature function, np.array 
##############################################################################
def JULES_RothC_T_func( T, minT=260. ):
    #
    # convert T to np.array
    T=np.array(T,dtype='float64')
    #
    # convert to K if max T less than 100
    if (np.max(T)<=100.):
        print( 'Max T is less than 100, I doubt this is Kelvin' )
        print( 'Converting T to Kelvin seen as you were too lazy!' )
        T=T+273.15
    #
    # Any points with T less than 254.85 to be set to 254.85 to remove negative exponents
    T_wk = T
    if len(np.where(T_wk<=minT)[0])>0:
        T_wk[np.where(T_wk<=minT)]=minT    # plus a teenytiny amount to remove infinites
    #
    #print( np.min(T_wk) )
    #
    exp_part = np.exp(  106./ (T_wk-254.85) )
    #
    F_T = 47.9 / (1 + exp_part)
    #print( F_T#[np.where(T_wk==minT)[0]] )
    #
    return F_T 
    #

###############################################################################
# Function: JULES_std_SM_func
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the standard JULES SM respiration function
# Required Inputs: 
#  SM - Soil Moisture as a fraction of saturation
#         list or np.array
#  SMw - Soil Moisture at wilting point as a fraction of saturation
#         list or np.array
#         Must be same dimensions as SM
# Optional Inputs:
#  SM_min_factor - Factor for SM min calculation
#                   Default = 1.0
#  SM_cap - Cap for soil moisture, set to value to cap SM at, 
#             as a fraction of saturation, 
#             e.g. SM_cap=1.0 as SM can exceed 1.0 in JULES             
#            
# Outputs:
#  F_SM - Standard JULES SM function, np.array 
#
##############################################################################
def JULES_std_SM_func( SM, SMw, \
                       SM_min_factor=1.7, SM_cap=False):
    #
    # convert SM and SMw to np.arrays
    SM  = np.array(SM)
    SMw = np.array(SMw)
    #
    # check SM and SMw are same shape
    if not (SM.shape==SMw.shape):
        print( 'SM and SMW are not same shape!' )
        print( 'Returning NaN, sort it out mate!' )
        return np.nan()
        #
    # calculate optimum and minimum soil moisture
    SM_opt = 0.5 * (1+SMw)
    SM_min = SM_min_factor * SMw
    #
    # Apply SM_cap if present
    if SM_cap:
        SM[np.where(SM>SM_cap)]=SM_cap
    #
    # Apply appropriate SM function to data
    # First create F_SM array same shape as SM with value for SM < SM_min (0.2)
    F_SM=np.zeros_like(SM)+0.2
    #
    # For SM vals between SM_min and SM_opt:
    index =  np.where((SM>SM_min)&(SM<=SM_opt)) 
    F_SM[ index ] = \
         0.2 + ( 0.8 * ((SM[index]-SM_min[index])/(SM_opt[index]-SM_min[index])) )
    del index
    #
    # For SM vals greater than SM_opt
    index = np.where(SM>SM_opt)
    F_SM[ index ] = \
         1. - ( 0.8 * (SM[index]-SM_opt[index]) )
    del index
    #
    return F_SM

###############################################################################
# Function: JULES_std_Vf_func
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Function to calculate the standard JULES Vf respiration function
#
# Required Inputs: 
#  Vf - Vegetation fraction, number list or np.array
#            
# Outputs:
#  F_Vf - Standard JULES Vf function, np.array 
#
##############################################################################
def JULES_std_Vf_func( Vf ):
    #
    # convert Vf to np.array
    Vf=np.array(Vf)
    #
    F_Vf = 0.6+ (0.4*(1.-Vf))
    #
    return F_Vf
#


###############################################################################
# Function: JULES_Cs_Exponential_Decay
# Author: Edward Comyn-Platt, Jan 2016
# Purpose: Function to calculate the exponential decay of CS 
#
# Required Inputs: 
#  t - time series of the decay, units must correspond to kappa units
#                                units = 1/kappa.units
#  C0 - Initial Soil Carbon
#  C_ss - Steady State soil carbon
#  alpha = a*kappa, where kapp has units 1/t.units
#            a is a modifier for the resp. rate to account for litterfall top-up,
#              this should be negative
#            
# Outputs:
#  Cs_t - Soil Carbon as a function of t
#
##############################################################################
def JULES_Cs_Exponential_Decay(t,C0,C_ss,alpha):
    #
    #  Calculate the exponential decay
    Cs_t = (C0-C_ss)*np.exp(alpha*t)
    #
    # Add on the steady state value
    Cs_t += C_ss
    #
    return Cs_t
#






