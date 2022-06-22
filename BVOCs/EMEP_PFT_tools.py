#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: To output an LAI seasonality for EMEP4UK based on the EMEP method
#           see Simpson et al 2012, The EMEP MSC-W chemical transport model
# 
##################################################################################
import numpy as np
import netCDF4 as nc
import sys

#########################################################################################
# Class - CH_functions
# Purpose  - Series of functions for calculating the latitudinal 
#            dependent Canopy Height for the various EMEP PFTs
#            The calculations are based on that described in the EMEP model documentation
#            Simpson et. al, 2012 - The EMEP MSC-W chemical transport model
#########################################################################################
class CH_functions:
    ####################################################################################
    # CF canopyt Height is 20m up to 60 degrees north 
    # whereupon it begins to reduce by 5% per latitude
    # down to a minimum of 6m at 74 degrees north.
    ###########################################################################
    # Input: latitudes - 1D or 2D numpy array of latitudes
    # Output: CHeight - array of canopy height vals same shape as latitudes
    ####################################################################################
    def CF(self, latitudes):
        lats_above_60 = latitudes-60.
        lats_above_60[lats_above_60<0.]=0.
        lats_above_60[lats_above_60>14.]=14.
        CHeight= 20. - ( 20*0.05*lats_above_60)
        return CHeight
    ####################################################################################
    # DF canopyt Height is 20m up to 60 degrees north 
    # whereupon it begins to reduce by 5% per latitude
    # down to a minimum of 6m at 74 degrees north.
    ###########################################################################
    # Input: latitudes - 1D or 2D numpy array of latitudes
    # Output: CHeight - array of canopy height vals same shape as latitudes
    ####################################################################################
    def DF(self,latitudes):
        lats_above_60 = latitudes-60.
        lats_above_60[lats_above_60<0.]=0.
        lats_above_60[lats_above_60>14.]=14.
        CHeight= 20. - ( 20*0.05*lats_above_60)
        return CHeight





#########################################################################################
# Class - LAI_functions
# Purpose  - Series of functions for calculating the seasonal and or latitudinal 
#            dependent LAI for the various EMEP PFTs
#            The calculations are based on that described in the EMEP model documentation
#            Simpson et. al, 2012 - The EMEP MSC-W chemical transport model
#            
#            All PFT functions take same inputs, dayofyear and lats array
#             even if they are not required in calculation
#########################################################################################
class LAI_functions:
    ####################################################################################
    # Function - SEAS_LAT_LAI_CALC
    # Purpose  - Calculate array of LAI based on day of year and latitude
    #             See EMEP documention for rules
    ###########################################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        SGS/EGE   - 1D or 2D arrays of start/end day of growth periods.
    #                     These should correspond to the latitude grid of the data
    #        LAI_min/max - scalar values of minimum values of LAI for PFT
    #        L_s,L_e   - scalar values of length of growth/decline periods
    #####################################################################################
    def SEAS_LAT_LAI_CALC(self,dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e):
        #
        out_dims=[dayofyear.shape[0]]
        for dim in SGS.shape:
            out_dims.append(dim)
        #
        LAI_range=LAI_max-LAI_min
        # create meshed arrays of day of years and flattened latitudes
        sgs_mesh, doy_mesh = np.meshgrid(SGS.flat,dayofyear)
        egs_mesh = np.meshgrid(EGS.flat,dayofyear)[0]
        # Create LAI array with dimensions equal to meshed arrays, and set all values to LAI_min
        LAI_array = np.zeros_like(doy_mesh,dtype='float32')+LAI_min
        #
        # growth period LAI = LAI_min+ growth_rate*growth_days
        LAI_array[(doy_mesh>sgs_mesh)&(doy_mesh<=sgs_mesh+L_s)] += \
                    (LAI_range/L_s) * ((doy_mesh-sgs_mesh)[(doy_mesh>sgs_mesh)&(doy_mesh<=sgs_mesh+L_s)])
        #
        # max period LAI = LAI_max, also set decline period to LAI max for simple subtraction in next step
        LAI_array[(doy_mesh>sgs_mesh+L_s)&(doy_mesh<=egs_mesh)] = LAI_max
        #
        # Decline period LAI=LAI_max- (decline_rate)*(decline_days)
        LAI_array[(doy_mesh>=egs_mesh-L_e)&(doy_mesh<=egs_mesh)] -= \
                    (LAI_range/L_e) * (doy_mesh - (egs_mesh-L_e) )[(doy_mesh>=egs_mesh-L_e)&(doy_mesh<=egs_mesh)]
        #
        return LAI_array.reshape(out_dims)
    #
    ##################################################################################
    # Function - self.CF
    # Purpose  - Calculate the LAI for coniferous forest,
    #    
    #     LAI is 5 m2/m2 up to 60 degrees north 
    #     whereupon it begins to reduce by 5% per latitude
    #     down to a minimum LAI at 74 degrees north.
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def CF(self, dayofyear, latitudes):
        lats_above_60 = latitudes-60.
        lats_above_60[lats_above_60<0.]=0.
        lats_above_60[lats_above_60>14.]=14.
        LAI = 5. - (5.*0.05*lats_above_60)
        # add time dimension to LAI array
        LAI= np.array([LAI for i in dayofyear])
        return LAI
    #
    ##################################################################################
    # Function - self.DF
    # Purpose  - Calculate the LAI for deciduous forest
    #    
    #     DF LAI is seasonally and latitudinally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def DF(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=0.
        LAI_max=4.
        # Length of growth/decline periods
        L_s = 20.
        L_e = 30.
        # calculate latitudinally dependent start and end of growing season
        lats_above_50=latitudes-50.
        SGS=100.+(lats_above_50*1.5)
        EGS=307.+(lats_above_50*-2.)
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        # Reduce by 5% for every latitude above 60 north up to 74 degree
        lats_above_60 = latitudes-60.
        lats_above_60[lats_above_60<0.]=0.
        lats_above_60[lats_above_60>14.]=14.
        LAI = LAI - (np.array([lats_above_60 for i in dayofyear]) * LAI * 0.05 )
        return LAI
    #
    ##################################################################################
    # Function - self.NF
    # Purpose  - Calculate the LAI for Med. Needleleaf forest
    #    
    #     NF LAI is constant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def NF(self, dayofyear, latitudes):
        LAI = np.array([np.zeros_like(latitudes)+4.0 for i in dayofyear])
        return LAI
    #
    ##################################################################################
    # Function - self.BF
    # Purpose  - Calculate the LAI for Med. Broadleaf forest
    #    
    #     BF LAI is constant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def BF(self, dayofyear, latitudes):
        LAI = np.array([np.zeros_like(latitudes)+4.0 for i in dayofyear])
        return LAI
    #
    ##################################################################################
    # Function - self.TC
    # Purpose  - Calculate the LAI for Temp/Bor crops
    #    
    #     TC LAI is seasonally and latitudinally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def TC(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=0.
        LAI_max=3.5
        # Length of growth/decline periods
        L_s = 70.
        L_e = 20.
        # calculate latitudinally dependent start and end of growing season
        lats_above_50=latitudes-50.
        SGS=123.+(lats_above_50*2.57)
        EGS=213.+(lats_above_50*2.57)
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    #
    ##################################################################################
    # Function - self.MC
    # Purpose  - Calculate the LAI for Med. crops
    #    
    #     MC LAI is seasonally and latitudinally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def MC(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=0.
        LAI_max=3.
        # Length of growth/decline periods
        L_s = 70.
        L_e = 44.
        # calculate latitudinally dependent start and end of growing season
        lats_above_50=latitudes-50.
        SGS=123.+(lats_above_50*2.57)
        EGS=237.+(lats_above_50*2.57)
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    #
    ########################################################################################
    # Function - self.RC
    # Purpose  - Calculate the LAI for Root crops
    #    
    #     RC LAI is seasonally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ########################################################################################
    def RC(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=0.
        LAI_max=4.2
        # Length of growth/decline periods
        L_s = 35.
        L_e = 65.
        # calculate latitudinally dependent start and end of growing season
        SGS=np.zeros_like(latitudes)+130.
        EGS=np.zeros_like(latitudes)+250.
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    #
    ########################################################################################
    # Function - self.SNL
    # Purpose  - Calculate the LAI for Moorland
    #    
    #     ML LAI is seasonally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ########################################################################################
    def SNL(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=2.
        LAI_max=3.
        # Length of growth/decline periods
        L_s = 192.
        L_e = 96.
        # calculate latitudinally dependent start and end of growing season
        SGS=np.zeros_like(latitudes)+0.
        EGS=np.zeros_like(latitudes)+366.
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    ##################################################################################
    # Function - self.GR
    # Purpose  - Calculate the LAI for Grass
    #    
    #     GR LAI is seasonally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #
    # Output: LAI - 1D array of LAI corresponding to day numbers
    ##################################################################################
    def GR(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=2.
        LAI_max=3.5
        # Length of growth/decline periods
        L_s = 140.
        L_e = 135.
        # calculate latitudinally dependent start and end of growing season
        SGS=np.zeros_like(latitudes)+0.
        EGS=np.zeros_like(latitudes)+366.
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    #
    ##################################################################################
    # Function - self.MS
    # Purpose  - Calculate the LAI for Med. Scrub
    #    
    #     MS LAI is constant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def MS(self, dayofyear, latitudes):
        LAI = np.array([np.zeros_like(latitudes)+2.5 for i in dayofyear])
        return LAI
    #
    ##################################################################################
    # Function - self.IAM_CR
    # Purpose  - Calculate the LAI for Generic Crop
    #    
    #     IAM_CR LAI is constant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def IAM_CR(self, dayofyear, latitudes):
        LAI = np.array([np.zeros_like(latitudes)+3.5 for i in dayofyear])
        return LAI
    #
    ##################################################################################
    # Function - self.IAM_DF
    # Purpose  - Calculate the LAI for Generic Decid. Forest
    #    
    #     IAM_DF LAI is seasonally and latitudinally dependant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def IAM_DF(self, dayofyear, latitudes):
        # 
        # LAI min and max values
        LAI_min=0.
        LAI_max=4.0
        # Length of growth/decline periods
        L_s = 15.
        L_e = 30.
        # calculate latitudinally dependent start and end of growing season
        lats_above_50=latitudes-50.
        SGS=105.+(lats_above_50*1.5)
        EGS=297.+(lats_above_50*-2.0)
        #
        LAI = self.SEAS_LAT_LAI_CALC(dayofyear,SGS,EGS,LAI_min,LAI_max,L_s,L_e)
        #
        return LAI
    #
    ##################################################################################
    # Function - self.IAM_MF
    # Purpose  - Calculate the LAI for Generic MF (Med. Forest???)
    #    
    #     IAM_MF LAI is constant
    #     See EMEP documentation for more details
    #    
    ######################################################
    # Input: dayofyear - 1D array of day of year numbers
    #        latitudes - 1D or 2D numpy array of latitudes
    #
    # Output: LAI - 2D or 3d array of LAI corresponding to day numbers and lat dimensions
    ##################################################################################
    def IAM_MF(self, dayofyear, latitudes):
        LAI = np.array([np.zeros_like(latitudes)+5.0 for i in dayofyear])
        return LAI
    #





