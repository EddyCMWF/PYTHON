#!/usr/bin/python2.7
##############################################################
#
# Python module for Flux Site Tools
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
##############################################################
import numpy as np
import scipy.ndimage as im
import maths_tools.MathsTools as MT
import maths_tools.SmoothTools as ST
from scipy.optimize import curve_fit


##############################################################
# Function: calculate_GPP_from_NEE_and_SW
# 
# Purpose: Remove respiration component from NEE,
#          Assume that negative NEE when no sunlight (SW=0)
#          is the respiration component.
#          This is best used for a single day or short time period 
#          as all night points are used in the dark_respiration calculation
# 
# Inputs: NEE - Net Ecosystem Exchange, numpy array
#         SW  - Short Wave Down, numpy array (same dims as NEE)
#                used as a day night flag
# 
# Outputs: GPP - Gross Primary Productivity, numpy array
#                                            (same dims as NEE)
#
##############################################################
def calculate_GPP_from_NEE_and_SW( NEE, SW ):
    
    if (NEE.shape != SW.shape):
        print('ERROR in "calculate_GPP_from_NEE_and_SW":')
        print('NEE and SW have incompatible dimensions')
        print('returning NaN')
        return NaN
    
    index = np.where(SW<=0)[0]
    
    if (len(index)>=1):
        dark_respiration=np.mean(NEE[index])
    else:
        print('WARNING in "calculate_GPP_from_NEE_and_SW":')
        print('No night time data, using min NEE as dark respiration')
        dark_respiration=np.min(NEE)
    
    GPP=NEE-dark_respiration
    
    return GPP

#
##############################################################
# Function: calculate_GPP_from_NEE_SW_T
# 
# Purpose: Remove respiration component from NEE,
#           Assume that negative NEE when no sunlight (SW=0)
#           is the respiration component.
# 
# Inputs: NEE - Net Ecosystem Exchange, numpy array 
#                Dimensions = [ nDays, time_steps_in_day ]
#         SW  - Short Wave Down, numpy array (same dims as NEE)
#                used as a day night flag
#         T   - Temperature - preferably soil temperature but 
#                             any temp will do
# 
# Outputs: GPP - Gross Primary Productivity, numpy array
#                                            (same dims as NEE)
#
##############################################################
def GPP_from_NEE_SW_T( NEEin, SWin, Tin,                         \
                       SWmin=1., SW_smooth_int=0,                \
                       spike_filter=None,fill_value=None,        \
                       Gsigma=None, Gaxis=-1,                    \
                       MED_filter_size=5,                        \
                       STDmask_window=20,STDmask_limit=None,     \
                       fit_function='Q10', FIT_maxfev=1500,      \
                       T_sf=0.1, T_0=10.,                        \
                       ReturnResp=False,ReturnNEEfiltered=False, \
                       ReturnFitParams=False,                    \
                       ):
    # Make copies of input arrays to work with
    # Unnecessary so removed
    NEEwk = NEEin#.copy()
    SWwk  = SWin#.copy()
    Twk   = Tin#.copy()
    
    # set fill value if not provided
    if fill_value==None:
        try:
            fill_value=NEEwk.fill_value
        except:
            fill_value=1e20

    #print(fill_value)
    
    # Ensure that arrays are masked arrays
    NEEwk = np.ma.masked_invalid(NEEwk)
    SWwk  = np.ma.masked_invalid(SWwk)
    Twk   = np.ma.masked_invalid(Twk)

    # Calculate Dark respiration and  GPP
    # Dark Respiration is calculated by fitting the Night-time data
    #    to a Q10 function
    # Night-time data is defined as data which is not within 2 timesteps
    #    of registered SWdown
    
    # Daytime mask for calculation of GPP
    DAYTIME_mask = SWwk>SWmin
    DAYTIME_mask_copy=DAYTIME_mask.copy()
    # add points within 2 timesteps of daytime to mask
    # copy mask to prevent masking all elements
    for iSMOOTH in range(SW_smooth_int,len(DAYTIME_mask)-SW_smooth_int):
        if any(DAYTIME_mask_copy[iSMOOTH-SW_smooth_int:iSMOOTH+SW_smooth_int]==True):
            DAYTIME_mask[iSMOOTH]=True
    del DAYTIME_mask_copy
    
    if spike_filter==None:
        NEEwk=NEEwk
    elif spike_filter.lower()=='gaussian':
        # Remove NEE spikes with a gaussian filter
        # create temporary copy of NEEwk
        temp = NEEwk.copy()
        temp[NEEwk.mask==True]=np.nan
        
        if Gsigma==None:
            Gsigma=np.std(NEEwk)
        
        NEE_filtered = im.gaussian_filter1d(temp,sigma=Gsigma,axis=Gaxis)
        del temp
        # Reapply original mask 
        NEE_filtered[NEEwk.mask==True]=fill_value
        # Mask out elements which became nan during filtering
        NEE_filtered[NEE_filtered!=NEE_filtered]=fill_value
        NEE_filtered=np.ma.masked_equal(NEE_filtered,fill_value)
        
        #update NEEwk
        NEEwk = NEE_filtered
        del NEE_filtered
    elif spike_filter.lower()=='median':
        # Remove NEE spikes with a median filter
        # create temporary copy of NEEwk
        temp = NEEwk.copy()
        temp[NEEwk.mask==True]=np.nan
        # 
        NEE_filtered = im.median_filter(temp,MED_filter_size)
        # Reapply original mask 
        NEE_filtered[NEEwk.mask==True]=fill_value
        # Mask out elements which became nan during filtering
        NEE_filtered[NEE_filtered!=NEE_filtered]=fill_value
        NEE_filtered=np.ma.masked_equal(NEE_filtered,fill_value)
        
        #update NEEwk
        NEEwk = NEE_filtered
        del NEE_filtered
    elif spike_filter.lower()=='stddev_mask':
        # Mask out spikes using the ECP stddev_mask filter
        stddev_mask = ST.STDDEV_mask(NEEwk.data,window=STDmask_window,std_limit=STDmask_limit)
        full_mask = stddev_mask|NEEwk.mask
        NEEwk=np.ma.masked_array(NEEwk,mask=full_mask)
        NEEwk.data[NEEwk.mask==True]=fill_value
    else:
        print('Unrecognised Filter, no filter applied')
        NEEwk=NEEwk
    
    # create masked arrays for dark_NEE and dark_Tair based on Daytime mask and data masks 
    dark_mask     = DAYTIME_mask|NEEwk.mask|Twk.mask
    dark_NEE_data = np.ma.masked_array(NEEwk.data.copy(),\
                                       mask=dark_mask,\
                                       fill_value=fill_value)
    dark_NEE_data.data[dark_mask==True]=fill_value    
    dark_T_data= np.ma.masked_array(Twk.copy(),\
                                    mask=dark_mask,\
                                    fill_value=fill_value)
    dark_T_data.data[dark_mask==True]=fill_value
    
    dark_Resp_data = dark_NEE_data.copy()
    full_mask = (dark_T_data.mask==False)&(dark_NEE_data.mask==False)

    # Optimise the Q10 fit parameter for local conditions by fitting:
    #    DarkResp = Q10^( m(Tair+c) ), 
    # the parameter fit is performed on rolling time frame covering the encompassing 48 hours
    if fit_function=='Q10':
        Flambda = lambda T,Q_10,Resp_10: MT.Q10_func(T,Q_10,Resp_10,T_sf,T_0)
        popt,pcov=curve_fit(Flambda,\
                            dark_T_data[full_mask],\
                            dark_Resp_data[full_mask],\
                            p0=(2.,np.mean(dark_Resp_data)),maxfev=FIT_maxfev )
        Resp_data = Flambda(Twk,popt[0],popt[1])
        FitParams = (popt,pcov)
    else:
        print('Unrecognised Fit Function, fitting to a Q10 function')
        Flambda = lambda T,Q_10,Resp_10: MT.Q10_func(T,Q_10,Resp_10,T_sf,T_0)
        popt,pcov=curve_fit(Flambda,\
                            dark_T_data[full_mask],\
                            dark_Resp_data[full_mask],\
                            p0=(2.,np.mean(dark_Resp_data)),maxfev=FIT_maxfev )
        Resp_data = Flambda(Twk,popt[0],popt[1])
        FitParams = (popt,pcov)
    
    Resp_data = np.ma.masked_invalid(Resp_data)
    Resp_data[Resp_data.mask==True]=fill_value
    Resp_data = np.ma.masked_equal(Resp_data.data,fill_value)
    #Resp_data.fill_value=fill_value

    GPP_data  = Resp_data-NEEwk
    GPP_data  = np.ma.masked_equal(GPP_data.data,fill_value)

    output = [GPP_data]
    if ReturnResp:
        output.append(Resp_data)

    if ReturnNEEfiltered:
        output.append(NEEwk)

    if ReturnFitParams:
        output.append(FitParams)

    return output


