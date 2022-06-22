#!/usr/bin/env python
##################################################################################
#
# Program: remap_scrip.py   
# Author: Edward Comyn-Platt, 02/2013
#
# Purpose: Set of tools relevent to remapping data
# 
##################################################################################

import numpy as np
import sys
#
##################################################################################
# SCRIP class: Tools relevent to SCRIP remapping
##################################################################################
#class SCRIP(object):
#########################################################################
# Consevative Remapping
#########################################################################
# INPUT
###########
# Required:
# map_wts - map weights, array of weights of src data on to dst array
#           Dimensions - [ num_links, orders ]
# src_array - source data, i.e. the data to be remapped onto a new grid
#              1 dimensional array of data
# src_add - address/index of source data for each map weight
# dst_add - address/index of destination data for each map weight
############
# Opitional:
# dst_grid_size - size of destination grid (dst_array), if not provided
#                 dst_array length is maximum value in dst_add
# order - order of remapping, default is 1st order, 
#          additional data is required for 2nd and third order remapping
# STDDEV - standard deviation flag, set to True to output a 
#            weighted standard deviation
# out_mv - value to assign to missing data in the dst_array
# 
#########################################################################
def conserv_remap(map_wts,src_array,src_add,dst_add, \
                  dst_grid_size=False, \
                  order=1,STDDEV=False, out_mv=-9999.0):
    
    # if dst_grid_size not set by user dst_grid_size is set to maximum dst_add
    if not (dst_grid_size):
        print 'setting dst_grid_size = ', str(np.amax(dst_add))
        dst_grid_size=np.amax(dst_add)
    
    # if dst_grid_size is less than maximum dst_add flag up and modifiy dst_grid_size
    if (dst_grid_size<np.amax(dst_add)):
        print 'dst_grid_size not long enough for dst_add values'
        print 'setting dst_grid_size = '+ str(np.amax(dst_add))
        dst_grid_size=np.amax(dst_add)
    
    dst_array=np.zeros(dst_grid_size)+out_mv
    if (STDDEV):
        dst_std_array = np.zeros(dst_grid_size)+out_mv
        
    # Array of unique destination points:
    uni_dst_add=set( dst_add )
    
    # Create list of indexes for each dst point
    uber_index = [ np.where(dst_add==idx) for idx in uni_dst_add ]
    
    if (order==1):
        print 'Performing conservative 1st order remapping\n'
        
        # loop through list of indexes for each uniques destination address and 
        # compute the area weighted mean based on that index of points and map_weights
        temp_array=  np.array([ np.sum( src_array[src_add[idx]]*map_wts[idx,0])  \
                               / np.sum( map_wts[idx,0] ) for idx in uber_index ])
        # put calculated values in appropriate slots of destination array
        dst_array[list(uni_dst_add)]=temp_array
        del temp_array
        
        # Repeat for standard deviation if neccessary
        if (STDDEV):
            # STD equation is quite complicated, refer to ECP for further info.
            # Basically, the methodology is not full proof and weighted STDs are largely undecided methods
            temp_std_array = np.array( [ np.sqrt( np.sum( map_wts[dst_add==n,0]    \
                                                          * ((src_array[src_add[dst_add==n]]-dst_array[n])**2) )  \
                                                  * np.sum( map_wts[dst_add==n,0] )  \
                                                  / ( (np.sum( map_wts[dst_add==n,0] )**2)-np.sum(map_wts[dst_add==n,0]**2) )  ) \
                                         for n in uni_dst_add ] )
            dst_std_array[uni_dst_add]=temp_std_array
    
    if (STDDEV):
        output=[dst_array,dst_std_array]
    else:
        output=dst_array
        
    return output
        
