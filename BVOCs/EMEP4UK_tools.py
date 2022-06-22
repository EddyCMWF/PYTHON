#
# Python Module containing tools for dealing with EMEP data
# 
# Edward Comyn-Platt
# CEH Wallingford
# Jan 2015
#
##########################################################
import numpy as np
import netCDF4 as nc
#

def lonlat_to_xy(lon,lat,EMEP_proj='50'):
    Re    = 6370.         # radius of Earth
    lat_0 = 60.           # defining latitude
    lon_0 = -32.          # rotation angle
    #
    if (EMEP_proj == '50'):
        xpol=8.           # x coord of pole
        ypol=110.         # y coord of pole
        d=50.             # grid length at defining latitude
    elif (EMEP_proj == '150'):
        xpol=3.           # x coord of pole
        ypol=37.          # y coord of pole
        d=150.            # grid length at defining latitude
    else:
        raise Exception('Unrecognised EMEP projection')
    
    M = (Re/d) * ( 1. + np.sin( np.deg2rad(lat_0) ) )
    
    x = xpol + M * np.tan( (np.pi/4.) - (np.deg2rad(lat)/2.) ) * np.sin(np.deg2rad(lon-lon_0))
    y = ypol - M * np.tan( (np.pi/4.) - (np.deg2rad(lat)/2.) ) * np.cos(np.deg2rad(lon-lon_0))
    
    return x, y
    

def xy_to_lonlat(x,y,EMEP_proj='50'):
    Re    = 6370.         # radius of Earth
    lat_0 = 60.           # defining latitude
    lon_0 = -32.          # rotation angle
    #
    if (EMEP_proj == '50'):
        xpol=8.           # x coord of pole
        ypol=110.         # y coord of pole
        d=50.             # grid length at defining latitude
    elif (EMEP_proj == '150'):
        xpol=3.            # x coord of pole
        ypol=37.           # y coord of pole
        d=150.             # grid length at defining latitude
    else:
        raise Exception('Unrecognised EMEP projection')
    
    M = (Re/d) * ( 1. + np.sin( np.deg2rad(lat_0) ) )
    r = np.sqrt( ((x-xpol)**2.) + ((y-ypol)**2.) )
    
    lat = 90. - ( (360./np.pi) * np.arctan( r/M ) )
    lon = lon_0 + ( (180./np.pi) * np.arctan( (x-xpol)/(ypol-y) ) )
    
    # Correct for points below the pole coord where the rotation angle should be lon_0+180
    lon[y>ypol]   = (lon[y>ypol]+180.)
    # Correct for points which are gt 180., i.e. change to negative
    lon[ lon>(180.) ] = lon[ lon>(180.) ] - 360.
    # correct for pole point
    lon[ (x==xpol) & (y==ypol) ] = 0.
    
    return lon,lat


