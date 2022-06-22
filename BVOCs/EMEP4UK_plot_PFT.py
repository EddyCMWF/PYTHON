#!/usr/bin/python
#
# Python module to plot the hwsd dat on EMEP grid
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import argparse
import netCDF4 as nc
import netcdftime as ncdt
import matplotlib.pyplot as plot
import plot_tools as PT
import EMEP4UK_tools as ET
#
###################################################################################################################
# Define class
###################################################################################################################
class plotLAIPFT:
    def parse_input(self):
        #
        parser=argparse.ArgumentParser(description='Extract a subset of a binary file to a netcdf file')
        # optional
        parser.add_argument('--infileUK',type=str,help='Input file UK',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_MPI_DougDiag.monthly_wetl.nc')
        parser.add_argument('--infileEU',type=str,help='Input file EU',required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/JULES_output/JULES_EMEP4UK_EUROPE_MPI_DougDiag.monthly_wetl.nc')
        parser.add_argument('--outdir',type=str,help='output directory' ,required=False, \
                             default='/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_Output/' )
        #
        # positional
        #
        # Parse the arguments
        args=parser.parse_args()
        #
        return args.infileUK, args.infileEU, args.outdir
#
###################################################################################################################
# Define functions
###################################################################################################################
# Round to given significant figures
def round2SignifFigs(vals,n):
    mags = 10.0**np.floor(np.log10(np.abs(vals)))  # order of mag's
    outvals = np.around(vals/mags,n-1)*mags             # round(val/omag)*omag
    try:
        outvals[np.where(np.isnan(vals))] = 0.0           # where order of mag = 0, set to zero
    except:
        if np.isnan(outvals):
            outvals=0.0
    #
    return outvals
#
###################################################################################################################
# Define Main Program
###################################################################################################################

if __name__=='__main__':
    #
    psm=plotLAIPFT()
    infileUK,infileEU,outdir =psm.parse_input()
    #
    PFT_infileUK = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_LandFrac.nc.noIAM.oneSOIL'
    PFT_infileEU = '/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_LandFrac.nc.noIAM.oneSOIL.binaryICE'
    outdir   = '/users/eow/edwcom/EMEP/EMEP4UK/plots/VEG_params/'
    EMEP4UK_2D_latlonfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_Landuse.nc'
    EMEP4UK_2D_latname='lat'
    EMEP4UK_2D_lonname='lon'
    EMEP4UK_EUROPE_2D_latlonfile='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_EUROPE_Landuse.nc'
    EMEP4UK_EUROPE_2D_latname='lat'
    EMEP4UK_EUROPE_2D_lonname='lon'
    #
    PFT_varname='Land_Frac'
    #
    # Fetch 2D lat lon grids:
    inf=nc.Dataset(EMEP4UK_2D_latlonfile,'r')
    lats1_2D=inf.variables[EMEP4UK_2D_latname][:]
    lons1_2D=inf.variables[EMEP4UK_2D_lonname][:]
    inf.close()
    inf=nc.Dataset(EMEP4UK_EUROPE_2D_latlonfile,'r')
    lats2_2D=inf.variables[EMEP4UK_EUROPE_2D_latname][:]
    lons2_2D=inf.variables[EMEP4UK_EUROPE_2D_lonname][:]
    inf.close()
    #
    # Read in PFT grids
    inf = nc.Dataset(PFT_infileUK)
    PFT_UK = inf.variables[PFT_varname][:]
    inf.close()
    inf = nc.Dataset(PFT_infileEU)
    PFT_EU = inf.variables[PFT_varname][:]
    inf.close()
    #
    # Calculate modal value for each pixel:
    PFT_UK_max = np.argmax(PFT_UK,axis=0)+1
    PFT_EU_max = np.argmax(PFT_EU,axis=0)+1
    
    #
    CLEVELS   = range(1,17)
    COLOURS   = [ '#228B22','#556B2F','#ADFF2F','#7CFC00', \
                  '#F0E68C','#FFD700','#DAA520',           \
                  '#98FB98','#00FA9A','#FFA500',           \
                  '#8B4513','#4169E1','#FFFFF0','#FF0000', \
                  '#808080']    
    #
    #
    tick_levels=list(np.arange(1.5,15,1.))
    # 
    #
    latrange = [41.,38.]
    lonrange = [-35.,60.]
    CU_lonrange=[-15,12]
    CU_latrange=[50.,56.]
    #
    #
    print 'plotting PFT maps:'
    cbar_title = 'PFT number'
    #
    plot_title = 'Maximum Grid-square Land-Cover'
    #
    # Plot Europe Data
    plotname=outdir+'EMEP4UK_EUROPE_max_PFT.png'
    PT.plot_map(PFT_EU_max,lons2_2D,lats2_2D,                             \
                CLEVELS=CLEVELS, COLOURS=COLOURS,                         \
                MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',             \
                WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],               \
                PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,             \
                iDISPLAY='N', FILE_PLOT=plotname,                         \
                LATDEL=15., LONDEL=15., RESOLUTION='h',                   \
                PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange )
    #
    plotname=outdir+'EMEP4UK_max_PFT.png'
    PT.plot_map(PFT_EU_max,lons2_2D,lats2_2D,                             \
                DATA2=PFT_UK_max,LONS2=lons1_2D,LATS2=lats1_2D,           \
                CLEVELS=CLEVELS, COLOURS=COLOURS,                         \
                MAP_TYPE='Mesh', CBAR_ORIENTATION='vertical',             \
                WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],               \
                PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,             \
                iDISPLAY='N', FILE_PLOT=plotname,                         \
                LATDEL=15., LONDEL=15., RESOLUTION='h',                   \
                PROJECTION='stere', LON_RANGE=CU_lonrange,LAT_RANGE=CU_latrange )
     #

