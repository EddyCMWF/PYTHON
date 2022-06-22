#!/bin/env python

import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from matplotlib.cm import get_cmap
from matplotlib.pyplot import sca

__author__ = "Phil Harris"
__email__ = "philip.p.harris@gmail.com"

MAPPING_FILE = "/users/global/ppha/swelter/um_amip/prep/ancils/amip_map.asc"

NC_AMIP = 192
NR_AMIP = 144
NP_AMIP = NC_AMIP*NR_AMIP
NL_AMIP = 10663

DLON_AMIP = 1.8750
DLAT_AMIP = 1.2500
LON0_AMIP = 0.9375
LAT0_AMIP = -89.375


NC_EURO = 23
NR_EURO = 20
NP_EURO = NC_EURO*NR_EURO
NL_EURO = 350

DLON_EURO = 1.8750
DLAT_EURO = 1.2500
LON0_EURO = 349.6875 - 360.0
LAT0_EURO = 35.6250



class Swmap:
    def readAmipLand(self):
        """Return the AMIP domain mapping file info as a NumPy array of
        void objects, one void per AMIP land point.
        """
        z = np.genfromtxt(MAPPING_FILE,
                          dtype=(int, int, int, int, float, float, int, int, int),
                          skip_header=1, delimiter=',')
        return z

    def readEuroLand(self):
        """Return the Euro sub domain mapping as a dictionary of
        AMIP index: (Euro column, Euro row) values.
        """
        amip = self.readAmipLand()
        euro = [p for p in amip if p[-3] > 0]
        return euro

    def unpack(self,x,landmap,world=False):
        """Returns data stored on a vector of land points, x, unpacked to a 2D
        array subject to the list of land point mapping info, land.
        """
        if world:
            c, r = zip(*[(l[2]-1,l[3]-1) for l in landmap])
            outShape = [NR_AMIP,NC_AMIP]
        else:
            c, r = zip(*[(l[-2]-1,l[-1]-1) for l in landmap])
            outShape = [NR_EURO,NC_EURO]
        xy = np.ma.masked_all(outShape,dtype=x.dtype)
        xy[r,c] = x[:]
        return xy

    def getMesh(self, nc=NC_EURO, nr=NR_EURO,
                lon0=LON0_EURO, lat0=LAT0_EURO, dlon=DLON_EURO, dlat=DLAT_EURO):
        """Return an nc+1 by nr+1 mesh of longitudes and latitudes.  The +1 is
        so that the extremes of the meshes can be passed to basemap() as
        lower-left and upper-right corners of the grid."""
        lon = (np.arange(nc+1)-0.5)*dlon + lon0
        lat = (np.arange(nr+1)-0.5)*dlat + lat0
        XX, YY = np.meshgrid(lon, lat)
        return XX, YY

    def getMap(self, nc=NC_EURO, nr=NR_EURO,
               lon0=LON0_EURO, lat0=LAT0_EURO, dlon=DLON_EURO, dlat=DLAT_EURO,
               meridians=np.arange(-180,180.1,5),parallels=np.arange( -90, 90.1,5),
               resolution='l',fontsize=10):
        """Return a regular lon/lat (cylindrical projection) Basemap instance
        covering the domain specified by the arguments.  By default this is
        the SWELTER Euro2 subsection of the AMIP domain."""
        XX, YY = self.getMesh(nc=nc, nr=nr,
                              lon0=lon0, lat0=lat0, dlon=dlon, dlat=dlat)
        M = bm(projection='cyl', resolution=resolution, area_thresh=1.0e4,
               llcrnrlon=XX[0,0], llcrnrlat=YY[0,0],
               urcrnrlon=XX[-1,-1], urcrnrlat=YY[-1,-1])
        M.drawcoastlines(linewidth=0.25, color='k')
        M.drawmeridians(meridians,linewidth=0.1,labels=[1,0,0,1],fontsize=fontsize)
        M.drawparallels(parallels,linewidth=0.1,labels=[1,0,0,1],fontsize=fontsize)
        return XX, YY, M

    def getGlobalMap(self,meridians=np.arange(0,360.1,45),
                     parallels=np.arange( -90, 90.1,15),
                     resolution='l',fontsize=10):
        """Return a regular lon/lat (cylindrical projection) Basemap instance
        covering global AMIP domain."""
        XX, YY, M = self.getMap(nc=NC_AMIP,nr=NR_AMIP,
                                lon0=-179.0625,lat0=LAT0_AMIP,
                                dlon=DLON_AMIP,dlat=DLAT_AMIP,
                                meridians=meridians,parallels=parallels,
                                resolution=resolution,fontsize=fontsize)
        return XX, YY, M

    def plotMap(self,ax,z,world=False,interpolation="nearest",**kwargs):
        """Plots a map to axes, ax, and then plots a 2D array of data to the 
        same axes using imshow().
        """
        keys_map = ["meridians","parallels","resolution","fontsize"]
        kw_map = dict([(k,v) for k, v in kwargs.iteritems() if k in keys_map])

        keys_imshow = ["cmap","vmin","vmax"]
        kw_imshow = dict([(k,v) for k, v in kwargs.iteritems() if k in keys_imshow])

        sca(ax)
        if world:
            X, Y, M = self.getGlobalMap(**kw_map)
        else:
            X, Y, M = self.getMap(**kw_map)
        I = M.imshow(z,interpolation=interpolation,**kw_imshow)
        return I

def main():
    print "Nothing to see here."

if __name__ == "__main__":
    main()
