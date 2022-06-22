#!/usr/bin/python

"""
Class for mapping between GSWP2 global land grid and the Euro sub-domain
used for SWELTER work.
"""

import numpy as np
from mpl_toolkits.basemap import Basemap as bm

__author__ = "Phil Harris"
__email__ = "philip.p.harris@gmail.com"

MAPPING_FILE = '/users/global/ppha/swelter/modis/L3_interim/princeton_map.asc'

NC_EURO = 25
NR_EURO = 19
NP_EURO = NC_EURO*NR_EURO

LON0_EURO = -9.5
LAT0_EURO = 41.5
DLON_EURO = 1.0
DLAT_EURO = 1.0

NC_GSWP = 360
NR_GSWP = 180
NP_GSWP = NC_GSWP*NR_GSWP

LON0_GSWP = -179.5
LAT0_GSWP = 90.0
DLON_GSWP = 1.0
DLAT_GSWP = -1.0

MDI = np.float64(-1.0E20)

class Swmap:
    def readGswpLand(self):
        """Return the GSWP2 domain mapping file info as a NumPy array of
        void objects, one void per GSWP2 land point.
        """
        z = np.genfromtxt(MAPPING_FILE,
                          dtype=(int, int, int, int, float, float, int, int, int),
                          skip_header=1, delimiter=',')
        return z

    def readEuroLand(self):
        """Return the Euro sub domain mapping as a dictionary of
        GSWP2 index: (Euro column, Euro row) values.
        """
        gswp = self.readGswpLand()
        euro = [p for p in gswp if p[-3] > 0]
        return euro

    def unpack(self,dgrid,dx,v=-1):
        """Unpack a dictionary of data, with pixel numbers as keys, to a
        2D numpy array for plotting as a map.
        """
        z = np.zeros((NR_EURO,NC_EURO))
        z[:,:] = MDI
        for p in dx.keys():
            c, r = dgrid[p]
            z[r-1,c-1] = dx[p][v]
        return np.ma.masked_less(z,-1.0E19)

    def getMesh(self, nc=NC_EURO, nr=NR_EURO,
                lon0=LON0_EURO, lat0=LAT0_EURO, dlon=DLON_EURO, dlat=DLAT_EURO):
        """Return an nc+1 by nr+1 mesh of longitudes and latitudes for use with pcolor."""
        lon = (np.arange(nc+1)-0.5)*dlon + lon0
        lat = (np.arange(nr+1)-0.5)*dlat + lat0
        XX, YY = np.meshgrid(lon, lat)
        return XX, YY

    def getMap(self, nc=NC_EURO, nr=NR_EURO,
               lon0=LON0_EURO, lat0=LAT0_EURO, dlon=DLON_EURO, dlat=DLAT_EURO):
        XX, YY = self.getMesh(nc=nc, nr=nr,
                              lon0=lon0, lat0=lat0, dlon=dlon, dlat=dlat)
        M = bm(projection='cyl', resolution='l', area_thresh=1.0e4,
               llcrnrlon=XX[0,0], llcrnrlat=YY[0,0],
               urcrnrlon=XX[-1,-1], urcrnrlat=YY[-1,-1])
        M.drawcoastlines(linewidth=0.5, color='k')
        M.drawmeridians(np.arange(-180,180,5),linewidth=0.1,labels=[1,0,0,1],fontsize=8)
        M.drawparallels(np.arange( -90, 90,5),linewidth=0.1,labels=[1,0,0,1],fontsize=8)
        return XX, YY, M

def main():
    print "Nothing to see here."

if __name__ == "__main__":
    main()

##
## End of script.
##
