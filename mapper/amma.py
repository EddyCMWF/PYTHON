#!/usr/bin/python

"""
Class of useful stuff for plotting AMMA assimilation inputs and outputs.
"""

__author__ = "Phil Harris"
__email__ = "philip.p.harris@gmail.com"

import numpy as np
from mpl_toolkits.basemap import Basemap as bm
import netCDF4 as nc4

class Epsat:
    NCOL = 500
    NROW = 250
    FIL_MAPPING = "/users/global/ppha/amma/assimilation/prep/ancil/epsat/assim_epsat_map.asc"
    DIR_DAILY = "/prj/amma/ppha/epsat"

    def getMapping(self):
        """Return the ASSIM-to-EPSAT mapping info."""
        z = np.loadtxt(self.FIL_MAPPING,dtype=np.int32)
        return z

    def getDaily(self, date):
        import os
        z = None
        d = self.DIR_DAILY + "/%4.4i/%2.2i/" % (date.year, date.month)
        f = "surf-rr_epsat-sg_multi-sat_010d_30min_%4.4i%2.2i%2.2i_v3.1-02c.nc" % (date.year, date.month, date.day)
        fnom = os.path.join(d,f)
        if os.path.exists(fnom):
           ncid = nc4.Dataset(fnom,'r',format='NETCDF3_CLASSIC')
           z = ncid.variables["data"][:,:,:]
        return z

class Assim:
    NCOL = 728
    NROW = 348
    FIL_ASSIM_LON = "/users/global/ppha/amma/assimilation/prep/ancil/assim_grid/assim_lon.bin"
    FIL_ASSIM_LAT = "/users/global/ppha/amma/assimilation/prep/ancil/assim_grid/assim_lat.bin"
    FIL_SOIL_ALB = "/users/global/ppha/amma/assimilation/prep/ancil/make_modis/assim_modis_alb.bin"
    FIL_TILE_FRAC = "/prj/amma/ppha/assim/ecoclimap/msg_map_pix_frac.gra"

    def getSoilAlb(self):
        """Return the assim grid soil albedo as a numpy array."""
        z = np.fromfile(self.FIL_SOIL_ALB,dtype=">f4").reshape(self.NROW,self.NCOL)
        z = np.ma.masked_less(z,0.0)
        return z

    def getTileFracs(self):
        """Return the JULES tile fractions as a numpy array."""
        z = np.fromfile(self.FIL_TILE_FRAC,dtype=">f4").reshape(9,self.NROW,self.NCOL)        
        return z

    def getMesh(self):
        """Return meshgrid like arrays of assim lon, lats for use with imshow()."""
        XX = np.fromfile(self.FIL_ASSIM_LON,dtype='>i4').reshape(self.NROW,self.NCOL) * 1.0e-2
        YY = np.fromfile(self.FIL_ASSIM_LAT,dtype='>i4').reshape(self.NROW,self.NCOL) * 1.0e-2
        return XX, YY

    def getMap(self, corners=(0,0,-1,-1), XX=None, YY=None):
        """
        Return a geostationary Basemap object suitable for plotting assim grid output.
        The corners argument is a tuple (llcol,llrow,urcol,urrow), which refers to
        columns and rows on the full ASSIM domain."""

        if 0 < corners[3] <= corners[1] or 0 < corners[2] <= corners[0]:
            raise ValueError("Your ASSIM sub-section corners are wrong.")

        if XX == None or YY == None:
            XXU, YYU = self.getMesh()
        else:
            XXU, YYU = XX, YY
        llcrnrlon = XXU[corners[1],corners[0]]
        llcrnrlat = YYU[corners[1],corners[0]]
        urcrnrlon = XXU[corners[3],corners[2]]
        urcrnrlat = YYU[corners[3],corners[2]]

        M = bm(projection='geos',resolution='i',area_thresh=0.1e3, 
               lon_0=0.0,
               llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
               urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        M.drawcountries()
        M.drawmapboundary()
        M.drawmeridians(np.arange(-180,180,2),linewidth=0.2,labels=[1,0,0,1])
        M.drawparallels(np.arange( -90, 90,2),linewidth=0.2,labels=[1,0,0,1])
        return M

if __name__ == "__main__":
    print "Hello."

