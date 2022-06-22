#!/bin/env python

import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from matplotlib.cm import get_cmap
from matplotlib.pyplot import sca

__author__ = "Phil Harris"
__email__ = "philip.p.harris@gmail.com"


class Mapper(object):
    def __init__(self,ncg,nrg,dlon,dlat,lon0,lat0,filename):
        self.NCG = ncg
        self.NRG = nrg
        self.DLON = dlon
        self.DLAT = dlat
        self.LON0 = lon0
        self.LAT0 = lat0
        self.land = self.readLand(filename)
        self.shape = (self.NRG,self.NCG)
        return

    def __len__(self):
        return len(self.land)

    def readLand(self,filename):
        """Return the global domain mapping file info as a NumPy array of
        void objects, one void per land point.
        """
        z = np.genfromtxt(filename,
                          dtype=(int, int, int, int, float, float, int, int, int),
                          skip_header=1, delimiter=',')
        return z

    def unpack(self,x,world=False):
        """Returns data stored on a vector of land points, x, unpacked to a 2D
        array subject to the list of land point mapping info, land.
        """
        if world:
            c, r = zip(*[(l[2]-1,l[3]-1) for l in self.land])
            outShape = self.shape
        else:
            c, r = zip(*[(l[-2]-1,l[-1]-1) for l in self.land])
            outShape = [max(r)-min(r)+1,max(c)-min(c)+1]
        xy = np.ma.masked_all(outShape,dtype=x.dtype)
        xy[r,c] = x[:]
        return xy

    def getMesh(self):
        """Return an nc+1 by nr+1 mesh of longitudes and latitudes.  The +1 is
        so that the extremes of the meshes can be passed to basemap() as
        lower-left and upper-right corners of the grid."""
        lon = (np.arange(self.NCG+1)-0.5)*self.DLON + self.LON0
        lat = (np.arange(self.NRG+1)-0.5)*self.DLAT + self.LAT0
        XX, YY = np.meshgrid(lon, lat)
        return XX, YY

    def getMap(self,subsection=None,
               meridians=np.arange(-180,180.1,30),
               parallels=np.arange(-90, 90.1,30),
               resolution='l',fontsize=10):
        """Return a regular lon/lat (cylindrical projection) Basemap instance
        covering the domain specified by the arguments.

        Keyword arguments:
        subsection -- if set, returns a map showing a subsection of the full
                      domain defined by this Mapper instance (default None).
                      The "subsection" argument is a list of indices on the
                      full domain, [x0,y0,x1,y1].
        meridians -- list of longitudes passed to Basemap.drawmerians().
        parallels -- list of latitudes passed to Basemap.drawparallels().
        resolution -- resolution of the map coastlines (default "l").
        fontsize -- size of the meridian/parallel tick labels (default 10 pt).
        """
        XX, YY = self.getMesh()
        if subsection is None:
            llcrnrlon, llcrnrlat = XX[0,0], YY[0,0]
            urcrnrlon, urcrnrlat = XX[-1,-1], YY[-1,-1]
        else:
            s = subsection
            llcrnrlon, llcrnrlat = XX[s[1],s[0]], YY[s[1],s[0]]
            urcrnrlon, urcrnrlat = XX[s[3],s[2]], YY[s[3],s[2]]
        M = bm(projection='cyl', resolution=resolution, area_thresh=1.0e4,
               llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
               urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        M.drawcoastlines(linewidth=0.25, color='k')
        M.drawmeridians(meridians,linewidth=0.1,labels=[1,0,0,1],fontsize=fontsize)
        M.drawparallels(parallels,linewidth=0.1,labels=[1,0,0,1],fontsize=fontsize)
        return XX, YY, M

    def plotMap(self,ax,z,world=False,subsection=None,**kwargs):
        """Plots a map to axes, ax, and then plots a 2D array of data to the 
        same axes using pcolormesh().
        """
        keys_map = ["meridians","parallels","resolution","fontsize"]
        kw_map = dict([(k,v) for k, v in kwargs.iteritems() if k in keys_map])

        keys_cmesh = ["cmap","vmin","vmax"]
        kw_cmesh = dict([(k,v) for k, v in kwargs.iteritems() if k in keys_cmesh])

        sca(ax)
        if world:
            X, Y, M = self.getMap(subsection=subsection,**kw_map)
        else:
            X, Y, M = self.getMap(subsection=subsection,**kw_map)
        I = M.pcolormesh(X,Y,z,**kw_cmesh)
        return I

    def addHatching(self,A,X,Y,M,dx,dy,hatch="/",alpha=1,lw=0,fill=False,snap=False):
        """Adds gridbox-based hatching to an existing map plot.  This is a
        manual hack for Matplotlib pre-v1.3.1, so don't expect it to be fast.

        A is a Matplotlib Axes object, which should have already been set up
        using Mapper.getMap(), Mapper.plotMap() or by calling Basemap()
        manually.

        X, Y are [nlat+1,nlon+1] arrays of gridbox edges, the same thing that
        is returned from Mapper.getmesh() and given to pcolormesh() or imshow().

        M is a [nlat,nlon] array of boolean, the sort you'd find attached to a
        Numpy masked array.  Hatching is applied where M==True, i.e., for masked
        grid boxes.
        """
        from matplotlib.patches import Rectangle
        for xll, yll, masked in zip(X[:-1,:-1].flat, Y[:-1,:-1].flat, M.flat):
            if masked:
                pp = Rectangle((xll,yll),dx,dy,hatch=hatch,alpha=alpha,lw=lw,
                               fill=fill,snap=snap)
                A.add_patch(pp)
        return

def main():
    print "Nothing to see here."

if __name__ == "__main__":
    main()
