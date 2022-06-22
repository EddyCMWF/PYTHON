#!/usr/bin/python2.7

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as col
import numpy as np



projection='cyl'
res='h'
latdel=15
londel=15
latrange=[-90,90]
lonrange=[-180,180]

lonrange=[-15,15]
latrange=[50.,55.]
latdel=5


londel=5

lat_cen = 90.
lon_cen = -32.


M = Basemap(projection=projection, \
            resolution=res,lat_0=lat_cen,lon_0=lon_cen, \
            llcrnrlat =latrange[0],urcrnrlat=latrange[1], \
            llcrnrlon =lonrange[0],urcrnrlon=lonrange[1]  )

#M.bluemarble()
M = Basemap(projection=projection, \
            resolution=res, \
            llcrnrlat =latrange[0],urcrnrlat=latrange[1], \
            llcrnrlon =lonrange[0],urcrnrlon=lonrange[1]  )
            
M.etopo()
M.drawrivers(linewidth =0.1,color='aqua')
M.drawcountries(linewidth=1.)
M.drawcoastlines(linewidth=1.)


parallels = np.arange(-90.,90.,latdel)
M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 1.,fontsize=10)

meridians = np.arange(-180.,180.,londel)
M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 1.,fontsize=10)



plt.show()


