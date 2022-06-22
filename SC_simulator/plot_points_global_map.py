#!/usr/bin/python2.7

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as col
import numpy as np



pt_file = '/users/eow/edwcom/SC_simulator/subdomains/SC_pointnames.dat'
pt_lat,pt_lon,pt_name = [], [], []
for line in open(pt_file).readlines():
    split_line=line.split()
    pt_lat.append(float(split_line[0]))
    pt_lon.append(float(split_line[1]))
    pt_name.append(' '.join(split_line[2:]))


projection='cyl'
res='c'
latdel=15
londel=15
latrange=[-90,90]
lonrange=[-180,180]


#M.bluemarble()
M = Basemap(projection=projection, \
            resolution=res, \
            llcrnrlat =latrange[0],urcrnrlat=latrange[1], \
            llcrnrlon =lonrange[0],urcrnrlon=lonrange[1]  )
            
#M.etopo()
#M.drawrivers(linewidth =0.1,color='aqua')
#M.drawcountries(linewidth=1.)
M.drawcoastlines(linewidth=1.)
M.fillcontinents(color='lightgrey')

parallels = np.arange(-90.,90.,latdel)
M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 1.,fontsize=10)

meridians = np.arange(-180.,180.,londel)
M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 1.,fontsize=10)

for pt in range(len(pt_name)):
    x,y=M(pt_lon[pt],pt_lat[pt])
    M.plot(x,y,marker='*',color='red',ms=10)
    plt.text(x+2,y,pt_name[pt],fontsize=15,weight='bold')

plt.show()


