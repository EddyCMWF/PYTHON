#!python

import netCDF4 as nc
import numpy as np
import pylab as plt
import plot_tools as PT
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle,Rectangle

file='/prj/wetlands_africa/jules/JASMIN/WFD_EI_global/MPI_WFD_EI_global.monthly.nc'

inf=nc.Dataset(file)

for var in inf.variables:
    print(str(var))

lats=inf.variables['latitude'][:].flatten()
lons=inf.variables['longitude'][:].flatten()
data=inf.variables['tl1'][:,:,:]
data_range=[0.,1.]
cbar='Dark2'
nlevels=11
cbar_title='ch4 fraction'
plot_title='Fraction of Methane from Wetlands'

data_plot=data[0,0,:]
lats_plot=lats[:]
lons_plot=lons[:]


cm=plt.cm.get_cmap('RdYlBu_r')


fig=plt.figure(1)
ax=fig.add_subplot(111)
M=Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=-180,urcrnrlat=90,urcrnrlon=180,resolution='i')

(x,y)=M(lons_plot,lats_plot)
circles=[]
for x1,x2 in zip(x,y):
    circles.append(Circle((x1,x2),0.1))
   # squares.append(Rectangle((x1,x2),0.1,0.1))


##M.scatter(x,y,marker='o', cmap=cm, c=data_plot, vmin=240, vmax=320 , s=40, *c*='face')

# Create patch collection
patchcoll=PatchCollection(circles,cmap=cm,alpha=0.4)
patchcoll.set_array(np.array(data_plot))
patchcoll.set_edgecolor('face')

M.drawmapboundary()
ax.add_collection(patchcoll)
M.drawcountries(zorder=0)
M.drawcoastlines()

plt.colorbar()

plt.show()





