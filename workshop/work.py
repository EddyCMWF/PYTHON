#!/usr/bin/python

from __future__ import print_function, division
import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


weather=pd.read_table("http://bioltfws1.york.ac.uk/cgi-bin/weather",na_values=["n/a"])

weather['Pressure (hPa)']=weather['Pressure (hPa)']*100.
weather.columns.values[9]='Pressure (Pa)'
weather.columns=weather.columns.values


m = Basemap(llcrnrlon=-12.5,llcrnrlat=48.5,urcrnrlon=6.0,urcrnrlat=61.0, \
            resolution='f',projection='cass',lon_0=-4.36,lat_0=54.7)

m.drawcoastlines(zorder=2)
m.fillcontinents(color='green',lake_color='blue',zorder=1)
m.drawparallels(np.arange(-40,61.,2.),zorder=3, color='white')
m.drawmeridians(np.arange(-20.,21.,2.),zorder=3,color='white')

m.drawmapboundary(fill_color='blue')
m.drawcountries()


lons=list(weather['Longitude'])
lats=list(weather['Latitude'])
(x,y) = m(lons,lats)
m.scatter(x,y,marker='*',color='white',s=40,zorder=5)


plt.show()



#index = weather['Temperature (C)'] >= 10
#print(weather[index])
#weather[weather['Country']=='Wales']
#x = np.arange(10)


