#!/usr/bin/env python

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as col
import numpy as np
import csv

PALS_DIR = '/prj/GREENHOUSE/PALS_comparison/'
site_metadata_file=PALS_DIR+'sites_meta_forlocplot.csv'
site_md_list = list(csv.reader(open(site_metadata_file,'r'),delimiter=','))
site_md_headers=site_md_list.pop(0)

site_names= [ site_md[0].strip()[:3] for site_md in site_md_list]
site_lats = np.array([ float(site_md[7]) for site_md in site_md_list])
site_lons = np.array([ float(site_md[8]) for site_md in site_md_list])
nSITEs=len(site_names)

projection='cyl'
res='c'
latdel=30
londel=30
latrange=[-60,90]
lonrange=[-180,180]

FIG,AX = plt.subplots(nrows=1,ncols=1,figsize=[22,12])
FIG.subplots_adjust(top=0.9)
M = Basemap(projection=projection, \
            resolution=res, \
            llcrnrlat =latrange[0],urcrnrlat=latrange[1], \
            llcrnrlon =lonrange[0],urcrnrlon=lonrange[1]  )
            
site_x,site_y=M(site_lons,site_lats)

M.drawcountries()
M.drawcoastlines(linewidth=1.)
M.drawmapboundary(fill_color='powderblue')
M.fillcontinents(color='darkkhaki',lake_color='powderblue',zorder=0)
parallels = np.arange(-90.,90.,latdel)
M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 1.,fontsize=20)
meridians = np.arange(-180.,180.,londel)
M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 1.,fontsize=20)
M.scatter(site_x,site_y,marker='*',color='red',s=60) #,lw=0)
for iSITE in range(nSITEs):
    AX.text(site_x[iSITE]-5,site_y[iSITE]+0.5,site_names[iSITE],color='k',fontsize=15,weight='bold')

FIG.suptitle('PALS site locations', fontsize=50, weight='bold')
FIG.savefig(PALS_DIR+'site_locations.png',bbox_inches='tight')

plt.show()


