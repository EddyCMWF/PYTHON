#!/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import sys,os
from PlotTools import plot_tools as PTs

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

nland=1631

gcm=optional_argparse('-gcm','CEN_MOHC_MOD_HadGEM2-ES')
climat_dir=optional_argparse('-climat_dir','/prj/CLIFFTOP/COMMON_DATA/WATCH/climatol/')
if climat_dir[-1:]!='/': climat_dir+='/'
pattern_dir=optional_argparse('-pattern_dir','/prj/CLIFFTOP/COMMON_DATA/CMIP/cmip5_patterns/')
if pattern_dir[-1:]!='/': pattern_dir+='/'
out_dir=optional_argparse('-out_dir','/users/eow/edwcom/CLIFFTOP/plots/DRIVE_DATA/')
if out_dir[-1:]!='/': out_dir+='/'

#
drive_months=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
# Climatology variables in file order:
climat_vars=['LON','LAT','T','RH15M','UWIND','VWIND','LW',
             'SW','DTEMP','RAINFALL','SNOWFALL','PSTAR_HA','F_WET']

# Read in climat data
CLIM_DATA_DICT={ var:np.zeros([12,nland]) for var in climat_vars }
for i_mnth in range(12):
    clim_file=climat_dir+drive_months[i_mnth]
    print(clim_file)
    inf=open(clim_file,'r')
    inlines=inf.readlines()
    inlines.pop(0)
    inf.close()
    for i_pt in range(nland):
        line=inlines[i_pt].replace('\n','')
        if len(line)>1:
            split=line.split()
            for var,val in zip(climat_vars,split):
                CLIM_DATA_DICT[var][i_mnth,i_pt]=float(val)

# Read in climat data
ANOM_DATA_DICT={ var:np.zeros([12,nland]) for var in climat_vars }
for i_mnth in range(12):
    anom_file=pattern_dir+gcm+'/'+drive_months[i_mnth]
    print(anom_file)
    inf=open(anom_file,'r')
    inlines=inf.readlines()
    inlines.pop(0)
    inf.close()
    for i_pt in range(nland):
        line=inlines[i_pt].replace('\n','')
        if len(line)>1:
            split=line.split()
            for var,val in zip(climat_vars,split):
                ANOM_DATA_DICT[var][i_mnth,i_pt]=float(val)

climat_vars=climat_vars[2:]

# build 2D imogen array:
minlat,maxlat,reslat=-55,90,2.5
minlon,maxlon,reslon=0,360,3.75

# 2D arrays:
lats=np.arange(minlat,maxlat+0.1,reslat)
lons=np.arange(minlon,maxlon+0.1,reslon)
lons_2D,lats_2D=np.meshgrid(lons,lats)
index_2D = np.zeros_like(lons_2D,dtype=int)-1

# get in lat/lon data and convert to index arrays
inlats=CLIM_DATA_DICT['LAT'][0,:]
latdex=(inlats-minlat)/reslat
latdex=latdex.astype('int')

inlons=CLIM_DATA_DICT['LON'][0,:]
londex=(inlons/reslon)-minlon
londex=londex.astype('int')

for i_pt in range(nland):
    index_2D[latdex[i_pt],londex[i_pt]]=int(i_pt)

index_2D=np.ma.masked_equal(index_2D,-1)


CLIMAT_VAR_CBAR_colours=[ [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          [ 'teal','white','orangered'],
                          ]

CLIMAT_VAR_limits = [ [ 240,320 ] ,
                      [ 0,100 ],
                      [ 0,5 ],
                      [ 0,10 ],
                      [ 100,500 ],
                      [ 100,400 ],
                      [ 0,20 ],
                      [ 0,20 ],
                      [ 800,1200 ],
                      [ 0,0.5 ],
                          ]

ANOM_VAR_limits = [ [ 0,4 ],
                    [-5,5 ],
                    [-0.5,0.5 ],
                    [-0.5,0.5 ],
                    [ 0,20 ],
                    [ -10,10 ],
                    [ -2,2 ],
                    [ -2,2 ],
                    [ -2,2 ],
                    [ -1.5,1.5 ],
                    [ -0.1,0.1 ],
                   ]

# plot the clim and anom data
for ivar in range(len(climat_vars)):
    var=climat_vars[ivar]
    CBAR_colours=CLIMAT_VAR_CBAR_colours[ivar]
    CLIM_RANGE=CLIMAT_VAR_limits[ivar]
    ANOM_RANGE=ANOM_VAR_limits[ivar]
    if (not os.path.isfile(out_dir+'climat/'+var+'.gif'))&(var not in ['VWIND','SNOWFALL']):
        print('Plotting Climatologies')
        for i_mnth in range(12):
            if (not os.path.isfile(out_dir+'climat/'+var+'_'+str('%02d'%(i_mnth+1,))+'.png')):
                plotdata=np.ma.masked_array(CLIM_DATA_DICT[var][i_mnth,index_2D],mask=index_2D.mask)
                PLOT_TITLE=str(var)+' - '+drive_months[i_mnth]
                print(PLOT_TITLE)
                PTs.plot_map(plotdata,lons_2D,lats_2D,RESOLUTION='c',
                             DATA_RANGE=CLIM_RANGE,
                             COLOURS=CBAR_colours,INTERPOLATE_COLOURS=True,
                             NLEVELS=250,NTICKS=11, 
                             PLOT_TITLE= PLOT_TITLE,
                             FILE_PLOT=out_dir+'climat/'+var+'_'+str('%02d'%(i_mnth+1,))+'.png' ) 
                plt.close()
        os.system('convert -delay 100 -loop 0 '+out_dir+'climat/'+var+'*.png '+out_dir+'climat/'+var+'.gif')

    if (not os.path.isfile(out_dir+'anom/'+var+'.gif'))&(var not in ['VWIND','SNOWFALL','F_WET']):
        print('Plotting Anomalies')
        for i_mnth in range(12):
            if (not os.path.isfile(out_dir+'anom/'+var+'_'+str('%02d'%(i_mnth+1,))+'.png')):
                plotdata=np.ma.masked_array(ANOM_DATA_DICT[var][i_mnth,index_2D],mask=index_2D.mask)
                PLOT_TITLE=str(var)+' - '+drive_months[i_mnth]
                print(PLOT_TITLE)
                PTs.plot_map(plotdata,lons_2D,lats_2D,RESOLUTION='c',
                             DATA_RANGE=ANOM_RANGE,
                             COLOURS=CBAR_colours,INTERPOLATE_COLOURS=True,
                             NLEVELS=250,NTICKS=11, 
                             PLOT_TITLE= PLOT_TITLE,
                             FILE_PLOT=out_dir+'anom/'+var+'_'+str('%02d'%(i_mnth+1,))+'.png' ) 
                plt.close()
        os.system('convert -delay 100 -loop 0 '+out_dir+'anom/'+var+'*.png '+out_dir+'anom/'+var+'.gif')


