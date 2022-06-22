#!/bin/env python

import netCDF4 as nc
import numpy as np
import netcdftime as nctime
import PlotTools.plot_tools as PT
import os,sys
import matplotlib.pyplot as plt

CHESS_dir = '/users/eow/edwcom/GREENHOUSE/CHESS_runs/'

jef_file = CHESS_dir+'jules_output/Jvn4.3-E-F_repara.monthly_mean.nc'
#jules_file=CHESS_dir+'jules_output/chess_4.2_trif.annual.nc'
jules_file=CHESS_dir+'jules_output/Jvn4.3.monthly_mean.nc'

plot_dir=CHESS_dir+'plots/'

chess_grid_file = '/users/eow/edwcom/CHESS/chess_jules_land_index.nc'

params_to_plot = ['GPP','NEE']

Carbon_ConFact = (86400*365) # * (44/12)  # *1e6
C_units = '$kg$C $yr^{-1}'

# Read in grid file
grinf = nc.Dataset(chess_grid_file,'r')
grindex=grinf.variables['index_2D'][:]
lats_2D=grinf.variables['lats_2D'][:]
lons_2D=grinf.variables['lons_2D'][:]
grinf.close()
grimask=grindex.mask
fill_value=-999.

#Read in J-E-F data
jinf=nc.Dataset(jef_file,'r')
# mean data of time vector then multiply by conversion factor to give total CO2 for the final year of simulation (2013)
JEF_gpp=np.mean(jinf.variables['gpp_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
JEF_npp=np.mean(jinf.variables['npp_nuptake_out_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
JEF_resp_s=np.mean(jinf.variables['co2_soil_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
jinf.close()

# Put data on to 2D grid
JEF_gpp_2D = np.ma.masked_array(JEF_gpp[grindex],mask=grimask,fill_value=fill_value)
JEF_npp_2D = np.ma.masked_array(JEF_npp[grindex],mask=grimask,fill_value=fill_value)
JEF_resp_s_2D = np.ma.masked_array(JEF_resp_s[grindex],mask=grimask,fill_value=fill_value)
JEF_nee_2D = JEF_npp_2D-JEF_resp_s_2D
JEF_ter_2D = JEF_gpp_2D-JEF_nee_2D

#Read in JULES data
jinf=nc.Dataset(jules_file,'r')
# mean data of time vector then multiply by conversion factor to give total CO2 for the final year of simulation (2013)
#J_gpp=jinf.variables['gpp_gb'][-1:,:].squeeze()*Carbon_ConFact
#J_npp=jinf.variables['npp_gb'][-1:,:].squeeze()*Carbon_ConFact
#J_resp_s=jinf.variables['resp_s'][-1:,0,:].squeeze()*Carbon_ConFact
J_gpp=np.mean(jinf.variables['gpp_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
J_npp=np.mean(jinf.variables['npp_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
J_resp_s=np.mean(jinf.variables['resp_s_gb'][-12:,:].squeeze(),axis=0)*Carbon_ConFact
jinf.close()

# Put data on to 2D grid
J_gpp_2D = np.ma.masked_array(J_gpp[grindex],mask=grimask,fill_value=fill_value)
J_npp_2D = np.ma.masked_array(J_npp[grindex],mask=grimask,fill_value=fill_value)
J_resp_s_2D = np.ma.masked_array(J_resp_s[grindex],mask=grimask,fill_value=fill_value)
J_nee_2D = J_npp_2D-J_resp_s_2D
J_ter_2D = J_gpp_2D-J_nee_2D


FIG=plt.figure(figsize=(30,10))
AX=FIG.add_subplot(1,3,1)
PT.plot_map(J_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,.3], \
            COLOURS=['brown','white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JULES', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )

AX=FIG.add_subplot(1,3,2)
PT.plot_map(JEF_nee_2D-J_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-.1,.1], \
            COLOURS=['blue','white','red'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF - JULES', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )    

AX.text(0.03,-0.2,'JULES',ha='left',fontsize=17,transform=AX.transAxes)
AX.text(0.97,-0.2,'J-E-F',ha='right',fontsize=17,transform=AX.transAxes)

AX=FIG.add_subplot(1,3,3)
PT.plot_map(JEF_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,.3], \
            COLOURS=['brown','white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )


FIG.suptitle('Net Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'NetCarbonUptake2013_maps.png',bbox_inches='tight')
plt.show()

quit()


FIG=plt.figure(figsize=(30,10))
AX=FIG.add_subplot(1,3,1)
PT.plot_map(J_gpp_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,2],                       \
            COLOURS=['white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JULES', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )

AX=FIG.add_subplot(1,3,2)
PT.plot_map(JEF_gpp_2D-J_gpp_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.5,0.5],                       \
            COLOURS=['blue','white','red'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF - JULES', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )    

AX.text(0.03,-0.2,'JULES',ha='left',fontsize=17,transform=AX.transAxes)
AX.text(0.97,-0.2,'J-E-F',ha='right',fontsize=17,transform=AX.transAxes)

AX=FIG.add_subplot(1,3,3)
PT.plot_map(JEF_gpp_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,2],                       \
            COLOURS=['white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )


FIG.suptitle('Gross Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'GrossCarbonUptake2013_maps.png',bbox_inches='tight')
plt.show()




FIG=plt.figure(figsize=(10,10))
PT.plot_map(J_gpp_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,2],                       \
            COLOURS=['white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            FIGURE=FIG, )

FIG.suptitle('JULES Gross Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'JULES_GrossCarbonUptake2013_map.png',bbox_inches='tight')
plt.show()

FIG=plt.figure(figsize=(10,10))
PT.plot_map(JEF_gpp_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,2],                       \
            COLOURS=['white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            FIGURE=FIG, )

FIG.suptitle('J-E-F Gross Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'JEF_GrossCarbonUptake2013_map.png',bbox_inches='tight')
plt.show()





FIG=plt.figure(figsize=(10,10))
PT.plot_map(J_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,0.3],                       \
            COLOURS=['brown','white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            FIGURE=FIG, )

FIG.suptitle('JULES Net Ecosystem Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'JULES_NetEcoCarbonUptake2013_map.png',bbox_inches='tight')
plt.show()

FIG=plt.figure(figsize=(10,10))
PT.plot_map(JEF_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,0.3],                       \
            COLOURS=['white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            FIGURE=FIG, )

FIG.suptitle('J-E-F Net Ecosystem Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'JEF_NetEcoCarbonUptake2013_map.png',bbox_inches='tight')
plt.show()





FIG=plt.figure(figsize=(30,10))
AX=FIG.add_subplot(1,3,1)
PT.plot_map(J_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,.3], \
            COLOURS=['brown','white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JULES', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )

AX=FIG.add_subplot(1,3,2)
PT.plot_map(JEF_nee_2D-J_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-.1,.1], \
            COLOURS=['blue','white','red'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF - JULES', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )    

AX.text(0.03,-0.2,'JULES',ha='left',fontsize=17,transform=AX.transAxes)
AX.text(0.97,-0.2,'J-E-F',ha='right',fontsize=17,transform=AX.transAxes)

AX=FIG.add_subplot(1,3,3)
PT.plot_map(JEF_nee_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-0.1,.3], \
            COLOURS=['brown','white','khaki','gold','forestgreen'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=9,  \
            CBAR_LABEL='Carbon Uptake ('+C_units+')', \
            PLOT_TITLE='JEF', \
            TICK_FORMAT='%0.2f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )


FIG.suptitle('Net Carbon Uptake, 2013',fontsize=30)
FIG.savefig(plot_dir+'NetCarbonUptake2013_maps.png',bbox_inches='tight')
plt.show()



FIG=plt.figure(figsize=(30,10))
AX=FIG.add_subplot(1,3,1)
PT.plot_map(J_ter_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,3],\
            COLOURS=['white','brown'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Respired Carbon ('+C_units+')', \
            PLOT_TITLE='JULES', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )

AX=FIG.add_subplot(1,3,2)
PT.plot_map(JEF_ter_2D-J_ter_2D,lons_2D,lats_2D, \
            DATA_RANGE=[-1.5,1.5],\
            COLOURS=['blue','white','red'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Respired Carbon ('+C_units+')', \
            PLOT_TITLE='JEF - JULES', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )    

AX.text(0.03,-0.2,'JULES',ha='left',fontsize=17,transform=AX.transAxes)
AX.text(0.97,-0.2,'J-E-F',ha='right',fontsize=17,transform=AX.transAxes)

AX=FIG.add_subplot(1,3,3)
PT.plot_map(JEF_ter_2D,lons_2D,lats_2D, \
            DATA_RANGE=[0,3],\
            COLOURS=['white','brown'], \
            INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,  \
            CBAR_LABEL='Respired Carbon ('+C_units+')', \
            PLOT_TITLE='JEF', \
            TICK_FORMAT='%0.1f',               \
            LONDEL=2,LATDEL=2, \
            LON_RANGE=[-8,2],LAT_RANGE=[50,60],         \
            RESOLUTION='i',FONTSIZES=[12,12,15,18],      \
            AXIS=AX, )


FIG.suptitle('Total Ecosystem Respiration, 2013',fontsize=30)
FIG.savefig(plot_dir+'TotEcoResp2013_maps.png',bbox_inches='tight')
plt.show()





quit()

plt.subplot(1,3,1)
plt.imshow(J_resp_s_2D,origin='bottom',vmax=1.5)
plt.title('JULES')
plt.colorbar()

plt.subplot(1,3,2)
plt.imshow(JEF_resp_s_2D,origin='bottom',vmax=1.5)
plt.title('JEF')
plt.colorbar()

plt.subplot(1,3,3)
plt.imshow(J_resp_s_2D-JEF_resp_s_2D,origin='bottom',vmax=1,vmin=-1)
plt.title('J - JEF')
plt.colorbar()

plt.show()

