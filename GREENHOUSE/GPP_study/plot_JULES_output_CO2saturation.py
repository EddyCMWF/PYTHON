#!/usr/bin/python2.7
import os
import netCDF4 as nc
import numpy as np
import netcdftime as nctime
import csv, glob
import matplotlib.pyplot as plt
import matplotlib as mpl
#
#####################################################################################
#
GH_base_dir='/users/eow/edwcom/GREENHOUSE/GPP_sim_forLuke/'
runs=glob.glob(GH_base_dir+'data_*')
print 'Available runs: '
for i in range(len(runs)):
    print i,':',runs[i]
iRUN=input('Select a run: ')
run=runs[iRUN]

GH_dir     = run+'/'
pheno_file = GH_dir+'downscaled_pheno_drivers.csv'
met_file   = GH_dir+'downscaled_met_drivers.csv'

#
#canradmod=raw_input('Choose a canrad model: ')
JULES_dir     = GH_dir # '/users/eow/edwcom/GREENHOUSE/GPP_sim_forLuke/data_201606/'
J_hourly_file = JULES_dir+'JULES_hourly_output.dat' #.canrad'+canradmod
J_daily_file  = JULES_dir+'JULES_daily_output.dat' #.canrad'+canradmod
#
# read in pheno data
pheno_csv=open(pheno_file,'r')
data=list(csv.reader(pheno_csv,delimiter=','))
pheno_hdrs=data.pop(0)
pheno_dict={}
for hdr in pheno_hdrs:
    pheno_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(pheno_hdrs,dat_lin[1:]):
        pheno_dict[hdr].append(float(dat))
#
pheno_csv.close()
del data
#
# read in the met data    
met_csv=open(met_file,'r')
data=list(csv.reader(met_csv,delimiter=','))
met_hdrs=data.pop(0)
met_hdrs=[hdr.replace(' ','') for hdr in met_hdrs]
met_dict={}
for hdr in met_hdrs:
    met_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(met_hdrs,dat_lin[1:]):
        met_dict[hdr].append(float(dat))
#
met_csv.close()
del data
#
met_dict['co2'] = [ co2 * ((28.97/44.01)**2.) for co2 in met_dict['co2'] ]
#
# read in JULES daily data
J_daily_csv=open(J_daily_file,'r')
data=list(csv.reader(J_daily_csv,delimiter=','))
J_daily_hdrs=data.pop(0)
J_daily_hdrs= [ hdr.replace('#','').replace(' ','') for hdr in J_daily_hdrs ]
J_daily_dict={}
for hdr in J_daily_hdrs:
    J_daily_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(J_daily_hdrs,dat_lin):
        J_daily_dict[hdr].append(float(dat))
#
J_daily_csv.close()
del data
#
#
# read in JULES hourly data
J_hourly_csv=open(J_hourly_file,'r')
data=list(csv.reader(J_hourly_csv,delimiter=','))
J_hourly_hdrs=data.pop(0)
J_hourly_hdrs= [ hdr.replace('#','').replace(' ','') for hdr in J_hourly_hdrs ]
J_hourly_dict={}
for hdr in J_hourly_hdrs:
    J_hourly_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(J_hourly_hdrs,dat_lin):
        J_hourly_dict[hdr].append(float(dat))
#
J_hourly_csv.close()
del data
#
#

n_d_points=len(J_daily_dict[J_daily_hdrs[-1]])
n_h_points=n_d_points*24.
print n_d_points, n_h_points

LAI=np.array(J_daily_dict['LAI'])[:n_d_points]
GPP=np.array(J_daily_dict['gpp'])[:n_d_points]*86.4e6
NPP=np.array(J_daily_dict['npp'])[:n_d_points]*86.4e6
TSOIL=np.array(J_daily_dict['t_soil'])[:n_d_points]
SM=np.array(J_daily_dict['SM'])[:n_d_points]
FSMC=np.array(J_daily_dict['fsmc'])[:n_d_points]
PFT=np.array(J_daily_dict['PFT'])[:n_d_points]

GPP_diurnal=np.array(J_hourly_dict['gpp'])[:n_h_points].reshape(n_d_points,24)*86.4e6

fsmc_diurnal=np.array(J_hourly_dict['fsmc'])[:n_h_points].reshape(n_d_points,24)


airt_diurnal  = np.array(met_dict['airt']).reshape(n_d_points,24)[:n_d_points]
co2_diurnal   = np.array(met_dict['co2']).reshape(n_d_points,24)[:n_d_points]
wind_diurnal  = np.array(met_dict['wind']).reshape(n_d_points,24)[:n_d_points]
SWRAD_diurnal = np.array(met_dict['SWRAD']).reshape(n_d_points,24)[:n_d_points]
LWRAD_diurnal = np.array(met_dict['LWRAD']).reshape(n_d_points,24)[:n_d_points]
Psurf_diurnal = np.array(met_dict['sf_pressure']).reshape(n_d_points,24)[:n_d_points]
SH_diurnal    = np.array(met_dict['sp_humidity']).reshape(n_d_points,24)[:n_d_points]

airt_dmean = np.mean(airt_diurnal,axis=1)
airt_dmax  = np.max(airt_diurnal,axis=1)
co2_dmax   = np.max(co2_diurnal,axis=1)
wind_dmax  = np.max(wind_diurnal,axis=1)
SWRAD_dmax = np.max(SWRAD_diurnal,axis=1)
LWRAD_dmax = np.max(LWRAD_diurnal,axis=1)
Psurf_dmax = np.max(Psurf_diurnal,axis=1)
SH_dmax    = np.max(SH_diurnal,axis=1)

lat = np.array(pheno_dict['lat'])[:n_d_points]
lon = np.array(pheno_dict['long'])[:n_d_points]
avgN= np.array(pheno_dict['avgN'])[:n_d_points]
LAI_pheno=np.array(pheno_dict['LAI'])[:n_d_points]
avg_airt=np.array(pheno_dict['avg_airt'])[:n_d_points]



N_xpts = 5
x_sp   = 100.
x_ep   = 1000.
x_inc  = (x_ep-x_sp)/N_xpts
x_linpts = np.arange(x_sp,x_ep,x_inc,dtype=float)

fig=plt.figure(figsize=[12,18])
cmap=plt.cm.get_cmap('RdYlBu_r')
# plot CO2 vs GPP with temperature colour bar
ax=fig.add_subplot(2,1,1)
y_data=GPP
x_data=co2_dmax
z_data=airt_dmean
for i in range(10):
    y_linpts = np.zeros_like(x_linpts)
    for j in range(N_xpts):
        jindex=np.where( (z_data>((i*5)-15)) & (z_data<=(i*5)-10) &   \
                         (x_data>=x_linpts[j]) & (x_data<x_linpts[j]+x_inc) )
        if len(jindex[0])>0:
            y_linpts[j] = np.mean(y_data[jindex])
        else:
            y_linpts[j] = float(np.nan)
    col_num = int( (i*18.)+70  )
    ax.plot(x_linpts+(x_inc/2.),y_linpts,color=cmap(col_num),linewidth=2,\
            label=str((i*5)-15)+' - '+str((i*5)-10)+' $^o$C'  )
ax.set_xlabel('co2')
ax.set_xlim( [100,800] )
ax.set_ylabel('GPP')
ax.set_ylim( [0,10] )
ax.legend(loc='upper left',ncol=2)

# plot CO2 vs GPP with SWRAD colour bar
ax=fig.add_subplot(2,1,2)
z_data=SWRAD_dmax
for i in range(10):
    y_linpts = np.zeros_like(x_linpts)
    for j in range(N_xpts):
        jindex=np.where( (z_data>500+(i*100)) & (z_data<=600+(i*100))  &   \
                         (x_data>=x_linpts[j]) & (x_data<x_linpts[j]+x_inc) )
        if len(jindex[0])>0:
            y_linpts[j] = np.mean(y_data[jindex])
        else:
            y_linpts[j] = float(np.nan)
    col_num = int( (i*18.)+70  )
    ax.plot(x_linpts+(x_inc/2.),y_linpts,color=cmap(col_num),linewidth=2,\
            label=str(500+(i*100))+' - '+str(600+(i*100))+' (W $m^{-2}$ $s^{-1}$)' )
ax.set_xlabel('co2')
ax.set_xlim( [100,800] )
ax.set_ylabel('GPP')
ax.set_ylim( [0,10] )
ax.legend(loc='upper left',ncol=2)

plt.show()




fig=plt.figure(figsize=[12,18])

# plot CO2 vs GPP with temperature colour bar
ax=fig.add_subplot(2,1,1)
x_data=co2_dmax
x_linpts = np.arange(100,1001,10,dtype=float)
z_data=airt_dmean
cmap=plt.cm.get_cmap('RdYlBu_r')
colbounds = np.arange(-10,35,5)
norm=mpl.colors.BoundaryNorm(colbounds,cmap.N)
sc=ax.scatter(x_data,y_data,c=z_data,cmap=cmap,norm=norm,vmin=-10,vmax=35)
for i in range(10,0,-1):
    index=np.where( (z_data>((i*5)-15)) & (z_data<=(i*5)-10) )
    ax.scatter(x_data[index],y_data[index],c=z_data[index],cmap=cmap,norm=norm,vmin=-10,vmax=35)
ax.set_xlabel('co2')
ax.set_xlim( [100,800] )
ax.set_ylabel('GPP')
ax.set_ylim( [0,15] )
cbar=plt.colorbar(sc)
cbar.set_label('Mean Daily Temperature ($^o$C)')


# plot CO2 vs GPP with SWRAD colour bar
ax=fig.add_subplot(2,1,2)
x_data=co2_dmax
z_data=SWRAD_dmax
cmap=plt.cm.get_cmap('RdYlBu_r')
colbounds = np.arange(500,1500,100)
norm=mpl.colors.BoundaryNorm(colbounds,cmap.N)
sc=ax.scatter(x_data,y_data,c=z_data,cmap=cmap,norm=norm,vmin=500,vmax=1500)
for i in range(0,10,1):
    index=np.where( (z_data>500+(i*100)) & (z_data<=600+(i*100)) )
    ax.scatter(x_data[index],y_data[index],c=z_data[index],cmap=cmap,norm=norm,vmin=500,vmax=1500)
ax.set_xlabel('co2')
ax.set_xlim( [100,800] )
ax.set_ylabel('GPP')
ax.set_ylim( [0,15] )
cbar=plt.colorbar(sc)
cbar.set_label('Maximum Daily SWRAD (W $m^{-2}$ $s^{-1}$)')

plt.show()



