#!/usr/bin/python2.7
import os
import netCDF4 as nc
import numpy as np
import netcdftime as nctime
import csv
import matplotlib.pyplot as plt
import glob
from PlotTools import plot_tools as PT
#
#####################################################################################
GH_base_dir='/users/eow/edwcom/GREENHOUSE/GPP_sim_forLuke/'
runs=glob.glob(GH_base_dir+'data_*')
print 'Available runs: '
for i in range(len(runs)):
    print i,':',runs[i]
iRUN=input('Select a run: ')
run=runs[iRUN][len(GH_base_dir):]

GH_dir=GH_base_dir+run+'/'
if 'uk' in run:
    pheno_file = GH_dir+'downscaled_pheno_drivers_uk.csv'
    met_file   = GH_dir+'downscaled_met_drivers_uk.csv'
else:
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
#
# read in JULES daily data
J_daily_csv=open(J_daily_file,'r')
data=list(csv.reader(J_daily_csv,delimiter=','))
J_daily_hdrs=data.pop(0)
J_daily_hdrs=[ hdr.strip() for hdr in J_daily_hdrs ]
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
J_hourly_hdrs=[ hdr.strip() for hdr in J_hourly_hdrs ]
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
print np.array(met_dict['airt']).shape

LAI=np.array(J_daily_dict['LAI'])[:n_d_points]
GPP=np.array(J_daily_dict['# gpp'])[:n_d_points]*86.4e6
NPP=np.array(J_daily_dict['npp'])[:n_d_points]*86.4e6
TSOIL=np.array(J_daily_dict['t_soil'])[:n_d_points]
SM=np.array(J_daily_dict['SM'])[:n_d_points]
FSMC=np.array(J_daily_dict['fsmc'])[:n_d_points]
PFT=np.array(J_daily_dict['PFT'])[:n_d_points]

airt_dmax  = np.max(np.array(met_dict['airt'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
co2_dmax   = np.max(np.array(met_dict['co2'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
wind_dmax  = np.max(np.array(met_dict['wind'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
SWRAD_dmax = np.max(np.array(met_dict['SWRAD'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
LWRAD_dmax = np.max(np.array(met_dict['LWRAD'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
Psurf_dmax = np.max(np.array(met_dict['sf_pressure'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]
SH_dmax    = np.max(np.array(met_dict['sp_humidity'])[:n_h_points].reshape(n_d_points,24),axis=1)[:n_d_points]

lat = np.array(pheno_dict['lat'])[:n_d_points]
lon = np.array(pheno_dict['long'])[:n_d_points]
avgN= np.array(pheno_dict['avgN'])[:n_d_points]
LAI_pheno=np.array(pheno_dict['LAI'])[:n_d_points]
avg_airt=np.array(pheno_dict['avg_airt'])[:n_d_points]

y_data=GPP

xdatas=[ SWRAD_dmax,LWRAD_dmax,SH_dmax,airt_dmax,co2_dmax,avgN,SM,TSOIL,LAI,NPP,lat,lon]
x_labels=['SWRAD','LWRAD','Spec. Hum.','air temp.','co2','Leaf Nitrogen','Soil Moist.','T soil (JULES)','LAI','NPP','Latitude','Longitude']

FIG,AXES=plt.subplots(3,4,figsize=(20,12))
for i in range(12):
    x_data=xdatas[i]
    AX=AXES.flat[i]
    AX.plot(x_data,y_data,ls='',marker='.')
    AX.set_ylabel('GPP')
    AX.set_xlabel(x_labels[i])
    AX.set_ylim(0,30)
FIG.suptitle('GPP vs state variables - '+run,fontsize=30)
plt.show()

#Plot by PFT color
FIG,AXES=plt.subplots(3,4,figsize=(20,12))
#if not 'uk' in run:
PFTcols=['forestgreen','springgreen','yellow','orange','saddlebrown']
PFTnames=['Broadleaf','Needleleaf','C3 Grass','C4 Grass','Shrub']
#else:
#    PFTcols=['forestgreen','springgreen','yellow','saddlebrown','c']
#    PFTnames=['Broadleaf','Needleleaf','C3 Grass','Shrub','Crops']
for i in range(12):
    x_data=xdatas[i]
    AX=AXES.flat[i]
    for iPFT in [3,2,0,1,4]:
        index=PFT==iPFT
        AX.scatter(x_data[index],y_data[index],marker='.',lw=0, label=PFTnames[iPFT],c=PFTcols[iPFT] )
    AX.set_ylabel('GPP')
    AX.set_xlabel(x_labels[i])
    AX.set_ylim(0,30)
handles,labels=AX.get_legend_handles_labels()
FIG.legend(handles,labels,loc=8,ncol=5)
FIG.suptitle('GPP vs state variables by PFT - '+run, fontsize=30)
plt.show()


CMAP=PT.custom_div_cmap(numcolors=255,colors=['yellow','peru','orangered'])
#xdatas=xdatas[:10]
#x_labels=x_labels[:10]
FIG,AXES=plt.subplots(3,4,figsize=(20,15))
for i in range(12):
    x_data=xdatas[i]
    AX=AXES.flat[i]
    AX.hist2d(x_data,y_data,(100,100),range=([x_data.min(),x_data.max()],[0.5,20]),cmap=CMAP, cmin=0.1)
    AX.set_ylabel('GPP')
    AX.set_xlabel(x_labels[i])
    AX.set_ylim(0,20)
FIG.suptitle('GPP vs state variables - '+run, fontsize=30)
plt.show()


