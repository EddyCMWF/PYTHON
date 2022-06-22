#!/usr/bin/python2.7
import os, glob
import netCDF4 as nc
import numpy as np
import netcdftime as nctime
import csv
import matplotlib.pyplot as plt
#
#####################################################################################
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

GPP_diurnal=np.array(J_hourly_dict['# gpp'])[:n_h_points].reshape(n_d_points,24)*86.4e6
fsmc_diurnal=np.array(J_hourly_dict['fsmc'])[:n_h_points].reshape(n_d_points,24)

airt_diurnal  = np.array(met_dict['airt'])[:n_h_points].reshape(n_d_points,24)
co2_diurnal   = np.array(met_dict['co2'])[:n_h_points].reshape(n_d_points,24)
wind_diurnal  = np.array(met_dict['wind'])[:n_h_points].reshape(n_d_points,24)
SWRAD_diurnal = np.array(met_dict['SWRAD'])[:n_h_points].reshape(n_d_points,24)
LWRAD_diurnal = np.array(met_dict['LWRAD'])[:n_h_points].reshape(n_d_points,24)
Psurf_diurnal = np.array(met_dict['sf_pressure'])[:n_h_points].reshape(n_d_points,24)
SH_diurnal    = np.array(met_dict['sp_humidity'])[:n_h_points].reshape(n_d_points,24)

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

BL_index=np.where(PFT==0)[0]
BL_hLAI_index=np.where( (PFT==0) & (LAI>6) & (LAI<6.5) )[0]

latlons_BLh=zip(lat[BL_hLAI_index],lon[BL_hLAI_index])
ulatlons_BLh=list(set(latlons_BLh))

site_index=[]
for i,latlon in zip(range(len(lat)),zip(lat,lon)):
    if latlon in ulatlons_BLh:
        site_index.append(i)

site_index=np.array(site_index)

#C3_index=np.where(PFT==2)[0]
#C3_hLAI_index=np.where( (PFT==2) & (LAI>5) )[0]

y_data=GPP
x_data=LAI


#plt.plot(x_data,y_data,ls='',marker='.')
#plt.plot(x_data[site_index],y_data[site_index],ls='',marker='.',color='red')
#plt.xlabel('LAI')
#plt.ylabel('GPP')
#plt.show()


plt.plot(x_data[BL_index],y_data[BL_index],ls='',marker='.',color='grey')

col_cnt=0
colours=['blue','red','green','yellow','cyan','darkorange','yellowgreen']
for site in ulatlons_BLh:
    index=np.where( (lat==site[0]) & (lon==site[1]) )[0]
    plt.plot(x_data[index],y_data[index],ls='',marker='.',label=str(site),c=colours[col_cnt])
    col_cnt+=1

plt.xlabel('LAI')
plt.ylabel('GPP')
plt.legend()
plt.show()




plt.subplot(2,4,1)
ydata=SWRAD_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.xlabel('hour of day')
plt.ylabel('SWRAD')


plt.subplot(2,4,2)
ydata=LWRAD_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.xlabel('hour of day')
plt.ylabel('LWRAD')


plt.subplot(2,4,3)
ydata=airt_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.xlabel('hour of day')
plt.ylabel('airt')


plt.subplot(2,4,4)
ydata=co2_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.xlabel('hour of day')
plt.ylabel('co2')


plt.subplot(2,4,5)
ydata=wind_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.ylabel('wind')
plt.xlabel('hour of day')


plt.subplot(2,4,6)
ydata=Psurf_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.ylabel('Psurf')
plt.xlabel('hour of day')


plt.subplot(2,4,7)
ydata=SH_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.ylabel('SH')
plt.xlabel('hour of day')


plt.subplot(2,4,8)
ydata=GPP_diurnal
for day in BL_index:
    plt.plot(ydata[day,:],ls='-',color='grey')

for day in BL_hLAI_index:
    plt.plot(ydata[day,:],ls='-')

plt.ylabel('GPP')
plt.xlabel('hour of day')


#plt.subplot(2,4,8)
#ydata=fsmc_diurnal
#for day in BL_index:
#    plt.plot(ydata[day,:],ls='-',color='grey')
#col_cnt=0
#for site in ulatlons_BLh:
#    index=np.where( (lat==site[0]) & (lon==site[1]) )[0]
#    for day in index:
#        plt.plot(ydata[day,:],ls='-',c=colours[col_cnt])
#    col_cnt+=1
#
#plt.xlabel('hour of day')
#plt.ylabel('fsmc')

plt.show()








