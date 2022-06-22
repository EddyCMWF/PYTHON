#!/usr/bin/env python
#
# 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

iDISPLAY='N'

K_DATA_DIR   ='/prj/GREENHOUSE/SC_simulator/JULES_output/WFDEI_KovenMethod/'
K_FILE_TAG_1 ='JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_KovenMethod.dump.spin'
K_FILE_TAG_2 ='JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_KovenMethod_spin2.dump.spin'

SC_DATA_DIR ='/prj/GREENHOUSE/SC_simulator/OPERATIONAL_output/work_WFDEI_GLOBAL_20yr/'
SC_FILE_TAG ='J4.3_WFDEI_GLOBAL_20yr_spin'

OUTPUT_DIR = '/users/eow/edwcom/SC_simulator/plots/'

# soil Params etc.
nPOOLs=4
fill_value=-999.
CS_units = '$kg$C $m^{-2}$'
CS_max = 15 

# Read lat and lons from first Koven Dump file
print K_DATA_DIR+K_FILE_TAG_1+'1.19800101.0.nc'
inf=nc.Dataset(K_DATA_DIR+K_FILE_TAG_1+'1.19800101.0.nc','r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]
inf.close()

# Define indexes
index_names = ['Global','Northern Lats','Tropical','Southern Lats']
indexes     = [ lats>-9999., lats>30 , (lats<15)&(lats>-15), lats<-30 ]
colours     = [ 'red' , 'blue', 'orange', 'darkgreen' ]

#Read in Koven Data
# Koven parameters
Syear    = 1980
Eyear    = 1990
SpinRange = Eyear-Syear

nSPINs_1 = 7
K_years_1  = np.arange(0,SpinRange*nSPINs_1)

nSPINs_2 = 3
K_years_2  = np.arange(0,SpinRange*nSPINs_2)+K_years_1[-1]


K_CS_1 = {}
K_CS_2 = {}
for zone in index_names:
    K_CS_1[zone]= [] 
    K_CS_2[zone]= [] 

print 'KovenSpin-1'
for iSPIN in range(nSPINs_1):
    for year in range(Syear,Eyear):
        file=K_DATA_DIR+K_FILE_TAG_1+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        inf=nc.Dataset(file,'r')
        cs = inf.variables['cs'][:]
        cs = np.sum(cs,axis=0)
        for zone,index in zip(index_names,indexes):
            K_CS_1[zone].append(np.mean(cs[index]))
        inf.close()

print 'KovenSpin-2'
for iSPIN in range(nSPINs_2):
    for year in range(Syear,Eyear):
        file=K_DATA_DIR+K_FILE_TAG_2+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        inf=nc.Dataset(file,'r')
        cs = inf.variables['cs'][:]
        cs = np.sum(cs,axis=0)
        for zone,index in zip(index_names,indexes):
            K_CS_2[zone].append(np.mean(cs[index]))
        inf.close()


# Read in SC simulator data
# SC params
SC_Syear = 1980
SC_Eyear = 2000
SC_range = SC_Eyear-SC_Syear
SC_years_1=np.arange(0,SC_range)
SC_years_2=SC_years_1+SC_range
SC_years_3=SC_years_2+SC_range


SC_CS_1 = {}
SC_CS_2 = {}
SC_CS_3 = {}

# SC spin 1
file=SC_DATA_DIR+SC_FILE_TAG+'1.monthly_mean.nc'
print file
inf=nc.Dataset(file,'r')
cs=inf.variables['cs'][::12].squeeze()
cs=cs.sum(axis=1)
for zone,index in zip(index_names,indexes):
    SC_CS_1[zone]=cs[:,index].mean(axis=1)
inf.close()


# SC spin 2
file=SC_DATA_DIR+SC_FILE_TAG+'2.monthly_mean.nc'
print file
inf=nc.Dataset(file,'r')
cs=inf.variables['cs'][::12].squeeze()
cs=cs.sum(axis=1)
for zone,index in zip(index_names,indexes):
    SC_CS_2[zone]=cs[:,index].mean(axis=1)
inf.close()

# SC spin 3
file=SC_DATA_DIR+SC_FILE_TAG+'3.monthly_mean.nc'
print file
inf=nc.Dataset(file,'r')
cs=inf.variables['cs'][::12].squeeze()
cs=cs.sum(axis=1)
for zone,index in zip(index_names,indexes):
    SC_CS_3[zone]=cs[:,index].mean(axis=1)
inf.close()





# Create figure to plot to
FIG = plt.figure(figsize=[18,12])

AX  = FIG.add_subplot(1,1,1)

for iZON in range(len(indexes)):
    index=indexes[iZON]
    zone=index_names[iZON]
    # print index.shape
    K1_CS_timeseries = np.array(K_CS_1[zone])
    K2_CS_timeseries = np.array(K_CS_2[zone])
    AX.plot(K_years_1,K1_CS_timeseries,color=colours[iZON],lw=2,ls=':')
    AX.plot(K_years_2,K2_CS_timeseries,lw=2,color=colours[iZON],ls=':')

    SC1_timeseries = np.array(SC_CS_1[zone])
    SC2_timeseries = np.array(SC_CS_2[zone])
    SC3_timeseries = np.array(SC_CS_3[zone])
    AX.plot(SC_years_1,SC1_timeseries,color=colours[iZON],lw=2,label=zone)
    AX.plot(SC_years_2,SC2_timeseries,lw=2,color=colours[iZON])
    AX.plot(SC_years_3,SC3_timeseries,lw=2,color=colours[iZON])
        

AX.plot(np.array([K_years_1[-1],K_years_1[-1]]),np.array([0,CS_max]),lw=2,ls='-',color='grey')
AX.plot(np.array([SC_years_1[-1],SC_years_1[-1]]),np.array([0,CS_max]),lw=2,ls='-',color='black')
AX.plot(np.array([SC_years_2[-1],SC_years_2[-1]]),np.array([0,CS_max]),lw=2,ls='-',color='black')
AX.set_ybound(lower=0,upper=CS_max)

AX.legend( ncol=4 )
FIG.text(0.03,0.6,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])

if (iDISPLAY=='Y'):
    plt.show()
else:
    FIG.savefig(OUTPUT_DIR+'WFDEI_SCversusKOVEN_Timeseries.png')
    plt.close()

   

    
