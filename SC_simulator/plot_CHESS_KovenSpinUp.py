#!/usr/bin/env python
#
# 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

iDISPLAY='N'

DATA_DIR   ='/prj/GREENHOUSE/SC_simulator/JULES_output/KovenMethod/'
FILE_TAG_1 ='JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_KovenMethod.dump.spin'
FILE_TAG_2 ='JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_KovenMethod_spin2.dump.spin'
OUTPUT_DIR = '/users/eow/edwcom/SC_simulator/plots/'

Syear    = 1980
Eyear    = 1990
SpinRange = Eyear-Syear

nSPINs_1 = 7
years_1  = np.arange(0,SpinRange*nSPINs_1)

nSPINs_2 = 3
years_2  = np.arange(0,SpinRange*nSPINs_2)+years_1[-1]

nPOOLs=4
fill_value=-999.
CS_pools = ['DPM', 'RPM', 'BIO', 'HUM']
CS_units = '$kg$C $m^{-2}$'
CS_maxes = [0.1,5,0.4,12]

print DATA_DIR+FILE_TAG_1+'1.19800101.0.nc'
inf=nc.Dataset(DATA_DIR+FILE_TAG_1+'1.19800101.0.nc','r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]
inf.close()

index_names = ['Global','Northern Lats','Tropical','Southern Lats']
indexes     = [ lats>-9999., lats>30 , (lats<15)&(lats>-15), lats<-30 ]
colours     = [ 'red' , 'blue', 'orange', 'darkgreen' ]

CS_1 = {}
CS_2 = {}
for zone in index_names:
    CS_1[zone]=[ [] for i in range(nPOOLs)]
    CS_2[zone]=[ [] for i in range(nPOOLs)]

print 'KovenSpin-1'
for iSPIN in range(nSPINs_1):
    #print iSPIN 
    for year in range(Syear,Eyear):
        #print str(year)
        file=DATA_DIR+FILE_TAG_1+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        inf=nc.Dataset(file,'r')
        cs = inf.variables['cs'][:]
        for zone,index in zip(index_names,indexes):
            for iPOOL in range(nPOOLs):
                CS_1[zone][iPOOL].append(np.mean(cs[iPOOL,index]))

        inf.close()

print 'KovenSpin-2'
for iSPIN in range(nSPINs_2):
    #print iSPIN 
    for year in range(Syear,Eyear):
        #print str(year)
        file=DATA_DIR+FILE_TAG_2+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        inf=nc.Dataset(file,'r')
        cs = inf.variables['cs'][:]
        for zone,index in zip(index_names,indexes):
            for iPOOL in range(nPOOLs):
                CS_2[zone][iPOOL].append(np.mean(cs[iPOOL,index]))

        inf.close()


# Create figure to plot to
FIG = plt.figure(figsize=[18,12])
nplts_wdth   = 2.
nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)


for iPOOL in range(nPOOLs):
    print 'Plotting: '+CS_pools[iPOOL]
    AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
    AX.set_title(CS_pools[iPOOL])

    for iZON in range(len(indexes)):
        index=indexes[iZON]
        zone=index_names[iZON]
        # print index.shape
        K1_CS_timeseries = np.array(CS_1[zone][iPOOL])
        K2_CS_timeseries = np.array(CS_2[zone][iPOOL])
        print years_1.shape,K1_CS_timeseries.shape
        AX.plot(years_1,K1_CS_timeseries,label=zone,color=colours[iZON],lw=2)
        print years_2.shape,K2_CS_timeseries.shape
        AX.plot(years_2,K2_CS_timeseries,label=zone,lw=2)
        
    AX.plot(np.array([years_1[-1],years_1[-1]]),np.array([0,12]),lw=2,ls=':',color='black')
    AX.set_ybound(lower=0,upper=12)

AX.legend( bbox_to_anchor=(-0.1,-0.11),loc=10,borderaxespad=0.,ncol=4 )
FIG.text(0.03,0.6,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])
FIG.suptitle('Standard Soil Carbon Spin-Up Time-series',fontsize=24)

if (iDISPLAY=='Y'):
    plt.show()
else:
    FIG.savefig(OUTPUT_DIR+'KoveSpinUp_Timeseries.png')
    plt.close()

   

    
