#!/usr/bin/env python
#
# 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob,os

iDISPLAY='N'

K_DATA_DIR   ='/prj/GREENHOUSE/SC_simulator/JULES_output/CHESS_KovenMethod/'
K_FILE_TAG_1 ='JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod.dump.spin'
K_FILE_TAG_2 ='JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod_spin2.dump.spin'

SC_DATA_DIR  = '/prj/GREENHOUSE/SC_simulator/OPERATIONAL_output/work_CHESS_20yr/'
SC_FILE_TAG  = 'J4.3_CHESS_spin'
OUTPUT_DIR = '/users/eow/edwcom/SC_simulator/plots/'

nPOOLs=4
fill_value=-999.
CS_pools = ['DPM', 'RPM', 'BIO', 'HUM']
CS_units = '$kg$C $m^{-2}$'
CS_maxes = [0.1,5,0.4,12]


# Read in the Koven data from spin dumps and create time-series
Syear    = 1961
Eyear    = 1971
SpinRange = Eyear-Syear
nSPINs_1 = 7
years_1  = np.arange(0,SpinRange*nSPINs_1)

nSPINs_2 = 3
years_2  = np.arange(0,SpinRange*nSPINs_2)+years_1[-1]

print K_DATA_DIR+K_FILE_TAG_1+'1.'+str(Syear)+'0101.0.nc'
inf=nc.Dataset(K_DATA_DIR+K_FILE_TAG_1+'1.'+str(Syear)+'0101.0.nc','r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]
inf.close()
index_names = ['Global','Northern Lats','Tropical','Southern Lats']
indexes     = [ lats>-9999., lats>30 , (lats<15)&(lats>-15), lats<-30 ]
colours     = [ 'red' , 'blue', 'orange', 'darkgreen' ]


K_CS_1 = {}
K_CS_2 = {}
for zone in index_names:
    K_CS_1[zone]= [] 
    K_CS_2[zone]= []

print 'KovenSpin-1'
for iSPIN in range(nSPINs_1):
    #print iSPIN 
    for year in range(Syear,Eyear):
        #print str(year)
        file=K_DATA_DIR+K_FILE_TAG_1+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        if os.path.isfile(file):
            inf=nc.Dataset(file,'r')
            cs = inf.variables['cs'][:]
            cs = np.sum(cs,axis=0)
            for zone,index in zip(index_names,indexes):
                K_CS_1[zone].append(np.mean(cs[index]))
            inf.close()
        else:
            for zone in index_names:
                K_CS_1[zone].append(np.nan)

print 'KovenSpin-2'
for iSPIN in range(nSPINs_2):
    #print iSPIN 
    for year in range(Syear,Eyear):
        #print str(year)
        file=K_DATA_DIR+K_FILE_TAG_2+str(iSPIN+1)+'.'+str(year)+'0101.0.nc'
        if os.path.isfile(file):
            inf=nc.Dataset(file,'r')
            cs = inf.variables['cs'][:]
            cs = np.sum(cs,axis=0)
            for zone,index in zip(index_names,indexes):
                K_CS_2[zone].append(np.mean(cs[index]))
            inf.close()
        else:
            for zone in index_names:
                K_CS_2[zone].append(np.nan)


SC_CS_1 = {}
SC_CS_2 = {}
SC_CS_3 = {}
for zone in index_names:
    SC_CS_1[zone]= [] 
    SC_CS_2[zone]= [] 
    SC_CS_3[zone]= [] 

# Read in the cs from the SC runs
# spin 1
inf=nc.Dataset(SC_DATA_DIR+SC_FILE_TAG+i+'.monthly_mean.nc','r')
cs=inf.variables['cs'][:].squeeze()
cs=np.sum(cs,axis=1)
for zone,index in zip(index_names,indexes):
    SC_CS_1[zone]=np.mean(cs[:,index],axis=1)
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

   

    
