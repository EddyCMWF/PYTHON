#!/usr/bin/env python
#
# 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

DATA_DIR   ='/prj/GREENHOUSE/SC_simulator/JULES_output/WFDEI_nocomp/'
FILE_TAG   ='JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin_nocomp.dump.spin'
OUTPUT_DIR = '/users/eow/edwcom/SC_simulator/plots/'

nSPINs=100
years=np.arange(1,1001,10)
iDISPLAY='N'

nPOOLs=4
fill_value=-999.
CS_pools = ['DPM', 'RPM', 'BIO', 'HUM']
CS_units = '$kg$C $m^{-2}$'
CS_maxes = [0.1,5,0.4,12]

inf=nc.Dataset(DATA_DIR+FILE_TAG+'1.19800101.0.nc','r')
lats=inf.variables['latitude'][:]
lons=inf.variables['longitude'][:]
inf.close()

index_names = ['Global','Northern Lats','Tropical','Southern Lats']
indexes     = [ lats>-9999., lats>30 , (lats<15)&(lats>-15), lats<-30 ]
colours     = [ 'red' , 'blue', 'orange', 'darkgreen' ]

CS = {}
for zone in index_names:
    CS[zone]=[ [] for i in range(nPOOLs)]


for iSPIN in range(nSPINs):
    print iSPIN
    file=DATA_DIR+FILE_TAG+str(iSPIN+1)+'.19800101.0.nc'
    inf=nc.Dataset(file,'r')
    cs = inf.variables['cs'][:]
    for zone,index in zip(index_names,indexes):
        for iPOOL in range(nPOOLs):
            CS[zone][iPOOL].append(np.mean(cs[iPOOL,index]))

    inf.close()



# Create figure to plot to
FIG = plt.figure(figsize=[18,12])
nplts_wdth   = 2.
nplts_hght = np.ceil(len(CS_pools)/nplts_wdth)


for iPOOL in range(nPOOLs):
    #print CS_pools[iPOOL]
    AX  = FIG.add_subplot(nplts_hght,nplts_wdth,iPOOL+1)
    AX.set_title(CS_pools[iPOOL])

    for iZON in range(len(indexes)):
        index=indexes[iZON]
        zone=index_names[iZON]
        # print index.shape
        CS_timeseries = np.array(CS[zone][iPOOL])
        
        AX.plot(years,CS_timeseries,label=zone,color=colours[iZON],lw=2)
        

    AX.set_ybound(lower=0,upper=12)

AX.legend( bbox_to_anchor=(-0.1,-0.11),loc=10,borderaxespad=0.,ncol=4 )
FIG.text(0.03,0.6,'Soil Carbon $kg$ C $m^{-2}$',rotation='vertical',fontsize=24)
FIG.tight_layout(rect=[0.05,0.05,1.0,0.96])
FIG.suptitle('Standard Soil Carbon Spin-Up Time-series',fontsize=24)

if (iDISPLAY=='Y'):
    plt.show()
else:
    FIG.savefig(OUTPUT_DIR+'BigSpinUp_Full_Timeseries.png')
    plt.close()

   

    
