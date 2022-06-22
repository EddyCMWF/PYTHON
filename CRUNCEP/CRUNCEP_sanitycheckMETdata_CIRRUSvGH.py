#!/usr/bin/python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import netcdftime as nctime


GH_DIR='/users/eow/edwcom/CRUNCEP/garr_data/'
garr_color = 'red'
garr_name = 'garr_data'

CIRRUS_DIR='/users/eow/edwcom/CRUNCEP/cirrus_data/cruncep_v3_red'
cirrus_color  = 'blue'
cirrus_name = 'cirrus_data'

diff_color = 'black'
diff_name = 'GH_minus_Cirrus'

params=['LWdown','SWdown','Precip','PSurf','Qair','Tair']  #,'Wind']

F=plt.figure()

date='201001'

pltcnt=0
for param in params:
    AX=F.add_subplot(6,2,(pltcnt*2)+1)
    
    inf_c4  = nc.Dataset(GH_DIR+'/cruncep_'+param+'_'+date+'.nc')
    data_c4 = inf_c4.variables[param][:]
    mean_c4 = np.mean(data_c4,axis=1)
    std_c4  = np.std(data_c4,axis=1)
    time_c4 = nctime.num2date(inf_c4.variables['time'][:],inf_c4.variables['time'].units)
    inf_c4.close()
    
    line1=AX.plot(time_c4,mean_c4,color=garr_color,ls='-',linewidth=1.0,label=garr_name)
    AX.plot(time_c4,mean_c4+std_c4,color=garr_color,ls='--',linewidth=0.5)
    AX.plot(time_c4,mean_c4-std_c4,color=garr_color,ls='--',linewidth=0.5)
    
    inf_c3  = nc.Dataset(CIRRUS_DIR+'/cruncep_'+param+'_'+date+'.nc')
    data_c3 = inf_c3.variables[param][:]
    mean_c3 = np.mean(data_c3,axis=1)
    std_c3  = np.std(data_c3,axis=1)
    time_c3 = nctime.num2date(inf_c3.variables['time'][:],inf_c3.variables['time'].units)
    inf_c3.close()
    
    line2=AX.plot(time_c3,mean_c3,color=cirrus_color,ls='-',linewidth=1.0,label=cirrus_name)
    AX.plot(time_c3,mean_c3+std_c3,color=cirrus_color,ls='--',linewidth=0.5)
    AX.plot(time_c3,mean_c3-std_c3,color=cirrus_color,ls='--',linewidth=0.5)
    AX.set_title(param)
    
    mean_diff=np.mean(data_c4-data_c3,axis=1)
    std_diff=np.std(data_c4-data_c3,axis=1)
    
    AX=F.add_subplot(6,2,(pltcnt*2)+2)
    line_diff=AX.plot(time_c3,mean_diff,color=diff_color,ls='-',linewidth=1.0,label=diff_name)
    AX.plot(time_c3,mean_diff+std_diff,color=diff_color,ls='--',linewidth=0.5)
    AX.plot(time_c3,mean_diff-std_diff,color=diff_color,ls='--',linewidth=0.5)
    AX.set_title(param)
    
   
    pltcnt+=1



#F.legend((line1,line2,line3),('cruncep_v4','cruncep_v3','NCEP_v4'),loc=(0.6,0.1))


plt.show()
