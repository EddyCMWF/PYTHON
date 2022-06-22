#!/usr/bin/python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import netcdftime as nctime


dir_cruncep_v4='/prj/ALANIS/jules_data/cru_ncep/0.5deg/meteo/v4/6hr/'
c_4_color = 'red'
c_4_name = 'cruncep_v4'

dir_cruncep_v3='/prj/ALANIS/jules_data/cru_ncep/0.5deg/meteo/v3.1/6hr/data/'
c_3_color  = 'blue'
c_3_name = 'cruncep_v3'

dir_NCEP_v4='/prj/wetlands_africa/NCEP_CRU_MET/v4_tstep/'
N_4_color  = 'green'
N_4_name = 'NCEP_CRU_v4'


params=['LWdown','Precip','PSurf','Qair','SWdown','Tair','Wind']

pltcnt=1

F=plt.figure()

date='200001'

for param in params:
    AX=F.add_subplot(4,2,pltcnt)
    
    inf_c4  = nc.Dataset(dir_cruncep_v4+param+'/cruncep_'+param+'_'+date+'.nc')
    data_c4 = inf_c4.variables[param][:]
    mean_c4 = np.mean(data_c4,axis=1)
    std_c4  = np.std(data_c4,axis=1)
    time_c4 = nctime.num2date(inf_c4.variables['time'][:],inf_c4.variables['time'].units)
    inf_c4.close()
    
    line1=AX.plot(time_c4,mean_c4,color=c_4_color,ls='-',linewidth=1.0,label='cruncep_v4')
    AX.plot(time_c4,mean_c4+std_c4,color=c_4_color,ls='--',linewidth=0.5)
    AX.plot(time_c4,mean_c4-std_c4,color=c_4_color,ls='--',linewidth=0.5)
    
    inf_c3  = nc.Dataset(dir_cruncep_v3+param+'/cruncep_'+param+'_'+date+'.nc')
    data_c3 = inf_c3.variables[param][:]
    mean_c3 = np.mean(data_c3,axis=1)
    std_c3  = np.std(data_c3,axis=1)
    time_c3 = nctime.num2date(inf_c3.variables['time'][:],inf_c3.variables['time'].units)
    inf_c3.close()
    
    line2=AX.plot(time_c3,mean_c3,color=c_3_color,ls='-',linewidth=1.0,label='cruncep_v3')
    AX.plot(time_c3,mean_c3+std_c3,color=c_3_color,ls='--',linewidth=0.5)
    AX.plot(time_c3,mean_c3-std_c3,color=c_3_color,ls='--',linewidth=0.5)
        
    #inf_N4  = nc.Dataset(dir_NCEP_v4+param+'/cruncep_'+param+'_'+date+'.nc')
    #data_N4 = inf_N4.variables[param][:]
    #mean_N4 = np.mean(data_N4,axis=1)
    #std_N4  = np.std(data_N4,axis=1)
    #time_N4 = nctime.num2date(inf_N4.variables['tstep'][:],inf_N4.variables['tstep'].units)
    #inf_N4.close()
    
    #line3=AX.plot(time_N4,mean_N4,color=N_4_color,ls='-',linewidth=1.0,label='NCEP_v4')
    #AX.plot(time_N4,mean_N4+std_N4,color=N_4_color,ls='--',linewidth=0.5)
    #AX.plot(time_N4,mean_N4-std_N4,color=N_4_color,ls='--',linewidth=0.5)
    
    AX.set_title(param)
    
    pltcnt+=1



#F.legend((line1,line2,line3),('cruncep_v4','cruncep_v3','NCEP_v4'),loc=(0.6,0.1))


plt.show()
