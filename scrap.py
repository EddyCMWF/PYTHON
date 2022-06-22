
PTs=reload(PTs)



fig=plt.figure(figsize=(18,12))

for i in range(4):
    ax=fig.add_subplot(2,2,i+1)
    PTs.plot_map(transcom_regions_2D,lons_2D,lats_2D,AXIS=ax, \
                 CBAR_SIZE=str((i*4)+4)+'%',CBAR_PAD=0.2,RESOLUTION='l')

plt.show()


import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

file1='Jv4.5_WFDEI_nti_TRIFFID_gridded_threehourly_co2flux_respFix.2014.nc'
inf1=nc.Dataset(file1,'r')

file2='Jv4.5_WFDEI_nti_TRIFFID_gridded_threehourly_co2flux.2014.nc'
inf2=nc.Dataset(file2,'r')

data_inname='resp_s_gb'
dataname='Soil-Respiration'

data1=inf1.variables[data_inname][-100:,:]*1e3*3600*24
data2=inf2.variables[data_inname][-100:,:]*1e3*3600*24

plt.plot(np.mean(data1.reshape(data1.shape[0],-1),axis=1),label='fixed')
plt.plot(np.mean(data2.reshape(data2.shape[0],-1),axis=1),label='original')
plt.show()


data1=np.mean(inf1.variables[data_inname][-100:,:],axis=0)*1e3*3600*24
data2=np.mean(inf2.variables[data_inname][-100:,:],axis=0)*1e3*3600*24
data_diff=data1-data2
data_mean=(data1+data2)/2.

global_max=10
uk_max=5

plt.figure(figsize=(20,10))
plt.subplot(2,2,1)
plt.imshow(data1,origin='bottom',vmax=global_max)
plt.title('old')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,2)
plt.imshow(data2,origin='bottom',vmax=global_max)
plt.title('new')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,3)
plt.imshow(data_diff,origin='bottom',vmin=0-(global_max/10.),vmax=(global_max/10.),cmap='RdBu_r')
plt.title('old - new')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,4)
plt.imshow(data1/data2,origin='bottom',vmin=0,vmax=2,cmap='RdBu_r')
plt.title('old / new')
plt.colorbar().set_label('unitless ratio')
plt.suptitle(dataname,fontsize=30)
plt.savefig('/home/users/ecomynplatt/Soil_Resp_Corr_Plots/Global_'+dataname+'_comparison_SoilRespFix.png')
plt.show()

xlims=[330,375]
ylims=[275,305]
plt.figure(figsize=(15,10))
plt.subplot(2,2,1)
plt.imshow(data1[ylims[0]:ylims[1],xlims[0]:xlims[1]],origin='bottom',vmax=uk_max)
plt.title('old')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,2)
plt.imshow(data2[ylims[0]:ylims[1],xlims[0]:xlims[1]],origin='bottom',vmax=uk_max)
plt.title('new')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,3)
plt.imshow(data_diff[ylims[0]:ylims[1],xlims[0]:xlims[1]],origin='bottom',vmin=0-(uk_max/10),vmax=uk_max/10,cmap='RdBu_r')
plt.title('old - new')
plt.colorbar().set_label('gC m$^{-2}$ per day')
plt.subplot(2,2,4)
plt.imshow((data1/data2)[ylims[0]:ylims[1],xlims[0]:xlims[1]],origin='bottom',vmin=0.8,vmax=1.2,cmap='RdBu_r')
plt.title('old / new')
plt.colorbar().set_label('unitless ratio')
plt.suptitle(dataname,fontsize=30)
plt.savefig('/home/users/ecomynplatt/Soil_Resp_Corr_Plots/UK_'+dataname+'_comparison_SoilRespFix.png')
plt.show()



