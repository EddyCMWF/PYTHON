import pylab as plt
import numpy as np
import netCDF4 as nc

params=['LWdown', 'Precip', 'PSurf', 'Qair', 'SWdown', 'Tair', 'Wind' ]

data_1981=[]
for param in params:
  file=param+'/cruncep_'+param+'_198107.nc'
  inf=nc.Dataset(file,'r')
  data_1981.append(inf.variables[param][:])
  inf.close()


for i in range(len(params)):
  plt.subplot(4,2,i+1)
  dmean=np.mean(data_1981[i],axis=1)
  plt.plot(dmean)
  dmax =np.amax(data_1981[i],axis=1)
  plt.plot(dmax)
  dmin =np.amin(data_1981[i],axis=1)
  plt.plot(dmin)
  plt.title=params[i]


plt.show()
