#/usr/bin/python2.7

import netCDF4 as nc
import numpy as np
import pylab as plt

soil_file='/users/eow/edwcom/EMEP/hwsd2emep/EMEP4UK_soilparams_hwsd_vg.nc'

LandFrac_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_LandFrac.nc.noIAM.oneSOIL'

LF_inf=nc.Dataset(LandFrac_file,'r')
LSmask=LF_inf.variables['LSmask'][:]
LF_inf.close()

soil_inf=nc.Dataset(soil_file,'r')
sathh=soil_inf.variables['oneoveralpha'][:]


temp=sathh[0,:,:].copy()
temp.mask=[LSmask==0]

plt.subplot(1,2,1)
plt.imshow(temp,origin='bottom')
plt.colorbar()
plt.subplot(1,2,2)
plt.imshow(sathh[0,:,:],origin='bottom')
plt.colorbar()
plt.show()


/users/eow/edwcom/EMEP/hwsd2emep/input_data/grid_files/EMEP4UK_gridfile.nc

LSMfile='/users/eow/edwcom/EMEP/EMEP4UK/LANDMASK.d3'
LSM_data       = np.loadtxt(LSMfile)
LSM_data=LSM_data[:,2].reshape(270,220)/100.

plt.imshow(LSM_data-LSmask,origin='bottom')

plt.colorbar()
plt.show()


