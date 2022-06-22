#!/bin/env python

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


data_dir='/users/eow/edwcom/WFD_EI/' 

infile=data_dir+'JULES_output/Jv4.5_WFDEI_nti_TRIFFID_nounf.monthly_soil.2000.nc'

gridfile=data_dir+'wfdei-land-mask.nc'
grinf=nc.Dataset(gridfile,'r')
grindex=grinf.variables['land_index'][:]
grlats=grinf.variables['latitude'][:]
grlons=grinf.variables['longitude'][:]
grinf.close()

land_file=data_dir+'qrparm.veg.fracNew.nc'
lfinf=nc.Dataset(land_file,'r')
lf=lfinf.variables['field1391'][:]
lfinf.close()
ice_mask= lf[-1,:]==1

scpool_resp_rates=np.array([3.22e-7,9.65e-9,2.12e-8,6.43e-10,])

inf=nc.Dataset(infile,'r')
cs_raw=inf.variables['cs'][:].squeeze()
inf.close()

ice_mask_raw=np.array([ [ice_mask for i in range(4)] for j in range(12)])
cs=np.ma.masked_array(cs_raw,mask=ice_mask_raw)

cs_total=np.sum(cs,axis=1)

cs_by_resp_total=np.sum(cs.transpose(0,2,1)*scpool_resp_rates,axis=2)

single_resp_rate=cs_by_resp_total/cs_total

print('Mean','Std. Dev.')
print(np.mean(single_resp_rate),np.std(single_resp_rate))

mean_single_RR=np.mean(single_resp_rate[:,grindex-1],axis=0)
mean_single_RR.mask= (mean_single_RR.mask==True)|(grindex.mask==True)

std_single_RR=np.std(single_resp_rate[:,grindex-1],axis=0)
std_single_RR.mask= (std_single_RR.mask==True)|(grindex.mask==True)

#fig,axes=plt.subplots(nrows=2,ncols=1,figsize=[15,15])
#plt.subplot(2,1,1)
plt.imshow(mean_single_RR,origin='bottom')
plt.colorbar()

#plt.subplot(2,1,2)
#plt.imshow(std_single_RR,origin='bottom')
#plt.colorbar()
plt.show()

