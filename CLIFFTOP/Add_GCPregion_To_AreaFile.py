#!/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats


ancil_dir='/prj/CLIFFTOP/COMMON_DATA/ANCILS/'
area_file_1d=ancil_dir+'Area_in_iris_format.nc'
area_file_2d=ancil_dir+'grid_info.nc'

GCP_dir='/prj/CLIFFTOP/'
GCP_file=GCP_dir+'GCP_methane_regions_1x1_ext.nc'

# In[67]:

#Read IMAGE data:
GCPinf = nc.Dataset(GCP_file,'r')
GCPlons = GCPinf.variables['LON'][:]
GCPlats = GCPinf.variables['LAT'][:]
GCPdata = GCPinf.variables['REG'][:]
GCPinf.close()

lons_2d, lats_2d = np.meshgrid(GCPlons, GCPlats)

GCP_data = np.ma.masked_equal(GCPdata,0)
#GCP_data = GCP_data[::-1,:]

plt.figure(figsize=[15,11])
plt.imshow(GCP_data) # , origin='bottom')
plt.colorbar()
plt.show()


#read in the 1d imogen data:
inf=nc.Dataset(area_file_1d,'a')
im_lat=inf.variables['latitude'][:].squeeze()
im_lon=inf.variables['longitude'][:].squeeze()



lat_res=2.5
lon_res=3.75
n_impoints=len(im_lat)
fill_value=-999.
im_GCP_1d = np.zeros_like(im_lat)

for ipt in range(n_impoints):
    lat,lon = im_lat[ipt],im_lon[ipt]
    
    cell_mask = (  (lats_2d>=lat)&(lats_2d<lat+lat_res)
                 & (lons_2d>=lon)&(lons_2d<lon+lon_res) 
                 & (GCP_data.mask==False)
                 )
    
    
    if len(GCP_data[cell_mask])>=1:
        print(ipt,lat,lon,stats.mode(GCP_data[cell_mask]))
        im_GCP_1d[ipt] = stats.mode(GCP_data[cell_mask])[0]
    else:
        print(ipt,lat,lon,fill_value)
        im_GCP_1d[ipt] = fill_value
    
    
badex=np.where(im_GCP_1d==0)[0]


print(badex)
print(im_lat[badex])
print(im_lon[badex])

# fix gap fills manually:

GCP_region_names="0:Oceans,1:Boreal North America,2:USA contiguous,3:Central North America,4:Tropical South America,5:Temperate South America,6:Northern Africa,7:Southern Africa,8:Russia,9:Oceania,10:Europe,11:China, 12:India,13:South East Asia,14:Central Eurasia & Japan"

#add transcomm to Area file
outvar=inf.createVariable('GCP_region','int',('time','y','x'))
#outvar=inf.variables['GCP_region']
outvar.long_name='GCP index'
outvar.region_names=GCP_region_names
outvar.index_values=list(range(0,15))


outvar[:] = im_GCP_1d


# In[124]:

inf.close()


# In[125]:

# Read in 2d grid info to check and store in 2d
inf_2d=nc.Dataset(area_file_2d,'a')
lat_2d = inf_2d.variables['latitude'][:]
lon_2d = inf_2d.variables['longitude'][:]
land_index = inf_2d.variables['land_index'][:]


# In[130]:

im_GCP_2d=np.ma.masked_array(im_GCP_1d[land_index],mask=land_index.mask,fill_value=fill_value)
im_GCP_2d.data[im_GCP_2d.mask==True]=fill_value
plt.figure(figsize=[20,15])
plt.imshow(im_GCP_2d[:,:],origin='bottom') #,vmin=-5)
plt.colorbar()
plt.show()

# In[133]:

#add transcomm to Area file
outvar=inf_2d.createVariable('GCP_region','float32',('y','x'),fill_value=fill_value)
outvar.long_name='GCP index'
outvar[:] = im_GCP_2d


# In[134]:

inf_2d.close()


# In[ ]:



